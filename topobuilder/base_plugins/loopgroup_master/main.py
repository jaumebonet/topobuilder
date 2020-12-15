# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from pathlib import Path
from typing import Dict, Optional, Union, List
import math
from ast import literal_eval
from subprocess import run, DEVNULL
import gzip
import itertools

# External Libraries
import pandas as pd
from pandas.compat import StringIO
from rstoolbox.io import parse_rosetta_fragments, write_rosetta_fragments, parse_master_file

# This Library
from topobuilder.case import Case
from topobuilder.workflow import Node, NodeDataError
from .core import core
import topobuilder.utils as TButil
import topobuilder.core as TBcore
import topobuilder.utils.plot as TBPlot
import topobuilder.utils.master as TBMaster
import topobuilder.utils.structure as TBstructure


__all__ = ['loopgroup_master']


class loopgroup_master( Node ):
    """Find loops capable of closing a super-secondary structure with :term:`MASTER` searches.

    The search will provide loops for each pair of SSE, together with fragment sets covering each
    super-secondary structure.

    .. note::
        Depends on the ``master.pds`` configuration option.
        Depends on the ``loop_master.abego`` configuration option.
        Depends on the ``loop_master.fragments`` configuration option.

    .. caution::
        In order to execute this :class:`.Node`, `MASTER <https://grigoryanlab.org/master/>`_ needs
        to be installed.

    .. admonition:: To Developers

        Due to its use in multiple :class:`.Node`, functions to deal with master are mostly located
        in the :mod:`.utils` module.

    :param loop_range: Expected loop length is calculated from the euclidian distance between two secondary
        structures. This attribute adds a window of ``loop_range`` residues under and over the calculated
        length.
    :param top_loops: Number of loops that are selected from the sorted matches to retrieve final loop candidates.
    :param hairpins_2: When :data:`True`, enforce 2-residue loops on beta hairpins.
    :param rmsd_cut: Threshold value to include loops as match candidates.
    :param filter: List of PDB identifiers to use for the search. If nothing is provided, the full database,
        as defined by the ``master.pds`` global option, will be used.

    :raises:
        :NodeDataError: On **initialization**. If the PDS database, the ABEGO or fragments cannot be found.
        :NodeDataError: On **check**. If the required fields to be executed are not there.
        :NodeDataError: On **execution**. If the :class:`.Case` contains anything other than one defined connectivity.

    """
    REQUIRED_FIELDS = ('topology.architecture', 'topology.connectivity')
    RETURNED_FIELDS = ('metadata.loop_fragments', 'metadata.loop_lengths')
    VERSION = 'v1.0'

    def __init__( self, tag: int,
                  loop_groups: List,
                  loop_range: Optional[int] = 3,
                  top_loops: Optional[int] = 20,
                  #hairpins_2: Optional[bool] = True,
                  rmsd_cut: Optional[float] = 5.0 ):
        super(loopgroup_master, self).__init__(tag)

        self.steps = loop_groups
        self.loop_range = loop_range
        self.top_loops = top_loops
        #self.hairpins_2 = hairpins_2
        self.rmsd_cut = rmsd_cut
        self.multiplier = 3 # applied to the top_loops selection to give the file some margin.

        # Force is set to True as another MASTER-dependent pluggin might
        # work with a different filter of the database.
        self.pdsdb, _ = TBMaster.pds_database(self.log)
        #self.pdsdb, _ = TBMaster.pds_database(self.log, self.filter, force=True)

        # Load ABEGO and FRAGMENT data specific of this plugin
        self.fragments = self.get_fragfiles()
        self.abegos = self.get_abegos()

    def single_check( self, dummy: Dict ) -> Dict:
        kase = Case(dummy)

        # Check what it needs
        for itag in self.REQUIRED_FIELDS:
            if kase[itag] is None:
                raise NodeDataError(f'Field "{itag}" is required')

        # Include what keywords it adds (in this instance, nothing)
        kase.data.setdefault('metadata', {}).setdefault('loop_fragments', [])
        kase.data.setdefault('metadata', {}).setdefault('loop_lengths', [])
        return kase.data

    def single_execute( self, data: Dict ) -> Dict:
        kase = Case(data)
        # Loop MASTER is only applied to a Case with one single connectivity
        if kase.connectivity_count != 1:
            err = f'{self.nodeID} can only be applied to one connectivity. '
            err += f'Current case contains a total of {kase.connectivity_count}.'
            raise NodeDataError(err)
        # And has to be reoriented
        if not kase.is_reoriented:
            self.log.debug('Topology was provided without oriented SSE -> orienting.')
            kase = kase.apply_topologies()[0]

        # Generate the folder tree for a single connectivity.
        folders = kase.connectivities_paths[0].joinpath('loopgroup_master')
        folders.mkdir(parents=True, exist_ok=True)

        # Global step distance
        loop_step = kase.cast_absolute()['configuration.defaults.distance.loop_step']

        self.log.debug('STEPS HERE')
        self.log.debug(self.steps)

        # Output keys
        kase.data.setdefault('metadata', {}).setdefault('loop_fragments', [])
        kase.data.setdefault('metadata', {}).setdefault('loop_lengths', [])

        # Find steps: Each pair of secondary structure.
        #it = kase.connectivities_str[0].split('.')
        #steps = [it[i:i + 2] for i in range(0, len(it) - 1)]
        lengths = kase.connectivity_len[0]
        start = 1

        for i, (group, infos) in enumerate(self.steps.items()):
            self.log.info(f'Search at group {group}')

            # 1. Make folders and files
            wfolder = folders.joinpath(f'loopgroup{i + 1:02d}')
            wfolder.mkdir(parents=True, exist_ok=True)
            outfile = wfolder.joinpath(f'loopgroup_master.iter{i + 1:02d}.pdb')
            outfilePDS = wfolder.joinpath(f'loopgroup_master.iter{i + 1:02d}.pds')
            masfile = outfile.with_suffix('.master')
            checkpoint = wfolder.joinpath('checkpoint.json')

            # 2. Check if checkpoint exists, retrieve and skip
            reload = TButil.checkpoint_in(self.log, checkpoint)
            if reload is not None:
                kase.data['metadata']['loop_fragments'].append(reload)
                kase.data['metadata']['loop_lengths'].append(int(reload['edges']['loop']))
                start += (int(reload['edges']['sse1']) + int(reload['edges']['loop']))
                continue

            # 3. Check hairpin
            # Get SSEs and identifiers
            sses = [kase.get_sse_by_id(sse) for sse in infos[0]]
            #sse1_name, sse2_name = sse1['id'], sse2['id']
            #is_hairpin = self.check_hairpin(sse1_name, sse2_name)

            # 4. Generate structures
            sses = TBstructure.build_pdb_object(self.log, sses, 5, concat=False, outfile=outfile)

            if not masfile.is_file():
                # 5. calculate expected loop length by loop_step
                #Mdis, mdis = TBstructure.get_loop_length(self.log, sse1, sse2, loop_step, self.loop_range)

                # 6. Run MASTER
                #outfilePDS = outfile if outfile is not None else Path(outfile).with_suffix('.pds')
                self.log.debug(f'FILE {outfilePDS}')
                # -> make PDS query
                cmd = TBMaster.createPDS(outfile, outfilePDS)
                self.log.debug(f'EXECUTE: {" ".join(cmd)}')
                run(cmd, stdout=DEVNULL)
                # -> run MASTER
                cmd = TBMaster.master_groupedgap(outfilePDS, self.pdsdb, masfile, infos[1], self.rmsd_cut)
                self.log.debug(f'EXECUTE: {" ".join(cmd)}')
                result = run(cmd, stdout=DEVNULL)

                # TODO: implement motif compability
                # if result.returncode: # no loop between that connection, e.g. a motif ranging over multiple sse with keeping the loops
                #     # 4. Generate structures
                #     self.log.debug('generate combined structure')
                #     sse = pd.concat([sse1, sse2], sort=False)
                #
                #     # 6. Run MASTER
                #     self.log.debug(Path(outfile))
                #     #outfilePDS = outfile if outfile is not None else Path(outfile).with_suffix('.pds')
                #     self.log.debug(f'FILE {outfilePDS}')
                #     # -> make PDS query
                #     cmd = TBMaster.createPDS(outfile, outfilePDS)
                #     self.log.debug(f'EXECUTE: {" ".join(cmd)}')
                #     run(cmd, stdout=DEVNULL)
                #     # -> run MASTER
                #     cmd = TBMaster.master_nogap(outfilePDS, self.pdsdb, masfile, self.rmsd_cut)
                #     self.log.debug(f'EXECUTE: {" ".join(cmd)}')
                #     run(cmd, stdout=DEVNULL)
                #
                #     # 6. Minimize master data (pick top_loopsx3 lines to read and minimize the files)
                #     match_count = self.minimize_master_file(masfile)
                #     self.log.debug(f'match count here {match_count}')
                #
                #     # 7. Retrieve MASTER data
                #     dfloop = self.process_master_data_no_gap(masfile, sse1_name, sse2_name)
                #     sse1l, loopl, sse2l = lengths[i], int(dfloop['loop_length'].values[0]), lengths[i + 1]
                #     total_len = sse1l + loopl + sse2l
                #     end_edge = total_len + start - 1
                #     edges = {'ini': int(start), 'end': int(end_edge), 'sse1': int(sse1l), 'loop': int(loopl), 'sse2': int(sse2l)}
                #     self.log.debug(f'INI: {start}; END: {end_edge}; SSE1: {sse1l}; LOOP: {loopl}; SSE2: {sse2l}')
                #     self.log.debug(dfloop.to_string())
                #
                #     # 8. Bring and Combine fragments from the different sources.
                #     loop_data = self.make_fragment_files(dfloop, edges, masfile)
                #     loop_data['match_count'] += match_count

                #else:

                # 6. Minimize master data (pick top_loopsx3 lines to read and minimize the files)
                match_count = self.minimize_master_file(masfile)
                # 7. Retrieve MASTER data
                df_container = self.process_master_data(masfile, infos[0], infos[1], infos[2])

                loop_datas = []
                for indx in list(df_container.order.drop_duplicates()):
                    dfloop = df_container[df_container.order == indx]
                    sse1l, loopl, sse2l = lengths[i], int(dfloop['loop_length'].values[0]), lengths[i + 1]
                    total_len = sse1l + loopl + sse2l
                    end_edge = total_len + start - 1
                    edges = {'ini': int(start), 'end': int(end_edge), 'sse1': int(sse1l), 'loop': int(loopl), 'sse2': int(sse2l)}
                    self.log.debug(f'INI: {start}; END: {end_edge}; SSE1: {sse1l}; LOOP: {loopl}; SSE2: {sse2l}')
                    self.log.debug(dfloop.to_string())

                    # 8. Bring and Combine fragments from the different sources.
                    loop_data, nfolder = self.make_fragment_files(dfloop, edges, masfile, wfolder, no_loop=True)
                    loop_data['match_count'] += match_count

                    # 9. Save data in the Case
                    kase.data['metadata']['loop_fragments'].append(loop_data)
                    kase.data['metadata']['loop_lengths'].append(int(loopl))
                    start += (sse1l + loopl)

                    # 10. Checkpoint save
                    #checkpoint = nfolder.joinpath('checkpoint.json')
                    loop_datas.append(loop_data)

                TButil.checkpoint_out(self.log, checkpoint, loop_datas)

        return kase

    def get_fragfiles( self ) -> pd.DataFrame:
        """Obtain the fragment files.
        """
        fragpath = Path(core.get_option('loop_master', 'fragments'))
        self.log.debug(f'Listing available fragment files at: {fragpath.name}')
        if not fragpath.is_dir():
            raise NodeDataError(f'{fragpath.name} is not a folder.')
        return pd.DataFrame([(x.name[:4], x.name[5:6], x, y) for x, y in zip(sorted(fragpath.glob('*/*3mers.gz')),
                                                                             sorted(fragpath.glob('*/*9mers.gz')))],
                              columns=['pdb', 'chain', '3mers', '9mers'])

    def get_abegos( self ) -> pd.DataFrame:
        """Load ABEGO data.
        """
        abegos = Path(core.get_option('loop_master', 'abego'))
        if not abegos.is_file():
            raise NodeDataError(f'ABEGO file {abegos.name} cannot be found.')

        self.log.debug(f'Loading ABEGO data from: {abegos.name}\n')
        doopen = gzip.open if abegos.suffix == '.gz' else open
        abegodata = []
        with doopen(abegos, 'rt') as fd:
            for line1, line2 in itertools.zip_longest(*[fd] * 2):
                line2 = line2 if len(line2.strip()) != 0 else 'NON\n'
                line1 = line1.strip().lstrip('>').split('_')
                abegodata.append(f'{line1[0]},{line1[1]},{line2}')
        abegodata = pd.read_csv(StringIO(''.join(abegodata)), names=['pdb', 'chain', 'abego'], header=None)
        abegodata = abegodata[abegodata['abego'] != 'NON']
        return abegodata

    def minimize_master_file( self, masfile: Path ) -> int:
        """Pick a limited number of matches.

        Changes the original file to only the selected number.
        """
        try:
            with open(masfile) as fd:
                num_lines = sum(1 for line in fd if line.rstrip())
            with open(masfile) as fd:
                head = [next(fd) for x in range(self.top_loops * self.multiplier)]
            with open(masfile, 'w') as fd:
                fd.write(''.join(head))
        except StopIteration:
            pass
        return num_lines

    def process_master_data( self, masfile: Path, names: List, loop_lengths: str, loop_orders: str ) -> pd.DataFrame:
        """Get length data from the MASTER matches.
        """
        def cutter(row, num):
            match = row['match']
            # MASTER starts match count at 0!
            self.log.debug(match[num][1] + 1)
            self.log.debug(match[num + 1][0])
            self.log.debug(row['abego'])
            loop = row['abego'][match[num][1] + 1: match[num + 1][0]]
            return row['abego'][match[num][0]: match[num + 1][1] + 1], loop, len(loop), match[num][0], match[num + 1][1] + 1

        if masfile.with_suffix('.csv').is_file():
            df = pd.read_csv(masfile.with_suffix('.csv'))
            df['match'] = df['match'].apply(literal_eval)
            return df

        llens = loop_lengths.split(';')
        pnames = [(names[i],names[i+1]) for i in range(len(names) - 1)]
        lorder = loop_orders.split(';')

        dfloop = parse_master_file(masfile)
        dfloop = dfloop.merge(self.abegos, on=['pdb', 'chain']).merge(self.fragments, on=['pdb', 'chain']).dropna()

        container = []
        for k, (pname, llen, lord) in enumerate(zip(pnames, llens, lorder)):
            dfloop_copy = dfloop.copy()
            if lord is 'x': # skip regions that are not of interest
                continue
            self.log.info(f'Current jump is {pname[0], pname[1]} with {k}')
            dfloop_copy[['abego', 'loop', 'loop_length', 'start', 'stop']] = dfloop_copy.apply(cutter, num=k, axis=1, result_type='expand')
            dfloop_copy = dfloop_copy.iloc[:self.top_loops]
            self.log.debug('LOOP LENTS')
            self.log.debug(dfloop_copy['loop_length'])
            dfloop_copy['length_count'] = dfloop_copy.loop_length.map(dfloop_copy.loop_length.value_counts())
            self.log.debug('LOOP LENTS')
            self.log.debug(dfloop_copy['length_count'])
            dfloop_copy.drop(columns=['pds_path']).to_csv(masfile.with_suffix('.all.csv'), index=False)
            finaldf = dfloop_copy.sort_values('rmsd').drop_duplicates(['loop'])

            #pick = 0
            #if hairpin and 2 in finaldf['loop_length'].values:
            #    pick = 2
            #else:
            pick = finaldf[finaldf['length_count'] == finaldf['length_count'].max()]['loop_length'].min()
            finaldf = finaldf[finaldf['loop_length'] == pick]

            TBPlot.plot_loop_length_distribution(self.log, dfloop_copy, pick, masfile.with_suffix(''), f'loop {pname[0]} <-> {pname[1]}')

            df = finaldf.drop(columns=['pds_path'])
            df = df.assign(order=[int(lord)]*len(df))
            masfile2 = str(masfile) + f'.jump{int(lord):02d}.csv'
            df.to_csv(masfile2, index=False)
            container.append(df)
        return pd.concat(container).sort_values('order')

    # def process_master_data_no_gap( self, masfile: Path, name1: str, name2: str) -> pd.DataFrame:
    #     """Get length data from the MASTER matches.
    #     """
    #     def cutter(row):
    #         match = row['match']
    #         # MASTER starts match count at 0!
    #         return row['abego'][match[0][0]: match[0][-1] + 1], '-', 0
    #
    #     if masfile.with_suffix('.csv').is_file():
    #         df = pd.read_csv(masfile.with_suffix('.csv'))
    #         df['match'] = df['match'].apply(literal_eval)
    #         return df
    #
    #     dfloop = parse_master_file(masfile)
    #     dfloop = dfloop.merge(self.abegos, on=['pdb', 'chain']).merge(self.fragments, on=['pdb', 'chain']).dropna()
    #     dfloop[['abego', 'loop', 'loop_length']] = dfloop.apply(cutter, axis=1, result_type='expand')
    #     dfloop = dfloop.iloc[:self.top_loops]
    #     dfloop['length_count'] = dfloop.loop_length.map(dfloop.loop_length.value_counts())
    #     dfloop.drop(columns=['pds_path']).to_csv(masfile.with_suffix('.all.csv'), index=False)
    #     finaldf = dfloop.sort_values('rmsd').drop_duplicates(['loop'])
    #
    #     df = finaldf.drop(columns=['pds_path'])
    #     df.to_csv(masfile.with_suffix('.csv'), index=False)
    #     return df

    def make_fragment_files( self, dfloop: pd.DataFrame, edges: Dict, masfile: Path, wfolder: Path, no_loop: Optional[bool] = True ) -> Dict:
        """Combin the fragments from the different matches.
        """
        data = {'loop_length': int(dfloop.iloc[0]['loop_length']), 'abego': list(dfloop['loop'].values),
                'edges': edges, 'fragfiles': [], 'match_count': 0}

        dfs3 = []
        dfs9 = []
        sample = math.ceil(200 / dfloop.shape[0])
        if not no_loop:
            for i, row in dfloop.iterrows():
                # Remember: MASTER match starts with 0!
                dfs3.append((parse_rosetta_fragments(str(row['3mers']), source=f'{row["pdb"]}_{row["chain"]}')
                             .slice_region(row['start'], row['stop'] + 1).sample_top_neighbors(sample)
                             .renumber(edges['ini']).top_limit(edges['end'])))
                dfs9.append((parse_rosetta_fragments(str(row['9mers']), source=f'{row["pdb"]}_{row["chain"]}')
                             .slice_region(row['start'], row['stop'] + 1).sample_top_neighbors(sample)
                             .renumber(edges['ini']).top_limit(edges['end'])))
        else:
            for i, row in dfloop.iterrows():
                # Remember: MASTER match starts with 0!
                dfs3.append((parse_rosetta_fragments(str(row['3mers']), source=f'{row["pdb"]}_{row["chain"]}')
                             .slice_region(row['start'], row['stop'] + 1).sample_top_neighbors(sample)
                             .renumber(edges['ini']).top_limit(edges['end'])))
                dfs9.append((parse_rosetta_fragments(str(row['9mers']), source=f'{row["pdb"]}_{row["chain"]}')
                             .slice_region(row['start'], row['stop'] + 1).sample_top_neighbors(sample)
                             .renumber(edges['ini']).top_limit(edges['end'])))

        # Merge Fragments
        dfs3all = dfs3[0]
        dfs9all = dfs9[0]
        for i in range(1, len(dfs3)):
            dfs3all = dfs3all.add_fragments(dfs3[i], ini=edges['ini'], how='append')
            dfs9all = dfs9all.add_fragments(dfs9[i], ini=edges['ini'], how='append')
        dfs3all = dfs3all.sample_top_neighbors(200)
        dfs9all = dfs9all.sample_top_neighbors(200)

        # set up
        lord = int(dfloop.order.drop_duplicates().values[0])
        nfolder = masfile.parent.absolute().joinpath(f'loop{int(lord):02d}')
        self.log.debug(str(nfolder))
        nfolder.mkdir(parents=True, exist_ok=True)
        masfile2 = str(nfolder.joinpath(f'jump{int(lord):02d}'))

        self.log.debug(str(masfile2))

        self.log.debug('Writing 3mers fragfile\n')
        #data['fragfiles'].append(write_rosetta_fragments(dfs3all, prefix=str(masfile.with_suffix('')), strict=True))
        data['fragfiles'].append(write_rosetta_fragments(dfs3all, prefix=masfile2, strict=True))
        self.log.debug(f'3mers fragfile: {data["fragfiles"][-1]}\n')

        self.log.debug('Writing 9mers fragfile\n')
        #data['fragfiles'].append(write_rosetta_fragments(dfs9all, prefix=str(masfile.with_suffix('')), strict=True))
        data['fragfiles'].append(write_rosetta_fragments(dfs9all, prefix=masfile2, strict=True))
        self.log.debug(f'9mers fragfile: {data["fragfiles"][-1]}\n')

        dfs3all.drop(columns=['pdb', 'frame', 'neighbors', 'neighbor',
                              'aa', 'sse', 'phi', 'psi', 'omega']).to_csv(data['fragfiles'][0] + '.csv', index=False)
        dfs9all.drop(columns=['pdb', 'frame', 'neighbors', 'neighbor',
                              'aa', 'sse', 'phi', 'psi', 'omega']).to_csv(data['fragfiles'][1] + '.csv', index=False)
        imageprefix = masfile.with_suffix('.fragprofile')
        TBPlot.plot_fragment_templates(self.log, dfs3all, dfs9all, imageprefix)

        return data, nfolder

    # def check_hairpin( self, name1: str, name2: str) -> bool:
    #     """Check if a SSE-SSE pair is or not a hairpin.
    #     """
    #     if name1[0] != name2[0]:
    #         return False
    #     if name1[-1] != 'E':
    #         return False
    #     if int(name1[1]) == int(name2[1]) + 1:
    #         return True
    #     if int(name1[1]) == int(name2[1]) - 1:
    #         return True
    #     return False
