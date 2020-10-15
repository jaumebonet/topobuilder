# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Dict, Optional
import sys

# External Libraries

# This Library
from topobuilder.workflow import Node, NodeDataError
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from .parametric import SSEArchitect


__all__ = ['builder']


class builder( Node ):
    """Builds 3D SKETCH from given FORM string using ideal SSE elements.

    If corrections are available and specified, these will be applied onto the sketch.

    .. caution::
        In order to apply secondary structure or per layer corrections, the :mod:`.corrector` plugin
        needs to be set in the pipeline.

    :param connectivity: Expected secondary structure connectivity. *Important*: at the moment only a single
                         connectivity supported.
    :param pick_aa: Desired amino acid type to use for the SKETCH sequence. If not specified, it will
                    use pseudorandomly assign amino acid types based on secondary structure propensity scores.

    :raises:
        :NodeDataError: On **check**. If the required fields to be executed are not there.
        :NodeDataError: On **execution**. If the :class:`.Case` contains anything other than one defined connectivity.
    """
    REQUIRED_FIELDS = ('topology.architecture', 'topology.connectivity')
    RETURNED_FIELDS = ()
    VERSION = 'v1.0'

    def __init__( self, tag: int,
                  connectivity: Optional[bool] = True,
                  pick_aa: Optional[str] = None,
                  write2disc: Optional[str] = True):
        super(builder, self).__init__(tag)

        self.connectivity = connectivity
        self.pick_aa = pick_aa
        self.write2disc = write2disc

    def single_check( self, dummy: Dict ) -> Dict:
        kase = Case(dummy)

        # Check what it needs
        for itag in self.REQUIRED_FIELDS:
            if kase[itag] is None:
                raise NodeDataError(f'Field "{itag}" is required')

        # Include what keywords it adds (in this instance, nothing)
        # Here Nothing
        return kase.data

    def single_execute( self, data: Dict ) -> Dict:
        case = Case(data)

        # Apply connectivity?
        if self.connectivity:
            if case.connectivity_count > 1:
                raise NodeDataError('Only single connectivity cases can be build.')
            case = case.cast_absolute().apply_topologies()[0]
            ofile = case.connectivities_paths[0].joinpath('directed_sketch.pdb')
        else:
            case = case.cast_absolute()
            ofile = case.main_path.joinpath('architecture').joinpath('undirected_sketch.pdb')

        for i, j, sse in case:
            if 'atoms' in sse['metadata'] and not TBcore.get_option('system', 'overwrite'):
                self.log.debug(f'{case.name}.{sse["id"]} already has atoms defined\n')
                continue

            self.log.debug(f'Building coordinates for {case.name}.{sse["id"]}\n')
            case.data['topology']['architecture'][i][j] = self.make_structure(sse, self.pick_aa)

        if self.write2disc:
            ofile.parent.mkdir(parents=True, exist_ok=True)
            structure, _ = TButil.build_pdb_object( self.log, case.ordered_structures, 2 )

            TButil.plugin_filemaker(f'Writing structure {ofile}')
            structure.write(output_file=str(ofile), format='pdb', clean=True,
                            force=TBcore.get_option('system', 'overwrite'))

        return case

    def make_structure( self, sse: Dict, pick_aa: Optional[str] = None ) -> Case:
        """
        """

        structure = SSEArchitect(sse, type=sse['type'], pick_aa=pick_aa).pdb
        sse['metadata'].setdefault('atoms', None)
        sse['metadata']['atoms'] = list(structure[['auth_comp_id', 'auth_atom_id', 'auth_seq_id',
                                                   'Cartn_x', 'Cartn_y', 'Cartn_z']].values)
        return sse

        # # Loop MASTER is only applied to a Case with one single connectivity
        # if kase.connectivity_count != 1:
        #     err = f'{self.nodeID} can only be applied to one connectivity. '
        #     err += f'Current case contains a total of {kase.connectivity_count}.'
        #     raise NodeDataError(err)
        # # And has to be reoriented
        # if not kase.is_reoriented:
        #     self.log.debug('Topology was provided without oriented SSE -> orienting.')
        #     kase = kase.apply_topologies()[0]
        #
        # # Generate the folder tree for a single connectivity.
        # folders = kase.connectivities_paths[0].joinpath('loop_master')
        # folders.mkdir(parents=True, exist_ok=True)
        #
        # # Global step distance
        # loop_step = kase.cast_absolute()['configuration.defaults.distance.loop_step']
        #
        # # Output keys
        # kase.data.setdefault('metadata', {}).setdefault('loop_fragments', [])
        # kase.data.setdefault('metadata', {}).setdefault('loop_lengths', [])
        #
        # # Find steps: Each pair of secondary structure.
        # it = kase.connectivities_str[0].split('.')
        # steps = [it[i:i + 2] for i in range(0, len(it) - 1)]
        # lengths = kase.connectivity_len[0]
        # start = 1
        #
        # for i, sse in enumerate(steps):
        #     # 1. Make folders and files
        #     wfolder = folders.joinpath('loop{:02d}'.format(i + 1))
        #     wfolder.mkdir(parents=True, exist_ok=True)
        #     outfile = wfolder.joinpath('loop_master.jump{:02d}.pdb'.format(i + 1))
        #     masfile = outfile.with_suffix('.master')
        #     checkpoint = wfolder.joinpath('checkpoint.json')
        #
        #     # 2. Check if checkpoint exists, retrieve and skip
        #     reload = self.checkpoint_in(checkpoint)
        #     if reload is not None:
        #         kase.data['metadata']['loop_fragments'].append(reload)
        #         kase.data['metadata']['loop_lengths'].append(int(reload['edges']['loop']))
        #         start += (int(reload['edges']['sse1']) + int(reload['edges']['loop']))
        #         continue
        #
        #     # 3. Check hairpin
        #     # Get SSEs and identifiers
        #     sse1, sse2 = kase.get_sse_by_id(sse[0]), kase.get_sse_by_id(sse[1])
        #     sse1_name, sse2_name = sse1['id'], sse2['id']
        #     is_hairpin = self.check_hairpin(sse1_name, sse2_name)
        #
        #     if not masfile.is_file():
        #         # 4. Generate structures
        #         sse1, sse2 = TBstructure.build_pdb_object(self.log, [sse1, sse2], 5,
        #                                                   concat=False, outfile=outfile)
        #
        #         # 5. calculate expected loop length by loop_step
        #         Mdis, mdis = TBstructure.get_loop_length(self.log, sse1, sse2, loop_step, self.loop_range)
        #
        #         # 6. Run MASTER
        #         outfilePDS = outfile if outfile is not None else Path(outfile).with_suffix('.pds')
        #         # -> make PDS query
        #         cmd = TBMaster.createPDS(outfile, outfilePDS)
        #         self.log.debug(f'EXECUTE: {" ".join(cmd)}')
        #         run(cmd, stdout=DEVNULL)
        #         # -> run MASTER
        #         cmd = TBMaster.master_fixedgap(outfilePDS, self.pdsdb, masfile, mdis, Mdis, self.rmsd_cut)
        #         self.log.debug(f'EXECUTE: {" ".join(cmd)}')
        #         run(cmd, stdout=DEVNULL)
        #
        #         # 6. Minimize master data (pick top_loopsx3 lines to read and minimize the files)
        #         match_count = self.minimize_master_file(masfile)
        #
        #     # 7. Retrieve MASTER data
        #     dfloop = self.process_master_data(masfile, sse1_name, sse2_name, is_hairpin and self.harpins_2)
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
        #
        #     # 9. Save data in the Case
        #     kase.data['metadata']['loop_fragments'].append(loop_data)
        #     kase.data['metadata']['loop_lengths'].append(int(loopl))
        #
        #     start += (sse1l + loopl)
        #
        #     # 10. Checkpoint save
        #     self.checkpoint_out(checkpoint, loop_data)
        #
        # return kase

# def apply( cases: List[Case],
#            prtid: int,
#            connectivity: Optional[bool] = True,
#            pick_aa: Optional[str] = None,
#            **kwargs ) -> List[Case]:
#     """Create coordinate entities from the Case data.
#     """
#     TButil.plugin_title(__file__, len(cases))
#
#     for i, case in enumerate(cases):
#         cases[i] = case_apply(case, connectivity, pick_aa, write2disc=True)
#         cases[i] = cases[i].set_protocol_done(prtid)
#     return cases


# def case_apply( case: Case,
#                 connectivity: bool,
#                 pick_aa: Optional[str] = None,
#                 write2disc: Optional[bool] = False
#                 ) -> Case:
#     """
#     """
#     # Apply connectivity?
#     if connectivity:
#         if case.connectivity_count > 1:
#             raise ValueError('Only single connectivity cases can be build.')
#         case = case.cast_absolute().apply_topologies()[0]
#         ofile = case.connectivities_paths[0].joinpath('directed_sketch.pdb')
#     else:
#         case = case.cast_absolute()
#         ofile = case.main_path.joinpath('architecture').joinpath('undirected_sketch.pdb')
#
#     for i, j, sse in case:
#         if 'atoms' in sse['metadata'] and not TBcore.get_option('system', 'overwrite'):
#             if TBcore.get_option('system', 'verbose'):
#                 sys.stdout.write('{0}.{1} already has atoms defined\n'.format(case.name, sse['id']))
#             continue
#
#         if TBcore.get_option('system', 'verbose'):
#             sys.stdout.write('Building coordinates for {0}.{1}\n'.format(case.name, sse['id']))
#         case.data['topology']['architecture'][i][j] = make_structure(sse, pick_aa)
#
#     if write2disc:
#         ofile.parent.mkdir(parents=True, exist_ok=True)
#         structure, _ = TButil.build_pdb_object( case.ordered_structures, 2 )
#
#         TButil.plugin_filemaker('Writing structure {0}'.format(ofile))
#         structure.write(output_file=str(ofile), format='pdb', clean=True,
#                         force=TBcore.get_option('system', 'overwrite'))
#
#     return case


# def make_structure( sse: Dict, pick_aa: Optional[str] = None ) -> Case:
#     """
#     """
#
#     structure = SSEArchitect(sse, type=sse['type'], pick_aa=pick_aa).pdb
#     sse['metadata'].setdefault('atoms', None)
#     sse['metadata']['atoms'] = list(structure[['auth_comp_id', 'auth_atom_id', 'auth_seq_id',
#                                                'Cartn_x', 'Cartn_y', 'Cartn_z']].values)
#     return sse
