# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os
import sys
from pathlib import Path
from typing import List, Dict
from subprocess import run, DEVNULL

# External Libraries
import pandas as pd

# This Library
from logbook import Logger
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil


def make_folder_structure( case: Case, statmode: str ) -> Dict:
    """
    """

    # Generate the folder tree for a single connectivity.
    wfolder = case.connectivities_paths[0].joinpath('statistic')
    wfolder.mkdir(parents=True, exist_ok=True)

    # Generate internal folder
    thisfolder = wfolder.joinpath(statmode)
    thisfolder.mkdir(parents=True, exist_ok=True)

    return {'global': wfolder, 'main': thisfolder}


def make_directed_sketch( log: Logger, case: Case, folder: Path ) -> Path:
    """
    """

    pdbfile = folder.joinpath('directed_sketch.pdb')
    structure, _ = TButil.build_pdb_object(case.apply_topologies()[0].ordered_structures, 3)
    log.notice(f'New File: Writing structure {pdbfile}')
    structure.write(output_file=str(pdbfile), format='pdb', clean=True, force=True)
    return pdbfile


def count_single_master_matches( pdbfile: Path, folder: Path ) -> pd.DataFrame:
    """
    """

    createpds = TButil.createPDS(pdbfile)
    TButil.plugin_bash(createpds)
    run(createpds, stdout=DEVNULL)
    masters = TButil.master_best_each(pdbfile.with_suffix('.pds'), folder.joinpath('_master'), 5)

    unimaster = folder.joinpath('match.master')

    if not TBcore.get_option('slurm', 'use'):
        no_slurm(cmd, current_case_file, current_sse, unimaster, imaster, unidata.with_suffix(''))
    else:
        with_slurm(cmd, current_case_file, current_sse, unimaster, imaster, unidata)

    return {'matches': unimaster, 'stats': unidata, 'corrections': None,
            'layers': list(set([x[0] for x in current_sse.split('.')]))}


    data = submit_searches(masters, stepfolder, current_case_file, '.'.join([x['id'] for x in sses]))
