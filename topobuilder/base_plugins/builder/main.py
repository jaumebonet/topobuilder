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
