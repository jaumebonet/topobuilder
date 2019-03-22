# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. class:: ParametricStructure
"""
# Standard Libraries
import sys
from random import random
from bisect import bisect
from typing import Dict, Union, Optional

# External Libraries
import numpy as np
import pandas as pd
from SBI.structure import PDB, Frame3D
import SBI.structure.geometry as SBIgeo
from SBI.data import alphabet

# This Library
import topobuilder.core as TBcore


__all__ = ['SSEArchitect']


class SSEArchitect( object ):
    """Decides the correct type of secondary structure to build.
    """
    def __new__( cls, *args, **kwargs ):
        sse_type = kwargs.pop('type', None)
        if sse_type is None:
            raise AttributeError('A secondary structure type must be provided.')
        sse_type = sse_type.upper()

        if sse_type == 'H':
            return AlphaHelixArchitect(*args, **kwargs)
        elif sse_type == 'G':
            return Helix310Architect(*args, **kwargs)
        elif sse_type == 'I':
            return HelixPiArchitect(*args, **kwargs)
        elif sse_type == 'E':
            return FlatBetaArchitect(*args, **kwargs)
        else:
            raise ValueError('Unrecognized secondary structure type {}.'.format(sse_type))


class ParametricStructure( object ):

    _MONO, _PERIODE, _ROTATION, _AA_STAT = None, None, None, None

    def __init__( self, indata: Union[Frame3D, Dict], pick_aa: Optional[str] = None ):
        """
        """
        if isinstance(indata, Frame3D):
            self.pdb = indata
            self.desc = None
            self.reverse()
        elif isinstance(indata, dict):
            self.pdb = []
            self.desc = indata
            self.build(pick_aa)

    def build( self, pick_aa: Optional[str] = None ):
        """
        """
        if self._MONO is None or self._PERIODE is None:
            raise NotImplementedError()

        # 1. Locate center point for each residue we need to build
        vector_module = float(self._PERIODE * (self.desc['length'] - 1))
        upper_bound = np.copy(np.array([0., 0., 0.], dtype='float64')) + np.array([0, vector_module / 2, 0])
        points = [np.copy(upper_bound) - np.array([0, self._PERIODE * x, 0]) for x in range(self.desc['length'])]

        # 2. Build. For each point, we build one periode at [0, 0, 0]. Then, we rotate and then shift.
        self.pdb = []
        _MONO = pd.DataFrame(self._MONO).T
        for i, p in enumerate(points):
            coords = SBIgeo.rotate_degrees(_MONO.values, y=self._ROTATION * i)
            coords = SBIgeo.translate(coords, p)
            self.pdb.append(coords)
        self.pdb = np.vstack(self.pdb)
        # We want undirected structures to start always looking up
        self.pdb = SBIgeo.rotate_degrees(self.pdb, x=180)

        # Apply the case-defined placements for each structure
        if TBcore.get_option('system', 'debug'):
            sys.stdout.write('tilt: ' + str(self.desc['tilt']) + '\n')
            sys.stdout.write('move: ' + str(self.desc['coordinates']) + '\n')

        self.pdb = SBIgeo.rotate_degrees(self.pdb, x=self.desc['tilt']['x'],
                                         y=self.desc['tilt']['y'],
                                         z=self.desc['tilt']['z'])
        self.pdb = SBIgeo.translate(self.pdb, [self.desc['coordinates']['x'],
                                               self.desc['coordinates']['y'],
                                               self.desc['coordinates']['z']])

        # Prepare other data to create a coordinate entity
        resis = np.repeat(list(range(1, i + 2)), _MONO.shape[0])
        atoms = np.asarray([_MONO.index.values, ] * (i + 1)).flatten()

        # Prepare sequence
        sequence = []
        if pick_aa is not None:
            pick_aa = pick_aa if len(pick_aa) == 3 else alphabet.aminoacids1to3(pick_aa)
            sequence = [pick_aa, ] * self.desc['length']
        else:
            for _ in range(self.desc['length']):
                sequence.append(alphabet.aminoacids1to3(weighted_choice(self._AA_STAT)))
        sequence = np.repeat(np.asarray(sequence), _MONO.shape[0])

        self.pdb = PDB(pd.DataFrame(self.pdb, columns=["Cartn_x", "Cartn_y", "Cartn_z"])
                         .assign(auth_comp_id=sequence)
                         .assign(auth_atom_id=atoms).assign(auth_seq_id=resis)
                         .assign(id=list(range(1, self.pdb.shape[0] + 1))))

    def reverse( self ):
        """
        """
        raise NotImplementedError()


class AlphaHelixArchitect( ParametricStructure ):
    """
    """
    _PERIODE = 1.5
    _ROTATION = 100
    _MONO = {'N': [1.321, 0.841, -0.711],
             'CA': [2.300, 0.000, 0.000],
             'C': [1.576, -1.029, 0.870],
             'O': [1.911, -2.248, 0.871]}
    # CHOP780201 alpha-helix propensity AAindex (Chou-Fasman, 1978b)
    # TO 0: G -> 0.57; P -> 0.57
    _AA_STAT = [("A", 1.42), ("L", 1.21), ("R", 0.98), ("K", 1.16), ("N", 0.67),
                ("M", 1.45), ("D", 1.01), ("F", 1.13), ("C", 0.70), ("P", 0.00),
                ("Q", 1.11), ("S", 0.77), ("E", 1.51), ("T", 0.83), ("G", 0.00),
                ("W", 1.08), ("H", 1.00), ("Y", 0.69), ("I", 1.08), ("V", 1.06)]


class Helix310Architect( ParametricStructure ):
    """
    """
    _PERIODE = 2.0


class HelixPiArchitect( ParametricStructure ):
    """
    """
    _PERIODE = 1.1


class FlatBetaArchitect( ParametricStructure ):
    """
    """
    _PERIODE = 3.2
    _ROTATION = -180
    _MONO = {'N':  [-0.440,  1.200, -0.330],
             'CA': [-0.000,  0.000, -1.210],
             'C':  [-0.550, -1.200, -0.330],
             'O':  [-2.090, -1.300, -0.220]}
    # CHOP780202 beta-sheet propensity AAindex (Chou-Fasman, 1978b)
    # TO 0: G -> 0.75; P -> 0.55
    _AA_STAT = [("A", 0.83), ("L", 1.30), ("R", 0.93), ("K", 0.74), ("N", 0.89),
                ("M", 1.05), ("D", 0.54), ("F", 1.38), ("C", 1.19), ("P", 0.00),
                ("Q", 1.10), ("S", 0.75), ("E", 0.37), ("T", 1.19), ("G", 0.00),
                ("W", 1.37), ("H", 0.87), ("Y", 1.47), ("I", 1.60), ("V", 1.70)]


def weighted_choice(choices):
    values, weights = zip(*choices)
    total = 0
    cum_weights = []
    for w in weights:
        total += w
        cum_weights.append(total)
    x = random() * total
    i = bisect(cum_weights, x)
    return values[i]
