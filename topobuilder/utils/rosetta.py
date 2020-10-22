"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import textwrap
import math
from typing import Dict, Union

# External Libraries
from jinja2 import Template
from bs4 import BeautifulSoup

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore


__all__ = ['rosettascript', 'get_weight_patches', 'funfoldes', 'hybridize', 'constraint_design']


class ScriptPieces( dict ):
    def __add__( self, other ):
        data = ScriptPieces()
        for k in ['scorefxns', 'residueselectors', 'packerpalette', 'taskoperations', 'movemapfactory',
                  'simplemetrics', 'filters', 'movers', 'protocols', 'output']:
            if k in self or k in other:
                data.setdefault(k, [x for x in self.get(k, ['', ]) + other.get(k, ['', ]) if len(x) > 0])
        return data

def get_weight_patches():
    """
    """
    wts0 = textwrap.dedent("""\
    # score0 from rosetta++, used in stage 1 of the
    # ClassicAbinitio protocol.
    # Score 0 has a vdw weight of 1, in R++, but then it divides
    # the vdw score by 10 and rounds down to nearest integer.
    # Mini does not round down to the nearest integer.
    env     0.0
    pair    0.0
    cbeta   0.0
    vdw     0.1
    rg      0.0
    cenpack 0.0
    hs_pair 0.0
    ss_pair 0.0
    rsigma  0.0 """)

    wts0_patch = textwrap.dedent("""\
    env     = 0.0
    pair    = 0.0
    cbeta   = 0.0
    vdw     = 0.1
    rg      = 0.0
    cenpack = 0.0
    hs_pair = 0.0
    ss_pair = 0.0
    rsigma  = 0.0
    hbond_sr_bb = 0.3
    hbond_lr_bb = 0.7
    rsigma  = 0.2
    sheet   = 0.2
    ss_pair = 0.2
    hs_pair = 0.2 """)

    wts1 = textwrap.dedent("""\
    # score1 from rosetta++, used in stage 2 of ClassicAbinitio
    # protocol from Rosetta++.
    env     1.0
    pair    1.0
    cbeta   0.0
    vdw     1.0
    rg      0.0
    cenpack 0.0
    hs_pair 1.0
    ss_pair 0.3
    rsigma  0.0
    sheet   1.0
    STRAND_STRAND_WEIGHTS 1 11 """)

    wts1_patch = textwrap.dedent("""\
    env     = 1.0
    pair    = 1.0
    cbeta   = 0.0
    vdw     = 1.0
    rg      = 0.0
    cenpack = 0.0
    rsigma  = 1.1
    sheet   = 1.0
    ss_pair = 1.0
    hs_pair = 1.0
    hbond_sr_bb = 1.17
    hbond_lr_bb = 2.0
    STRAND_STRAND_WEIGHTS 1 3 """)

    wts2 = textwrap.dedent("""\
    # score2 from rosetta++, used in stage 3 of ClassicAbinitio protocol.
    env     1.0
    pair    1.0
    cbeta   0.25
    cenpack 0.5
    vdw     1.0
    rg      0.0
    hs_pair 1.0
    ss_pair 1.0
    rsigma  0.0
    sheet   1.0
    STRAND_STRAND_WEIGHTS 1 6 """)

    wts2_patch = textwrap.dedent("""\
    env     = 1.0
    pair    = 1.0
    cbeta   = 0.25
    cenpack = 0.5
    vdw     = 1.0
    rg      = 0.0
    rsigma  = 0.7
    sheet   = 0.7
    ss_pair = 0.7
    hs_pair = 0.7
    hbond_sr_bb = 1.17
    hbond_lr_bb = 2.0
    angle_constraint = 0.3
    dihedral_constraint = 0.3
    STRAND_STRAND_WEIGHTS 1 3 """)

    wts3 = textwrap.dedent("""\
    # score3 from rosetta++, used in stage 4 of the ClassicAbinitio protocol.
    env     1.0
    pair    1.0
    cbeta   1.0
    vdw     1.0
    rg      3.0
    cenpack 1.0
    hs_pair 1.0
    ss_pair 1.0
    rsigma  1.0
    sheet   1.0
    STRAND_STRAND_WEIGHTS 1 6 """)

    wts3_patch = textwrap.dedent("""\
    env     = 1.0
    pair    = 1.0
    cbeta   = 1.0
    vdw     = 1.0
    rg      = 2.0
    cenpack = 1.0
    hbond_sr_bb = 1.17
    rsigma  = 1.17
    sheet   = 1.17
    ss_pair = 1.17
    hs_pair = 1.17
    hbond_lr_bb = 2.0
    angle_constraint = 0.3
    dihedral_constraint = 0.3
    STRAND_STRAND_WEIGHTS 1 2 """)

    wts5 = textwrap.dedent("""\
    # score5.wts, used in stage 3 of ClassicAbinitio protocol
    env     1.0
    pair    1.0
    cbeta   0.25
    cenpack 0.5
    hs_pair 1.0
    ss_pair 1.0
    rsigma  0.0
    sheet   1.17
    rg      0.0
    vdw     1.0
    STRAND_STRAND_WEIGHTS 1 6 """)

    wts5_patch = textwrap.dedent("""\
    env     = 1.0
    pair    = 1.0
    cbeta   = 0.25
    cenpack = 0.5
    rg      = 0.0
    vdw     = 1.0
    rsigma  = 1.17
    sheet   = 1.17
    ss_pair = 1.17
    hs_pair = 1.17
    hbond_sr_bb = 1.17
    hbond_lr_bb = 2.0
    angle_constraint = 0.3
    dihedral_constraint = 0.3
    STRAND_STRAND_WEIGHTS 1 3 """)

    wts_pieces = [wts0, wts0_patch,
                  wts1, wts1_patch,
                  wts2, wts2_patch,
                  wts3, wts3_patch,
                  wts5, wts5_patch,]
    return wts_pieces


def constraint_minimization(  case: Case, natbias: float ) -> ScriptPieces:
    """
    """
    scorefxns = [textwrap.dedent("""\
    <ScoreFunction name="sfxn_cstmin" weights="ref2015">
        <Reweight scoretype="atom_pair_constraint" weight="1.0" />
        <Reweight scoretype="angle_constraint" weight="1.0" />
        <Reweight scoretype="dihedral_constraint" weight="1.0" />
        <Reweight scoretype="hbond_lr_bb" weight="2.0" />
    </ScoreFunction>
    """), ]

    residueselectors = [SELECTOR_SecondaryStructure('sse_cstmin', case), ]

    filters = [textwrap.dedent("""\
    <RmsdFromResidueSelectorFilter name="rmsd_cstmin" reference_selector="sse_cstmin"
            reference_name="eminPose_cstmin" query_selector="sse_cstmin" confidence="0." />
    """), ]

    movers = [textwrap.dedent("""\
    <SavePoseMover name="spose_cstmin" reference_name="eminPose_cstmin" restore_pose="0" />
    <AddConstraints name="cst_cstmin" >
        <SegmentedAtomPairConstraintGenerator name="cst_seg_cstmin" residue_selector="sse_cstmin" >
            <Outer sd="2.0" weight="1." ca_only="1"
             use_harmonic="1" unweighted="0" max_distance="40" />
        </SegmentedAtomPairConstraintGenerator>
        <AutomaticSheetConstraintGenerator name="cst_sheet_cstmin" sd="2.0" distance="6.1" />
    </AddConstraints>
    <MinMover name="fast_cstmin" scorefxn="sfxn_cstmin" chi="1" bb="1" />
    """), MOVER_SetSecStructEnergies( 'ssse_cstmin', 'sfxn_cstmin', natbias, case )]

    protocols = [textwrap.dedent("""\
    <Add mover="spose_cstmin" />
    <Add mover="cst_cstmin" />
    <Add mover="ssse_cstmin" />
    <Add mover="fast_cstmin" />
    <Add filter="rmsd_cstmin" />
    """), ]

    bf = PROTOCOL_BasicFilters(case, '_cstmin')
    return ScriptPieces({'scorefxns': scorefxns, 'movers': movers, 'filters': filters,
                         'residueselectors': residueselectors, 'protocols': protocols}) + bf


def constraint_design( case: Case, natbias: float, layer_design: bool = True ) -> ScriptPieces:
    """
    """
    scorefxns = [textwrap.dedent("""\
    <ScoreFunction name="sfxn_cstdes" weights="ref2015">
        <Reweight scoretype="atom_pair_constraint" weight="1.0" />
        <Reweight scoretype="angle_constraint" weight="1.0" />
        <Reweight scoretype="dihedral_constraint" weight="1.0" />
        <Reweight scoretype="hbond_lr_bb" weight="2.0" />
    </ScoreFunction>
    <ScoreFunction name="sfxn_cstdes_cart" weights="ref2015_cart">
        <Reweight scoretype="atom_pair_constraint" weight="1.0" />
        <Reweight scoretype="angle_constraint" weight="1.0" />
        <Reweight scoretype="dihedral_constraint" weight="1.0" />
        <Reweight scoretype="hbond_lr_bb" weight="2.0" />
    </ScoreFunction>
    """), ]

    residueselectors = [SELECTOR_SecondaryStructure('sse_cstdes', case), ]

    filters = [textwrap.dedent("""\
    <RmsdFromResidueSelectorFilter name="rmsd_cstdes" reference_selector="sse_cstdes"
            reference_name="eminPose_cstdes" query_selector="sse_cstdes" confidence="0." />
    """), ]

    movers = [textwrap.dedent("""\
    <SavePoseMover name="spose_cstdes" reference_name="eminPose_cstdes" restore_pose="0" />
    <AddConstraints name="cst_cstdes" >
        <SegmentedAtomPairConstraintGenerator name="cst_seg_cstdes" residue_selector="sse_cstdes" >
            <Outer sd="2.0" weight="1.0" ca_only="1" use_harmonic="1" unweighted="0" max_distance="40" />
        </SegmentedAtomPairConstraintGenerator>
        <!--AutomaticSheetConstraintGenerator name="cst_sheet_cstdes" sd="2.0" distance="6.1" /-->
    </AddConstraints>
    <RemoveConstraints name="rm_cstdes" constraint_generators="cst_seg_cstdes" />"""),
    MOVER_SetSecStructEnergies( 'ssse_cstdes', 'sfxn_cstdes', natbias, case ), MOVER_SetSecStructEnergies( 'ssse_cstdes_cart', 'sfxn_cstdes_cart', natbias, case ),
              textwrap.dedent("""\
    <FastDesign name="design_cstdes" scorefxn="sfxn_cstdes_cart" relaxscript="MonomerDesign2019" task_operations="layer_design" ramp_down_constraints="false" repeats="5" dualspace="true"/>
    """) if layer_design else textwrap.dedent("""\
    <FastDesign name="design_cstdes" scorefxn="sfxn_cstdes_cart" relaxscript="MonomerDesign2019" ramp_down_constraints="false" repeats="5" dualspace="true"/>"""),
              textwrap.dedent("""\
    <FastRelax name="cst_cstrel" scorefxn="sfxn_cstdes_cart" repeats="5" cartesian="true"/>""")]

    protocols = [textwrap.dedent("""\
    <Add mover="ssse_cstdes" />
    <Add mover="spose_cstdes" />
    <Add mover="cst_cstdes" />
    <Add mover="design_cstdes" />
    <Add mover="rm_cstdes" />
    <Add mover="ssse_cstdes_cart" />
    <Add mover="cst_cstrel" />
    <Add filter="rmsd_cstdes" />
    """), ]

    ld = PROTOCOL_LayerDesign(case) if layer_design else ScriptPieces()
    bf = PROTOCOL_BasicFilters(case, '_cstdes')
    return ScriptPieces({'scorefxns': scorefxns, 'movers': movers, 'filters': filters,
                         'residueselectors': residueselectors, 'protocols': protocols}) + ld + bf


def funfoldes( case: Case ) -> str:
    """
    """
    mid = 2 # math.floor(len(case.secondary_structure) / 2)
    residueselectors = [textwrap.dedent("""<Index name="piece_ffd" resnums="{}-{}" />""").format(mid - 1, mid + 1),
                        SELECTOR_SecondaryStructure('sse_ffd', case)]

    filters = [textwrap.dedent("""\
        <RmsdFromResidueSelectorFilter name="rmsd_ffd" reference_selector="sse_ffd"
            reference_name="sketchPose_ffd" query_selector="sse_ffd" confidence="0." />""")]

    movers = [MOVER_PeptideStubMover('add_loops_ffd', case), textwrap.dedent("""\
        <SavePoseMover name="save_ffd" reference_name="sketchPose_ffd" restore_pose="0" />
        <StructFragmentMover name="makeFrags_ffd" prefix="frags" small_frag_file="{}" large_frag_file="{}" />
        <AddConstraints name="foldingCST_ffd" >
            <SegmentedAtomPairConstraintGenerator name="foldCST" residue_selector="sse_ffd" >
                <!--Inner sd="1.2" weight="1." ca_only="1"
                    use_harmonic="true" unweighted="false" min_seq_sep="4" /-->
                <Outer sd="2" weight="2." ca_only="1"
                    use_harmonic="true" unweighted="false"  max_distance="40" />
            </SegmentedAtomPairConstraintGenerator>
            <!--AutomaticSheetConstraintGenerator name="sheetCST" sd="2.0" distance="6.1" /-->
        </AddConstraints>
        <NubInitioMover name="FFL_ffd" fragments_id="frags" template_motif_selector="piece_ffd" rmsd_threshold="10" correction_weights="0">
        """).format(*case['metadata.fragments.files']), textwrap.dedent("""\
            <Nub reference_name="sketchPose_ffd" residue_selector="piece_ffd" >
                <Segment order="1" n_term_flex="2" c_term_flex="1" editable="1,2,3"/></Nub>
            """), textwrap.dedent("""</NubInitioMover>""")]

    protocols = [textwrap.dedent("""\
        <Add mover="add_loops_ffd" />
        <Add mover="save_ffd" />
        <Add mover="makeFrags_ffd" />
        <Add mover="foldingCST_ffd" />
        <Add mover="FFL_ffd" />
        <Add filter="rmsd_ffd" />""")]

    with TBcore.on_option_value('psipred', 'script', None):
        bf = PROTOCOL_BasicFilters(case, '_ffd')
    return ScriptPieces({'movers': movers, 'filters': filters,
                         'residueselectors': residueselectors, 'protocols': protocols}) + bf


def hybridize( case: Case, template: str, natbias: float ) -> str:
    """
    """
    #mid = math.floor(len(case.secondary_structure) / 2)
    mid = 2
    scorefxns = [textwrap.dedent("""\
    <ScoreFunction name="fa" weights="ref2015">
      <Reweight scoretype="coordinate_constraint" weight="1" />
      <Reweight scoretype="atom_pair_constraint" weight="1" />
      <Reweight scoretype="dihedral_constraint" weight="1" />
      <Reweight scoretype="angle_constraint" weight="1" />
      <Reweight scoretype="rsigma" weight="1.17" />
      <Reweight scoretype="sheet" weight="2.0" />
      <Reweight scoretype="ss_pair" weight="2.0" />
      <Reweight scoretype="hs_pair" weight="2.0" />
      <Reweight scoretype="hbond_sr_bb" weight="1.17" />
      <Reweight scoretype="hbond_lr_bb" weight="2.0" />
    </ScoreFunction>
    <ScoreFunction name="stage1" weights="score3">
      <Reweight scoretype="coordinate_constraint" weight="1" />
      <Reweight scoretype="atom_pair_constraint" weight="1" />
      <Reweight scoretype="dihedral_constraint" weight="1" />
      <Reweight scoretype="angle_constraint" weight="1" />
      <Reweight scoretype="rsigma" weight="1.17" />
      <Reweight scoretype="sheet" weight="2.0" />
      <Reweight scoretype="ss_pair" weight="2.0" />
      <Reweight scoretype="hs_pair" weight="2.0" />
      <Reweight scoretype="hbond_sr_bb" weight="1.17" />
      <Reweight scoretype="hbond_lr_bb" weight="2.0" />
    </ScoreFunction>
    <ScoreFunction name="stage2" weights="score4_smooth_cart">
      <Reweight scoretype="coordinate_constraint" weight="1" />
      <Reweight scoretype="atom_pair_constraint" weight="1" />
      <Reweight scoretype="dihedral_constraint" weight="1" />
      <Reweight scoretype="angle_constraint" weight="1" />
      <Reweight scoretype="rsigma" weight="1.17" />
      <Reweight scoretype="sheet" weight="2.0" />
      <Reweight scoretype="ss_pair" weight="2.0" />
      <Reweight scoretype="hs_pair" weight="2.0" />
      <Reweight scoretype="hbond_sr_bb" weight="1.17" />
      <Reweight scoretype="hbond_lr_bb" weight="2.0" />
    </ScoreFunction>
    """), ]

    residueselectors = [textwrap.dedent("""\
    <Index name="piece_ffd" resnums="{}-{}" />
    """).format(mid - 1, mid + 1), SELECTOR_SecondaryStructure('sse_ffd', case)]

    filters = [textwrap.dedent("""\
        <RmsdFromResidueSelectorFilter name="rmsd_ffd" reference_selector="sse_ffd"
            reference_name="sketchPose_ffd" query_selector="sse_ffd" confidence="0." />""")]

    movers = [MOVER_PeptideStubMover('add_loops_ffd', case, 'VAL'),
    #MOVER_SetSecStructEnergies('sse_energies', 'fa', natbias, case),
    textwrap.dedent("""\
        <SavePoseMover name="save_ffd" reference_name="sketchPose_ffd" restore_pose="0" />
        <AddConstraints name="foldingCST_ffd" >
            <SegmentedAtomPairConstraintGenerator name="foldCST" residue_selector="sse_ffd" >
                <!--Inner sd="1.2" weight="1." ca_only="1"
                    use_harmonic="true" unweighted="false" min_seq_sep="4" /-->
                <Outer sd="2" weight="2." ca_only="1"
                    use_harmonic="true" unweighted="false"  max_distance="40" />
            </SegmentedAtomPairConstraintGenerator>
            <AutomaticSheetConstraintGenerator name="sheetCST" sd="2.0" distance="6.1" />
        </AddConstraints>"""),
        MOVER_SetSecStructEnergies( 'ssse_stage1', 'stage1', natbias, case ),
        MOVER_SetSecStructEnergies( 'ssse_stage2', 'stage2', natbias, case ),
        MOVER_SetSecStructEnergies( 'ssse_fa', 'fa', natbias, case ),
    textwrap.dedent("""\
        <Hybridize name="hybridize" stage1_scorefxn="stage1" stage2_scorefxn="stage2" fa_scorefxn="fa"
                   batch="1" stage1_increase_cycles="1.0" stage2_increase_cycles="1.0" linmin_only="0"
                   realign_domains="0" keep_pose_constraint="1" csts_from_frags="0" max_registry_shift="1">
                   <Fragments three_mers="{}" nine_mers="{}"/>
                   <Template pdb="{}" cst_file="AUTO" weight="1.000" />
        </Hybridize>""").format(*case['metadata.fragments.files'], template)]

    protocols = [textwrap.dedent("""\
        <Add mover="add_loops_ffd" />
        <Add mover="save_ffd" />
        <Add mover="foldingCST_ffd" />
        <Add mover="ssse_stage1" />
        <Add mover="ssse_stage2" />
        <Add mover="ssse_fa" />
        <Add mover="hybridize" />
        <Add filter="rmsd_ffd" />""")]

    with TBcore.on_option_value('psipred', 'script', None):
        bf = PROTOCOL_BasicFilters(case, '_ffd')
    return ScriptPieces({'scorefxns': scorefxns, 'movers': movers, 'filters': filters,
                         'residueselectors': residueselectors,
                         'protocols': protocols}) + bf


def rosettascript( data: Dict ) -> str:
    """
    """
    if 'output' in data and isinstance(data['output'], list):
        data['output'] = data['output'][0]
    content = BeautifulSoup(Template(textwrap.dedent("""\
    <ROSETTASCRIPTS>
        {% if scorefxns %}<SCOREFXNS>{% for item in scorefxns %}{{ item }}{% endfor %}</SCOREFXNS>{% endif %}
        {% if residueselectors %}
            <RESIDUE_SELECTORS>{% for item in residueselectors %}{{ item }}{% endfor %}</RESIDUE_SELECTORS>
        {% endif %}
        {% if packerpalettes %}<PACKER_PALETTES>{% for item in packerpalettes %}{{item}}{% endfor %}</PACKER_PALETTES>{% endif %}
        {% if taskoperations %}<TASKOPERATIONS>{% for item in taskoperations %}{{ item }}{% endfor %}</TASKOPERATIONS>{% endif %}
        {% if movemapfactory %}
            <MOVE_MAP_FACTORIES>{% for item in movemapfactory %}{{ item }}{% endfor %}</MOVE_MAP_FACTORIES>
        {% endif %}
        {% if simplemetrics %}<SIMPLE_METRICS>{% for item in simplemetrics %}{{ item }}{% endfor %}</SIMPLE_METRICS>{% endif %}
        {% if filters %}<FILTERS>{% for item in filters %}{{ item }}{% endfor %}</FILTERS>{% endif %}
        {% if movers %}<MOVERS>{% for item in movers %}{{ item }}{% endfor %}</MOVERS>{% endif %}
        {% if protocols %}<PROTOCOLS>{% for item in protocols %}{{ item }}{% endfor %}</PROTOCOLS>{% endif %}
        {% if output %}<OUTPUT scorefxn="{{ output }}"/>{% endif %}
    </ROSETTASCRIPTS>
    """)).render(data), 'xml').contents
    return '\n'.join(x.prettify() for x in content)


def SELECTOR_SecondaryStructure( name: str,
                                 sse: Union[str, Case],
                                 sse_type: str = 'HE',
                                 terminal_loops: bool = False
                                 ) -> str:
    """
    """
    if isinstance(sse, Case):
        sse = sse.secondary_structure

    return textwrap.dedent("""\
        <SecondaryStructure name="{}" overlap="0" minE="1"  minH="1" ss="{}" include_terminal_loops="{}"
        use_dssp="0" pose_secstruct="{}" />""").format(name, sse_type, int(terminal_loops), sse)


def MOVER_SetSecStructEnergies( name: str, score: str, natbias: float, case: Case ) -> str:
    """
    """
    data = dict(zip(['ss_pair', 'hh_pair', 'hss_triplets'], case.sse_pairing))
    data['sse'] = case.secondary_structure
    data['score'] = score
    data['name'] = name
    data['natbias'] = natbias
    return Template(textwrap.dedent("""\
        <SetSecStructEnergies name="{{name}}" scorefxn="{{score}}"
            secstruct="{{sse}}" use_dssp="0"
            {% if hh_pair|length > 0 %}hh_pair="{{hh_pair}}"{% endif %}
            {% if ss_pair|length > 0 %}ss_pair="{{ss_pair}}"{% endif %}
            {% if hss_triplets|length > 0 %}hss_triplets="{{hss_triplets}}"{% endif %}
            {% if ss_pair|length > 0 %}natbias_ss="{{natbias}}"{% endif %}
            {% if hh_pair|length > 0 %}natbias_hh="{{natbias}}"{% endif %}
            {% if hss_triplets|length > 0 %}natbias_hs="{{natbias}}"{% endif %}
        />""")).render(data)


def MOVER_PeptideStubMover( name: str, case: Case, residue: str = 'VAL' ) -> str:
    """
    """
    return Template(textwrap.dedent("""\
        <PeptideStubMover name="{{name}}" reset="false">
        {% for item in insert %}
        <Insert resname="{{residue}}" repeat="1" jump="false" anchor_rsd="{{item}}" anchor_atom="C" connecting_atom="N" />
        {% endfor %}
        </PeptideStubMover>""")).render({'insert': [i for i, ltr in enumerate(case.secondary_structure) if ltr == 'L'],
                                         'name': name, 'residue': residue})


def PROTOCOL_LayerDesign( case: Case ) -> ScriptPieces:
    """
    """
    residueselectors = [textwrap.dedent("""\
        <Layer name="surface" select_core="0" select_boundary="0" select_surface="1" use_sidechain_neighbors="1"/>
        <Layer name="boundary" select_core="0" select_boundary="1" select_surface="0" use_sidechain_neighbors="1"/>
        <Layer name="core" select_core="1" select_boundary="0" select_surface="0" use_sidechain_neighbors="1"/>"""),
                        SELECTOR_SecondaryStructure('sheet', case, 'E'),
                        SELECTOR_SecondaryStructure('entire_helix', case, 'H'),
                        SELECTOR_SecondaryStructure('entire_loop', case, 'L', True), textwrap.dedent("""\
        <And name="helix_cap" selectors="entire_loop">
            <PrimarySequenceNeighborhood lower="1" upper="0" selector="entire_helix"/></And>
        <And name="helix_start" selectors="entire_helix">
            <PrimarySequenceNeighborhood lower="0" upper="1" selector="helix_cap"/></And>
        <And name="helix" selectors="entire_helix"><Not selector="helix_start"/></And>
        <And name="loop" selectors="entire_loop"><Not selector="helix_cap"/></And>
        """)]

    taskoperations = [textwrap.dedent("""\
        <DesignRestrictions name="layer_design">
            <Action selector_logic="surface AND helix_start" aas="DEHKPQR"/>
            <Action selector_logic="surface AND helix" aas="EHKQR"/>
            <Action selector_logic="surface AND sheet" aas="EHKNQRST"/>
            <Action selector_logic="surface AND loop" aas="DEGHKNPQRST"/>
            <Action selector_logic="boundary AND helix_start" aas="ADEHIKLMNPQRSTVWY"/>
            <Action selector_logic="boundary AND helix" aas="ADEHIKLMNQRSTVWY"/>
            <Action selector_logic="boundary AND sheet" aas="DEFHIKLMNQRSTVWY"/>
            <Action selector_logic="boundary AND loop" aas="ADEFGHIKLMNPQRSTVWY"/>
            <Action selector_logic="core AND helix_start" aas="AFILMPVWY"/>
            <Action selector_logic="core AND helix" aas="AFILMVWY"/>
            <Action selector_logic="core AND sheet" aas="FILMVWY"/>
            <Action selector_logic="core AND loop" aas="AFGILMPVWY"/>
            <Action selector_logic="helix_cap" aas="DNST"/>
        </DesignRestrictions>"""), ]

    return ScriptPieces({'residueselectors': residueselectors, 'taskoperations': taskoperations})


def PROTOCOL_BasicFilters( case: Case, suffix: str = '' ) -> ScriptPieces:
    """
    """
    sse = case.secondary_structure

    scorefxns = textwrap.dedent("""\
        <ScoreFunction name="bb_only" weights="empty.wts" >
          <Reweight scoretype="fa_rep" weight="0.1" />
          <Reweight scoretype="fa_atr" weight="0.2" />
          <Reweight scoretype="hbond_sr_bb" weight="2.0" />
          <Reweight scoretype="hbond_lr_bb" weight="2.0" />
          <Reweight scoretype="rama_prepro" weight="0.45" />
          <Reweight scoretype="omega" weight="0.4" />
          <Reweight scoretype="p_aa_pp" weight="0.6" />
        </ScoreFunction>
    """)

    residueselectors = textwrap.dedent("""\
        <True name="full_pose" />
        <Layer name="surface{suffix}" select_core="0" select_boundary="0" select_surface="1" use_sidechain_neighbors="1"/>
        <Layer name="boundary{suffix}" select_core="0" select_boundary="1" select_surface="0" use_sidechain_neighbors="1"/>
        <Layer name="core{suffix}" select_core="1" select_boundary="0" select_surface="0" use_sidechain_neighbors="1"/>
        """).format(suffix=suffix)

    filters = textwrap.dedent("""\
    <PackStat name="pack{suffix}" confidence="0." />
    <CavityVolume name="cav_vol{suffix}" confidence="0." />
    <SecondaryStructure name="sse_match{suffix}" ss="{sse1}" compute_pose_secstruct_by_dssp="true" confidence="0." />
    <ScorePoseSegmentFromResidueSelectorFilter name="bbscore{suffix}" confidence="0"
    residue_selector="full_pose" scorefxn="bb_only" />
    """).format(sse1=sse, suffix=suffix)

    movers = [textwrap.dedent("""\
    <LabelPoseFromResidueSelectorMover name="labelcore{suffix}" property="CORE" residue_selector="core{suffix}" />
    <LabelPoseFromResidueSelectorMover name="labelboundary{suffix}" property="BOUNDARY" residue_selector="boundary{suffix}" />
    <LabelPoseFromResidueSelectorMover name="labelsurface{suffix}" property="SURFACE" residue_selector="surface{suffix}" />
    <DisplayPoseLabelsMover name="labeldump{suffix}" use_dssp="1" write="1" />
    """).format(suffix=suffix), ]
    if TBcore.get_option('psipred', 'script', in_path_none=True) is not None:
        movers.append(textwrap.dedent("""\
        <WriteSSEMover name="sse_report{suffix}" cmd="{psipred}" dssp="1" write_phipsi="1" />
        """).format(suffix=suffix, psipred=TBcore.get_option('psipred', 'script')))
    else:
        movers.append(textwrap.dedent("""\
        <WriteSSEMover name="sse_report{suffix}" dssp="1" write_phipsi="1" />
        """).format(suffix=suffix))

    protocols = textwrap.dedent("""\
    <Add mover="labelcore{suffix}"/>
    <Add mover="labelboundary{suffix}"/>
    <Add mover="labelsurface{suffix}"/>
    <Add mover="labeldump{suffix}"/>
    <Add mover="sse_report{suffix}"/>
    <Add filter="pack{suffix}" />
    <Add filter="cav_vol{suffix}" />
    <Add filter="sse_match{suffix}" />
    """).format(suffix=suffix)

    return ScriptPieces({'scorefxns': [scorefxns, ], 'residueselectors': [residueselectors, ], 'filters': [filters, ],
                         'movers': movers, 'protocols': [protocols, ]})
