"""The opt-in closed-shell check: ``validate_smiles_detailed(..., closed_shell=True)``.

Covers detection (warnings), the add-H / add-(negative)charge repair
suggestion, the exemptions (metals, aromatic rings, disconnected ions), that
``valid`` is never touched, and that the option defaults off and is
byte-identical when off.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

import pytest
from rdkit import Chem, RDLogger

from src import validate_smiles_detailed

RDLogger.DisableLog("rdApp.*")

_DATA = os.path.join(os.path.dirname(__file__), "../../data")


# --------------------------------------------------------------------------- #
# Default is off, and off is byte-identical to the option not existing
# --------------------------------------------------------------------------- #

RADICAL_SMILES = [
    "[CH3]", "C[N]C", "C[O]", "CC(=O)[O]", "c1ccc([O])cc1",
    "C[N+](=O)[O]", "[SH]", "CN([O])C",
]


@pytest.mark.parametrize("smiles", RADICAL_SMILES)
def test_default_is_off(smiles):
    default = validate_smiles_detailed(smiles)
    explicit_off = validate_smiles_detailed(smiles, closed_shell=False)
    assert default.warnings == [] and default.repaired_smiles is None
    assert default == explicit_off


@pytest.mark.parametrize("smiles", RADICAL_SMILES)
def test_off_is_byte_identical_to_on(smiles):
    """closed_shell=False never warns/repairs on inputs where closed_shell=True does."""
    off = validate_smiles_detailed(smiles, closed_shell=False)
    on = validate_smiles_detailed(smiles, closed_shell=True)
    assert off.warnings == []
    assert on.warnings != []
    assert off.valid == on.valid  # the option never changes the verdict


# --------------------------------------------------------------------------- #
# Never changes `valid`
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("smiles", RADICAL_SMILES)
def test_closed_shell_never_changes_valid(smiles):
    assert validate_smiles_detailed(smiles, closed_shell=True).valid


def test_closed_shell_does_not_mask_a_real_chemistry_failure():
    # genuinely over-valent AND, in a disconnected fragment, a radical.
    result = validate_smiles_detailed("[NH4][Pt]([NH4])(Cl)Cl.[CH3]", closed_shell=True)
    assert not result.valid and result.tier == "chemistry"
    # the radical still gets its warning...
    assert any("open-shell" in w for w in result.warnings)
    # ...but the repair is the real (over-valence) one, not conflated with it.
    assert result.repaired_smiles == "[NH3][Pt]([NH3])(Cl)Cl.[CH3]"


# --------------------------------------------------------------------------- #
# Missing H on C / N / O: detected, warned, repaired
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("smiles,repaired", [
    ("[CH3]", "C"),            # methyl radical -> methane
    ("C[N]C", "CNC"),          # dialkylaminyl radical -> dimethylamine
    ("C[O]", "CO"),            # alkoxyl radical -> methanol (aliphatic O: pKa 16, prefers H)
    ("[SH]", "S"),             # thiyl radical -> H2S
])
def test_missing_h_detected_and_repaired(smiles, repaired):
    result = validate_smiles_detailed(smiles, closed_shell=True)
    assert result.valid  # never an error
    assert any("open-shell" in w for w in result.warnings)
    assert result.repaired_smiles == repaired
    assert result.repairs
    # the repair is itself a real, closed-shell, RDKit-accepted molecule
    assert validate_smiles_detailed(repaired).valid
    assert Chem.MolFromSmiles(repaired) is not None


# --------------------------------------------------------------------------- #
# Missing negative charge on carboxylate / phenolate / nitro O
# --------------------------------------------------------------------------- #

def test_carboxylate_prefers_negative_charge():
    # pKa(carboxylic acid) ~4.2 << 7.4: deprotonated (anionic) at body pH.
    result = validate_smiles_detailed("CC(=O)[O]", closed_shell=True)
    assert result.valid
    assert result.repaired_smiles == "CC(=O)[O-]"
    assert any("set charge -1" in r and "carboxylate" in r for r in result.repairs)


def test_nitro_prefers_negative_charge():
    result = validate_smiles_detailed("C[N+](=O)[O]", closed_shell=True)
    assert result.valid
    assert result.repaired_smiles == "C[N+](=O)[O-]"
    assert any("set charge -1" in r and "nitro" in r for r in result.repairs)


def test_phenolate_oxygen_is_detected_and_repaired():
    # pKa(phenol) ~10 > 7.4: mostly *protonated* (neutral) at body pH, unlike
    # carboxylate/nitro O above -- physiological_ph picks the H-repair here,
    # the same way it already keeps aniline-type N neutral (docs/validation.md).
    result = validate_smiles_detailed("c1ccc([O])cc1", closed_shell=True)
    assert result.valid
    assert any("open-shell" in w and "O" in w for w in result.warnings)
    assert result.repaired_smiles == "c1ccc(O)cc1"
    assert any("phenolic" in r for r in result.repairs)
    assert Chem.MolFromSmiles(result.repaired_smiles) is not None


def test_physiological_ph_off_prefers_h_for_oxygen_too():
    # Without physiological_ph, there is no pKa table consulted: H wins.
    result = validate_smiles_detailed("CC(=O)[O]", closed_shell=True, physiological_ph=False)
    assert result.valid
    assert result.repaired_smiles == "CC(=O)O"


# --------------------------------------------------------------------------- #
# Aromatic atoms: no false radicals
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("smiles", [
    "c1ccccc1",            # benzene
    "c1ccncc1",            # pyridine
    "c1cc[nH+]cc1",        # pyridinium
    "c1ccc(O)cc1",         # phenol, written neutral
    "c1ccc2ccccc2c1",      # naphthalene (ring-fusion carbons)
])
def test_aromatic_atoms_are_not_false_radicals(smiles):
    result = validate_smiles_detailed(smiles, closed_shell=True)
    assert result.valid
    assert result.warnings == []
    assert result.repaired_smiles is None


# --------------------------------------------------------------------------- #
# Metals and d/f-block atoms are exempt
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("smiles", [
    "[Fe]",                        # bare d-block atom
    "[NH3][Pt]([NH3])(Cl)Cl",      # dative N/O donors to a d-block centre
    "Cl[Pt]1(Cl)[NH2]CC[NH2]1",
    "Cl[Cu]Cl",
])
def test_metal_complexes_are_exempt(smiles):
    """d/f-block atoms get no valence model (``_is_d_f_block``), so they can
    never be open-shell under this check -- matching valence_violations."""
    result = validate_smiles_detailed(smiles, closed_shell=True)
    assert result.valid
    assert result.warnings == []
    assert result.repaired_smiles is None


@pytest.mark.parametrize("smiles", [
    "CC(=O)O[Hg]c1ccccc1",                      # bare ester O donor -> Hg
    "N[Pt](N)(Cl)Cl",                            # bare amine N donor -> Pt (x2)
    "O=[N+]([O-])O[Cd]O[N+](=O)[O-]",            # bare O donor -> Cd (x2)
    "N#C[Fe-2](C#N)(C#N)(C#N)(C#N)N=O",          # nitrosyl N donor -> Fe
    "O[Fe]=O.O[Fe]=O.[Fe]",                      # bare O donor -> Fe (x2)
    "CCOC(=O)/C=C(\\C)O[Ti](O/C(C)=C/C(=O)OCC)(OC(C)C)OC(C)C",  # 4 bare O donors -> Ti
])
def test_bare_donor_atoms_bonded_to_metals_are_not_false_radicals(smiles):
    """Regression (found during the real-dataset false-positive scan,
    2026-10-17): an organic-subset (bracket-less) N/O donor bonded to a
    d/f-block metal was wrongly flagged as open-shell. fill_valence (run by
    _perceive, before this check ever sees the molecule) fills that donor's
    implicit H by counting *every* bond at face value, including the one to
    the metal -- but our own dative-bond convention excludes that same bond
    from valence. Trusting fill_valence's hcount therefore understated the
    atom's true total by exactly the (wrongly-counted) metal bond,
    manufacturing an odd deficit out of nothing. open_shell_atoms now
    recomputes an organic atom's implicit H itself, on top of the
    dative-aware bond total, instead of trusting fill_valence for it."""
    result = validate_smiles_detailed(smiles, closed_shell=True)
    assert result.valid
    assert result.warnings == []
    assert result.repaired_smiles is None


def test_main_group_metal_radical_is_still_flagged():
    # [Na] is a bare main-group atom (not d/f-block) with a genuine unpaired
    # 3s electron; valence_violations already documents it as "accepted
    # under the permissive policy (radicals/under-valent OK)" -- the
    # closed-shell check is expected to name that same radical, consistent
    # with main-group atoms keeping their ordinary valence model elsewhere
    # (docs/validation.md: "Main-group atoms of period 3+ use pysmiles'
    # electron-count maximum").
    result = validate_smiles_detailed("[Na]", closed_shell=True)
    assert result.valid
    assert any("open-shell" in w for w in result.warnings)


# --------------------------------------------------------------------------- #
# Genuine radicals: warned, but stay valid (explicit radical notation)
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("smiles", [
    "[CH3]",                               # methyl radical
    "CN([O])C",                            # nitroxide (aminoxyl radical)
    "c1ccccc1[O]",                         # phenoxyl radical
    "CC1(C)CC(=O)CC(C)(C)N1[O]",           # TEMPO -- a real tox21 spin-label record
    "[N]=O",                               # nitric oxide -- a real clintox/sider record
    "[As]", "[In]",                        # bare main-group atoms -- real sider records
])
def test_genuine_radicals_warned_but_valid(smiles):
    result = validate_smiles_detailed(smiles, closed_shell=True)
    assert result.valid
    assert any("open-shell" in w for w in result.warnings)


def test_even_deficit_is_not_a_forced_radical():
    # A carbene's 2 non-bonding electrons can pair up as a lone pair (singlet):
    # the simple electron-count model used here cannot tell that apart from a
    # triplet biradical, so an even shortfall is deliberately not flagged.
    result = validate_smiles_detailed("[CH2]", closed_shell=True)
    assert result.valid
    assert result.warnings == []


# --------------------------------------------------------------------------- #
# Disconnected free ions
# --------------------------------------------------------------------------- #

def test_free_closed_shell_ions_are_not_radicals():
    result = validate_smiles_detailed("[Na+].[Cl-]", closed_shell=True)
    assert result.valid
    assert result.warnings == []


# --------------------------------------------------------------------------- #
# Default-off equivalence on a sample of the real datasets
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("dataset", ["tox21", "sider", "zinc250k", "moses_test"])
def test_default_off_equivalence_on_dataset_sample(dataset):
    """``closed_shell`` defaults to off: a 300-record sample of each dataset
    must validate identically whether or not the kwarg is passed at all."""
    import pandas as pd

    path = os.path.join(_DATA, f"{dataset}.parquet")
    if not os.path.exists(path):
        pytest.skip(f"{path} not available")
    df = pd.read_parquet(path)
    col = "smiles" if "smiles" in df.columns else "SMILES"
    sample = df[col].dropna().astype(str).str.strip().head(300)

    for smi in sample:
        before = validate_smiles_detailed(smi)
        after = validate_smiles_detailed(smi, closed_shell=False)
        assert before.valid == after.valid
        assert before.warnings == after.warnings
        assert before.repaired_smiles == after.repaired_smiles
        assert before.repairs == after.repairs
        assert before.failed_fragment == after.failed_fragment
        assert before.parent_smiles == after.parent_smiles
        assert before.removed_fragments == after.removed_fragments


# --------------------------------------------------------------------------- #
# Anion classes beyond oxygen (2026-10-20): a stripped negative charge comes
# back as the charge when the X-H acid is deprotonated at pH 7.4, and as a
# hydrogen when it is not.  Each class has a positive case and a negative
# control from the same family that must keep the hydrogen.
# --------------------------------------------------------------------------- #

def _repair(smiles, **kwargs):
    result = validate_smiles_detailed(smiles, closed_shell=True, **kwargs)
    assert result.valid  # the check never changes the verdict
    assert any("open-shell" in w for w in result.warnings)
    return result


ANION_POSITIVES = [
    # (radical as written, expected repair, class word in the repair note)
    ("CCN(CC)C(=S)[S]", "CCN(CC)C(=S)[S-]", "thioacid"),             # dithiocarbamate
    ("CCOC(=S)[S]", "CCOC(=S)[S-]", "thioacid"),                     # xanthate
    ("CC(=O)[S]", "CC(=O)[S-]", "thioacid"),                         # thioacetate
    ("N#C[S]", "N#C[S-]", "thioacid"),                               # thiocyanate
    ("[S]c1ccccc1", "[S-]c1ccccc1", "aromatic thiol"),               # thiophenolate
    ("[S]c1nc2ccccc2s1", "[S-]c1nc2ccccc2s1", "aromatic thiol"),     # 2-mercaptobenzothiazole
    ("CC(=O)[N]S(=O)(=O)c1ccccc1", "CC(=O)[N-]S(=O)(=O)c1ccccc1", "two acceptors"),  # acylsulfonamide
    ("CNC(=O)[N]S(=O)(=O)c1ccc(C)cc1", "CNC(=O)[N-]S(=O)(=O)c1ccc(C)cc1", "two acceptors"),  # sulfonylurea
    ("Nc1ccc(cc1)S(=O)(=O)[N]c1ncccn1", "Nc1ccc(cc1)S(=O)(=O)[N-]c1ncccn1", "two acceptors"),  # sulfadiazine
    ("Cc1ccc(S(=O)(=O)[N]Cl)cc1", "Cc1ccc(S(=O)(=O)[N-]Cl)cc1", "two acceptors"),   # chloramine-T
    ("N#C[N]C#N", "N#C[N-]C#N", "two acceptors"),                    # dicyanamide
    ("CCc1ccc(nc1)CCOc2ccc(cc2)CC3C(=O)[N]C(=O)S3",
     "CCc1ccc(nc1)CCOc2ccc(cc2)CC3C(=O)[N-]C(=O)S3", "thiazolidinedione"),  # pioglitazone
    ("S=C1SCC(=O)[N]1", "S=C1SCC(=O)[N-]1", "thiazolidinedione"),    # rhodanine
    ("CN=[N+]=[N]", "CN=[N+]=[N-]", "azide/diazo"),                  # methyl azide
    ("[N]=[N+]=CC(=O)OC", "[N-]=[N+]=CC(=O)OC", "azide/diazo"),      # diazoacetate
    ("c1cc[n+]([O])cc1", "c1cc[n+]([O-])cc1", "N-oxide"),            # pyridine N-oxide
    ("C[N+](C)(C)[O]", "C[N+](C)(C)[O-]", "N-oxide"),                # trimethylamine N-oxide
    ("CC(C)=[N+]([O-])[O]", "CC(C)=[N+]([O-])[O-]", "N-oxide"),      # nitronate
    ("[O][Cl+3]([O-])([O-])[O-]", "[O-][Cl+3]([O-])([O-])[O-]", "N-oxide"),  # perchlorate, charge-separated
    ("O=Cl(=O)[O]", "O=Cl(=O)[O-]", "halogen oxyacid"),              # chlorate, hypervalent
    ("O=C1OC(C(O)CO)C([O])=C1O", "O=C1OC(C(O)CO)C([O-])=C1O", "vinylogous"),  # ascorbate
    ("[O]c1cc(=O)oc2ccccc12", "[O-]c1cc(=O)oc2ccccc12", "vinylogous"),         # 4-hydroxycoumarin
    ("O=C1CC(C)(C)CC([O])=C1", "O=C1CC(C)(C)CC([O-])=C1", "vinylogous"),       # dimedone enol
    ("O=c1[nH]c([O])cc(=O)[nH]1", "O=c1[nH]c([O-])cc(=O)[nH]1", "vinylogous"), # 6-hydroxyuracil
    ("c1ccn(c(=S)c1)[O]", "c1ccn(c(=S)c1)[O-]", "N-hydroxy-pyridinone"),       # pyrithione
    # ciclopirox: literature pKa 6.8-7.2, the closest call of all these positives (documented)
    ("Cc1cc(n(c(=O)c1)[O])C2CCCCC2", "Cc1cc(n(c(=O)c1)[O-])C2CCCCC2", "N-hydroxy-pyridinone"),
    ("[Cl]", "[Cl-]", "halogen atom"),                               # chloride
    ("[I]", "[I-]", "halogen atom"),                                 # iodide
]


@pytest.mark.parametrize("smiles,repaired,kind", ANION_POSITIVES)
def test_anion_class_repairs_to_negative_charge(smiles, repaired, kind):
    result = _repair(smiles)
    assert result.repaired_smiles == repaired
    assert any("set charge -1" in r and kind in r for r in result.repairs)
    assert validate_smiles_detailed(repaired).valid
    assert Chem.MolFromSmiles(repaired) is not None


ANION_NEGATIVE_CONTROLS = [
    # same families, but the X-H acid is protonated at pH 7.4: H is added
    ("CC[S]", "CCS", "aliphatic thiol"),                             # ethanethiol (pKa 10.6)
    ("OC(=O)C[S]", "OC(=O)CS", "aliphatic thiol"),                   # thioglycolic acid's SH
    ("C[N]S(=O)(=O)c1ccccc1", "CNS(=O)(=O)c1ccccc1", "sulfonamide/imide"),      # N-alkyl sulfonamide
    ("c1ccccc1[N]S(=O)(=O)C", "c1ccccc1NS(=O)(=O)C", "sulfonamide/imide"),      # N-aryl sulfonamide
    ("O=C1[N]C(=O)c2ccccc21", "O=C1NC(=O)c2ccccc21", "sulfonamide/imide"),      # phthalimide (8.3)
    ("CC(=O)[N]C(=O)OC", "CC(=O)NC(=O)OC", "sulfonamide/imide"),        # acyclic N-acyl carbamate
    ("CC(=O)[N]C(=O)SC", "CC(=O)NC(=O)SC", "sulfonamide/imide"),        # acyclic N-acyl thiocarbamate: no TZD ring
    ("O=C1CCSC(=O)[N]1", "O=C1CCSC(=O)N1", "sulfonamide/imide"),        # six-membered analogue: not the TZD ring
    ("CC(=O)[N]C", "CC(=O)NC", "amide"),                             # N-methylacetamide
    ("CN=[N]", "CN=N", "amine/imine"),                               # azo N next to a neutral N
    ("CN([O])C", "CN(O)C", "aliphatic alcohol"),                     # O on a neutral N: N-hydroxy amine
    ("C[Se][O]", "C[Se]O", "aliphatic alcohol"),                     # O on Se with no =O
    ("CC(=O)C=C(C)[O]", "CC(=O)C=C(C)O", "aliphatic alcohol"),       # acyclic 1,3-diketone enol (8.9)
    ("[O]c1ccccc1C=O", "Oc1ccccc1C=O", "phenolic"),                  # salicylaldehyde (H-bonded)
    ("O=C1c2ccccc2C(=O)c2c([O])cccc12", "O=C1c2ccccc2C(=O)c2c(O)cccc12", "phenolic"),  # peri-OH anthraquinone
    ("[O]c1ccc([N+](=O)[O-])cc1", "Oc1ccc([N+](=O)[O-])cc1", "phenolic"),  # 4-nitrophenol: kept neutral
    ("[O]c1c(Cl)cc(Cl)cc1Cl", "Oc1c(Cl)cc(Cl)cc1Cl", "phenolic"),          # trichlorophenol: kept neutral
    ("O=[Se]([O])[O-]", "O=[Se](O)[O-]", "aliphatic alcohol"),       # selenite: second deprotonation (8.3)
]


@pytest.mark.parametrize("smiles,repaired,kind", ANION_NEGATIVE_CONTROLS)
def test_anion_negative_controls_keep_the_hydrogen(smiles, repaired, kind):
    result = _repair(smiles)
    assert result.repaired_smiles == repaired
    assert any("added 1 H" in r and kind in r for r in result.repairs)
    assert Chem.MolFromSmiles(repaired) is not None


@pytest.mark.parametrize("smiles,repaired,_", ANION_POSITIVES)
def test_anion_classes_follow_physiological_ph_off(smiles, repaired, _):
    """Without physiological_ph no pKa table is consulted: H wins whenever it
    closes the shell, for the new classes too."""
    result = _repair(smiles, physiological_ph=False)
    assert result.repaired_smiles != repaired
    assert any("added 1 H" in r and "alternative: set charge -1" in r for r in result.repairs)


@pytest.mark.parametrize("smiles,_,__", ANION_POSITIVES + ANION_NEGATIVE_CONTROLS)
def test_anion_classes_are_off_by_default(smiles, _, __):
    """The new classes only act inside closed_shell=True: off stays silent."""
    result = validate_smiles_detailed(smiles)
    assert result.warnings == [] and result.repaired_smiles is None


def test_bonded_halogen_is_never_open_shell():
    for smiles in ("CCl", "ClC(Cl)Cl", "c1ccccc1Br"):
        result = validate_smiles_detailed(smiles, closed_shell=True)
        assert result.warnings == [] and result.repaired_smiles is None
