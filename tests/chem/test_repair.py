"""Valence policy, repair suggestions and the salt/pH settings.

The SMILES below are the RDKit-rejected records of public datasets (Tox21,
ClinTox, BBBP and others).  Each is invalid as written but describes a real
molecule; the tests pin both the verdict and the repair.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

import pytest
from rdkit import Chem, RDLogger

from src import validate_smiles_detailed

RDLogger.DisableLog("rdApp.*")


# --------------------------------------------------------------------------- #
# Valence policy
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("smiles,expected", [
    # second-row atoms cannot expand their octet
    ("[NH4]", False),
    ("[NH4]C", False),
    ("[OH3]C", False),
    ("C(N)=[NH2]", False),
    ("O=N([O-])C", False),
    ("[NH4+]", True),
    ("[NH2+]=C(N)N", True),
    # conventional pentavalent N shorthands stay accepted (RDKit accepts them)
    ("CN(=O)=O", True),
    ("c1ccn(=O)cc1", True),
    ("C=N#N", True),
    ("N=N#N", True),
    # closed-shell ions cannot bond
    ("[Cl-]C", False),
    ("[F-][B+3]([F-])([F-])[F-]", False),
    ("[Na+].[Cl-]", True),
    # d/f-block metals: no valence limit; neutral N/O donors bond datively
    ("[NH3][Pt]([NH3])(Cl)Cl", True),
    ("Cl[Pt]1(Cl)[NH2]CC[NH2]1", True),
    ("[NH4+].[NH4+].Cl[Pt-2](Cl)(Cl)Cl", True),
    ("Cl[Cu]Cl", True),
    ("[NH4][Pt]([NH4])(Cl)Cl", False),
    # main-group metals keep their valence
    ("O[AlH3](O)O", False),
    ("O[Al](O)O", True),
    # aromatic atoms: one pi bond, not 1.5 per aromatic bond
    ("c1ccc2ccccc2c1", True),
    ("c1cc[nH+]cc1", True),
])
def test_valence_policy(smiles, expected):
    result = validate_smiles_detailed(smiles)
    assert result.valid is expected, result.message


# --------------------------------------------------------------------------- #
# cis/trans marks on a non-single bond
# --------------------------------------------------------------------------- #

DIAZO = [
    r"NNC(=O)CNC(=O)\C=N\#N",
    r"O=C1NC(=O)\C(=N/#N)\C=N1",
    r"NC(=O)CNC(=O)\C=N\#N",
    r"CCCCN(CC(O)C1=C\C(=N/#N)\C(=O)C=C1)N=O",
    r"NC(COC(=O)\C=N/#N)C(=O)O",
    r"CCN(CC(O)C1=CC(=O)\C(=N\#N)\C=C1)N=O",
]


@pytest.mark.parametrize("smiles", DIAZO)
def test_direction_mark_on_triple_bond_is_dropped(smiles):
    result = validate_smiles_detailed(smiles)
    assert result.valid, result.message
    assert any("cis/trans mark" in w for w in result.warnings)


@pytest.mark.parametrize("smiles", [r"C\-C", r"c1cccc/:c1", r"C/=1CCCCC1"])
def test_direction_mark_before_any_bond_symbol(smiles):
    assert validate_smiles_detailed(smiles).valid


def test_plain_direction_marks_give_no_warning():
    result = validate_smiles_detailed("F/C=C/F")
    assert result.valid and result.warnings == []


def test_doubled_backslash_reads_as_one():
    # ClinTox stores some SMILES with escaped backslashes.
    result = validate_smiles_detailed(r"C[C@H](/C=C/C)CC/C=C\\C=C/C")
    assert result.valid, result.message
    assert any("doubled backslash" in w for w in result.warnings)


def test_two_direction_marks_still_rejected():
    assert validate_smiles_detailed("C//C").tier == "grammar"


# --------------------------------------------------------------------------- #
# Repairs
# --------------------------------------------------------------------------- #

# (input, repaired at pH 7.4, repaired without pH)
REPAIRS = [
    ("NC(=O)NC1N=C(O[AlH3](O)O)NC1=O", "NC(=O)NC1N=C(O[Al](O)O)NC1=O", None),
    ("CCCCO[AlH3](OCCCC)OCCCC", "CCCCO[Al](OCCCC)OCCCC", None),
    ("[NH4][Pt]([NH4])(Cl)Cl", "[NH3][Pt]([NH3])(Cl)Cl", None),
    ("O=N([O-])C1=C(CN=C1NCCSCc2ncccc2)Cc3ccccc3",
     "O=[N+]([O-])C1=C(CN=C1NCCSCc2ncccc2)Cc3ccccc3", None),
    ("c1(nc(NC(N)=[NH2])sc1)CSCCNC(=[NH]C#N)NC",
     "c1(nc(NC(N)=[NH2+])sc1)CSCCNC(=NC#N)NC",
     "c1(nc(NC(N)=N)sc1)CSCCNC(=NC#N)NC"),
    (r"Cc1nc(sc1)\[NH]=C(\N)N", r"Cc1nc(sc1)\[NH+]=C(\N)N", r"Cc1nc(sc1)\N=C(\N)N"),
    (r"c1c(c(ncc1)CSCCN\C(=[NH]\C#N)NCC)Br", r"c1c(c(ncc1)CSCCN\C(=N\C#N)NCC)Br", None),
    (r"c1(cc(N\C(=[NH]\c2cccc(c2)CC)C)ccc1)CC",
     r"c1(cc(N\C(=[NH+]\c2cccc(c2)CC)C)ccc1)CC",
     r"c1(cc(N\C(=N\c2cccc(c2)CC)C)ccc1)CC"),
    ("[F-][B+3]([F-])([F-])[F-].CC[N+]1(C)CCCC1", "F[B-](F)(F)F.CC[N+]1(C)CCCC1", None),
    ("[NH4+].[NH4+].[Cl-][Pt++]([Cl-])([Cl-])[Cl-]", "[NH4+].[NH4+].Cl[Pt-2](Cl)(Cl)Cl", None),
    ("[Cl-][Pt]1([Cl-])[NH2+]CC[NH2+]1", "Cl[Pt-2]1(Cl)[NH2+]CC[NH2+]1", None),
    ("[NH4]", "[NH4+]", "N"),
]


@pytest.mark.parametrize("smiles,at_ph,without_ph", REPAIRS)
def test_repair(smiles, at_ph, without_ph):
    result = validate_smiles_detailed(smiles)
    assert not result.valid and result.tier == "chemistry"
    assert result.repaired_smiles == at_ph
    assert result.repairs
    # the repair is valid for us and for RDKit
    assert validate_smiles_detailed(at_ph).valid
    assert Chem.MolFromSmiles(at_ph) is not None

    plain = validate_smiles_detailed(smiles, physiological_ph=False)
    assert plain.repaired_smiles == (without_ph or at_ph)


def test_repair_keeps_stereo():
    smiles = "O=S(=O)(OC[C@H]1O[C@H](O[C@]2(COS(=O)(=O)O[AlH3](O)O)OC2)CC1)O"
    result = validate_smiles_detailed(smiles)
    assert result.repaired_smiles == smiles.replace("[AlH3]", "[Al]")


def test_ph_explanation_is_reported():
    result = validate_smiles_detailed("c1(nc(NC(N)=[NH2])sc1)C")
    assert any("guanidine" in line and "pH 7.4" in line for line in result.repairs)


@pytest.mark.parametrize("smiles", [
    "CCCCc1c(=O)n(n(c1=O)c2ccccc2)c3ccccc3",  # aromatic ring that cannot exist
    "C(C)(C)(C)(C)C",                          # pentavalent carbon
])
def test_no_repair_when_not_a_single_atom_fix(smiles):
    result = validate_smiles_detailed(smiles)
    assert not result.valid
    assert result.repaired_smiles is None


# --------------------------------------------------------------------------- #
# Salt removal
# --------------------------------------------------------------------------- #

def test_salt_parent_and_counterions():
    result = validate_smiles_detailed("CC(=O)[O-].[Na+]")
    assert result.valid
    assert result.parent_smiles == "CC(=O)[O-]"
    assert result.removed_fragments == ["[Na+]"]


def test_salt_removal_off():
    result = validate_smiles_detailed("CC(=O)[O-].[Na+]", remove_salts=False)
    assert result.parent_smiles is None and result.removed_fragments == []


def test_broken_counterion_is_reported_not_hidden():
    result = validate_smiles_detailed("[F-][B+3]([F-])([F-])[F-].CC[N+]1(C)CCCC1")
    assert not result.valid
    assert result.failed_fragment == "[F-][B+3]([F-])([F-])[F-]"
    assert result.message.startswith("counter-ion")
    assert result.parent_smiles == "CC[N+]1(C)CCCC1"
    assert result.removed_fragments == ["F[B-](F)(F)F"]


def test_single_component_has_no_parent():
    assert validate_smiles_detailed("CCO").parent_smiles is None
