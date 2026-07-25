"""Tests for the MolecularGraph -> backend adapter and the pysmiles driver.

These cover three things:
  1. ``to_descriptor`` maps every atom/bond attribute correctly (the hand-off
     contract).
  2. The pysmiles driver returns the right chemistry verdicts.
  3. The backend's SMILES *string* parser is never on the path (bypass guards).
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

import pytest

from src import parse_smiles, validate_smiles
from chem.adapter import to_descriptor, valid_pysmiles


def _descriptor(smiles):
    ok, graph, exc = parse_smiles(smiles)
    assert ok, f"{smiles} failed to parse: {exc}"
    return to_descriptor(graph)


# --------------------------------------------------------------------------- #
# 1. to_descriptor attribute mapping
# --------------------------------------------------------------------------- #

def test_descriptor_charge_and_explicit_h():
    atoms, bonds = _descriptor("CC(=O)[O-]")
    charged = [a for a in atoms if a["charge"] == -1]
    assert len(charged) == 1 and charged[0]["element"] == "O"


def test_descriptor_isotope_and_hcount():
    atoms, _ = _descriptor("[13CH4]")
    assert len(atoms) == 1
    assert atoms[0]["element"] == "C"
    assert atoms[0]["isotope"] == 13
    assert atoms[0]["hcount"] == 4


def test_descriptor_main_chain_bond_order_preserved():
    # Regression guard for the main-chain bond-order fix: O=C=O must carry two
    # double bonds, not be flattened to single bonds.
    _, bonds = _descriptor("O=C=O")
    orders = sorted(o for _, _, o in bonds)
    assert orders == [2, 2]


def test_descriptor_aromatic_bonds_are_1_5():
    _, bonds = _descriptor("c1ccccc1")
    assert len(bonds) == 6
    assert all(o == 1.5 for _, _, o in bonds)


def test_descriptor_dot_is_disconnected():
    atoms, bonds = _descriptor("[Na+].[Cl-]")
    assert len(atoms) == 2
    assert bonds == []


def test_descriptor_branch_connectivity():
    # N must end up bonded to three carbons (degree 3), not a mis-wired hub.
    atoms, bonds = _descriptor("CCN(CC)CC")
    deg = [0] * len(atoms)
    for i, j, _ in bonds:
        deg[i] += 1
        deg[j] += 1
    assert sorted(deg) == [1, 1, 1, 2, 2, 2, 3]


# --------------------------------------------------------------------------- #
# 2. pysmiles driver verdicts
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("smiles,expected", [
    ("c1ccccc1C", True),       # toluene (was wrongly rejected by the regex hack)
    ("CC(=O)[O-]", True),      # acetate
    ("F[P-](F)(F)(F)(F)F", True),  # hexafluorophosphate (octet rule rejected it)
    ("O=S(=O)(O)O", True),     # sulfuric acid (expanded octet)
    ("[Na+].[Cl-]", True),     # salt
    ("[Li+]", True),           # free ion (adapter free-ion rule)
    ("c1ccc2ccccc2c1", True),  # naphthalene
    ("C1=CC=CC=C1", True),     # Kekulé benzene
    ("C(C)(C)(C)(C)C", False), # pentavalent carbon
    ("O=C(=O)=O", False),      # hexavalent carbon
    ("c1cccc1", False),        # 5 aromatic carbons - cannot kekulize
])
def test_driver_verdicts(smiles, expected):
    atoms, bonds = _descriptor(smiles)
    ok, _ = valid_pysmiles(atoms, bonds)
    assert ok is expected, f"{smiles}: expected {expected}, got {ok}"


# --------------------------------------------------------------------------- #
# 3. Bypass guards: the backend's SMILES string parser is never called
# --------------------------------------------------------------------------- #

def test_pysmiles_read_smiles_is_never_called(monkeypatch):
    import pysmiles

    def _boom(*args, **kwargs):
        raise AssertionError("read_smiles must not be called - parser is ours")

    monkeypatch.setattr(pysmiles, "read_smiles", _boom)
    # A full validation still works with the backend string parser disabled.
    ok, _ = validate_smiles("c1ccccc1C")
    assert ok is True


def test_validate_smiles_contract():
    result = validate_smiles("CCO")
    assert isinstance(result, tuple) and len(result) == 2
    assert result[0] is True and result[1] is None
    bad = validate_smiles("C(C)(C)(C)(C)C")
    assert bad[0] is False and bad[1] is not None
