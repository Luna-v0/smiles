import os
import sys

import pytest

from chem import chemistry
from exceptions import ParserException

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))
from src import validate_smiles


# --------------------------------------------------------------------------- #
# Atom factory + symbol validation + graph identity
# --------------------------------------------------------------------------- #

def test_atom():
    # Atoms are compared by atom_id for graph purposes, not symbol: each Atom
    # instance is unique even with the same symbol.
    atom1 = chemistry.Atom("C")
    atom2 = chemistry.Atom("C")
    assert atom1 != atom2, "Different Atom instances should not be equal (graph uniqueness)"
    assert atom1.symbol == atom2.symbol, "Atoms with the same symbol should have same symbol attribute"
    assert chemistry.Atom("C").symbol != chemistry.Atom("N").symbol
    assert chemistry.Atom("C") != chemistry.BracketAtom("C", hidrogens=0), \
        "Atom and BracketAtom should not be equal"


def test_bracket_atom():
    # BracketAtoms are also compared by atom_id, and carry their fields.
    atom1 = chemistry.BracketAtom("He", hidrogens=1)
    atom2 = chemistry.BracketAtom("He", hidrogens=1)
    assert atom1 != atom2, "Different BracketAtom instances should not be equal (graph uniqueness)"
    assert atom1.symbol == atom2.symbol and atom1.hcount == atom2.hcount
    atom3 = chemistry.BracketAtom("He", hidrogens=2)
    assert atom1.hcount != atom3.hcount


def test_atom_invalid_symbol():
    with pytest.raises(ParserException) as exc_info:
        chemistry.Atom("X")
    assert exc_info.value.message == "Invalid Atom Symbol: X"


def test_bracket_atom_invalid_symbol():
    with pytest.raises(ParserException) as exc_info:
        chemistry.BracketAtom("Xx")
    assert exc_info.value.message == "Invalid Atom Symbol: Xx"


def test_atom_lowercase_symbol():
    atom = chemistry.Atom("c")
    assert atom.symbol == "C"


def test_bracket_atom_lowercase_symbol():
    atom = chemistry.BracketAtom("c")
    assert atom.symbol == "C"
    assert atom.aromatic is True


# --------------------------------------------------------------------------- #
# Valence verdicts via the live pipeline (pysmiles backend, permissive policy).
# These replaced the old octet/duet unit tests; see
# docs/chemistry_backend_study.md.
# --------------------------------------------------------------------------- #

def test_bracket_atom_valency_li_alkyl_valid():
    """CC[Li] (ethyllithium) is valid: a one-bonded Li satisfies valence 1."""
    result, _ = validate_smiles("CC[Li]")
    assert result is True


def test_bracket_atom_valency_li_cation_stable():
    """[Li+] is valid - Li+ has a duet configuration."""
    result, _ = validate_smiles("[Li+]")
    assert result is True


def test_bracket_atom_valency_na_valid():
    """[Na] is accepted under the permissive policy (radicals/under-valent OK)."""
    result, _ = validate_smiles("[Na]")
    assert result is True


def test_bracket_atom_valency_na_cation_stable():
    """[Na+] is valid - Na+ has a noble-gas configuration."""
    result, _ = validate_smiles("[Na+]")
    assert result is True


def test_bracket_atom_valency_ch4_stable():
    """[CH4] is valid - carbon with four hydrogens."""
    result, _ = validate_smiles("[CH4]")
    assert result is True


def test_bracket_atom_valency_nh3_stable():
    """[NH3] is valid - nitrogen with three hydrogens."""
    result, _ = validate_smiles("[NH3]")
    assert result is True


def test_bracket_atom_valency_oh2_stable():
    """[OH2] is valid - oxygen with two hydrogens."""
    result, _ = validate_smiles("[OH2]")
    assert result is True


def test_bracket_atom_valency_alh3_valid():
    """[AlH3] (alane) is valid: aluminium's valence is 3."""
    result, _ = validate_smiles("[AlH3]")
    assert result is True
