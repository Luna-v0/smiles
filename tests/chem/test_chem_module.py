import pytest
from chem import chem


def test_number_of_electrons_per_bond():
    assert chem.number_of_electrons_per_bond("=") == 2
    assert chem.number_of_electrons_per_bond("#") == 3
    with pytest.raises(Exception):
        chem.number_of_electrons_per_bond("?")


def test_atom_creation_and_aromatic_flag():
    atom = chem.Atom("C")
    assert atom.symbol == "C"
    assert atom.aromatic is False
    with pytest.raises(Exception):
        chem.Atom("Zz")
    aromatic_atom = chem.Atom("c")
    assert aromatic_atom.symbol == "C"


def test_bracket_atom_and_validation():
    br = chem.BracketAtom("C", hidrogens=0, charge=0)
    assert br.symbol == "C"
    assert isinstance(br.compute_valency(), bool)
    assert chem.validate_valency_bracket(None, "He", None, None, None, None)
