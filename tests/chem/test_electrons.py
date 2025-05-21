from chem import Atom, chem


def test_atom():
    assert Atom("C") == Atom("C"), "Atoms with the same symbol should be equal"
    assert Atom("C") != Atom("N"), "Atoms with different symbols should not be equal"
    assert Atom("C") != chem.BracketAtom(
        "C", 0
    ), "Atom and BracketAtom should not be equal"


def test_bracket_atom():
    assert chem.BracketAtom("He", hidrogens=1) == chem.BracketAtom(
        "He", hidrogens=1
    ), "BracketAtoms with the same symbol and hydrogens should be equal"
    assert chem.BracketAtom("He", hidrogens=1) != chem.BracketAtom(
        "He", hidrogens=2
    ), "BracketAtoms with different hydrogens should not be equal"

    assert (
        chem.BracketAtom("He", hidrogens=1).compute_valency() == False
    ), "BracketAtom with 1 hydrogen should not be valency satisfied"
    assert (
        chem.BracketAtom("Li", charge=1).compute_valency() == True
    ), "BracketAtom with charge should be valency satisfied"


def test_valency():
    pass
