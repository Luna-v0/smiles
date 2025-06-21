from chem import chem


def test_atom_equality():
    a = chem.Atom("C")
    b = chem.Atom("C")
    c = chem.Atom("N")
    assert a == b
    assert a != c


def test_bracket_atom_equality():
    b1 = chem.BracketAtom("He", hidrogens=1)
    b2 = chem.BracketAtom("He", hidrogens=1)
    b3 = chem.BracketAtom("He", hidrogens=2)
    assert b1 == b2
    assert b1 != b3
