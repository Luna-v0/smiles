from chem import chem, Atom, BracketAtom

def test_atom():
    assert Atom("C") == Atom("C")
    assert Atom("C") != Atom("N")

    assert Atom("C") != BracketAtom("C", 1)


def test_valency():
    pass

