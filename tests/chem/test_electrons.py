from chem import chem, Atom, BracketAtom

def test_atom():
    assert Atom("C") == Atom("C")
    assert Atom("C") != Atom("N")
    assert Atom("C") != BracketAtom("C", 0)
    
def test_bracket_atom():
    assert BracketAtom("He",hidrogens=1) == BracketAtom("He", hidrogens=1)
    assert BracketAtom("He",hidrogens=1) != BracketAtom("He", hidrogens=2)
    
    assert BracketAtom("He", hidrogens=1).compute_valency() == False
    assert BracketAtom("Li", charge=1).compute_valency() == True

def test_valency():
    pass

