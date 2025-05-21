from chem import chem, Atom

def test_atom():
    assert Atom("C") == Atom("C")
    assert Atom("C") != Atom("N")
    assert Atom("C") != chem.BracketAtom("C", 0)
    
def test_bracket_atom():
    assert chem.BracketAtom("He",hidrogens=1) == chem.BracketAtom("He", hidrogens=1)
    assert chem.BracketAtom("He",hidrogens=1) != chem.BracketAtom("He", hidrogens=2)
    
    assert chem.BracketAtom("He", hidrogens=1).compute_valency() == False
    assert chem.BracketAtom("Li", charge=1).compute_valency() == True

def test_valency():
    pass

