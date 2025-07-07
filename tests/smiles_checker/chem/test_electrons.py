from smiles_checker.chem import Atom, chemistry


def test_atom():
    assert Atom("C") == Atom("C"), "Atoms with the same symbol should be equal"
    assert Atom("C") != Atom("N"), "Atoms with different symbols should not be equal"
    assert Atom("C") != chemistry.BracketAtom(
        "C", 0
    ), "Atom and BracketAtom should not be equal"


def test_bracket_atom():
    assert chemistry.BracketAtom("He", hidrogens=1) == chemistry.BracketAtom(
        "He", hidrogens=1
    ), "BracketAtoms with the same symbol and hydrogens should be equal"
    assert chemistry.BracketAtom("He", hidrogens=1) != chemistry.BracketAtom(
        "He", hidrogens=2
    ), "BracketAtoms with different hydrogens should not be equal"

    assert (
        chemistry.BracketAtom("He", hidrogens=1).compute_valency() == False
    ), "BracketAtom with 1 hydrogen should not be valency satisfied"
    assert (
        chemistry.BracketAtom("Li", charge=1).compute_valency() == True
    ), "BracketAtom with charge should be valency satisfied"


def test_valency():
    pass


def test_electron_distribution_and_subshells():
    # Test adding electrons (negative charge)
    # Oxygen normally: 1s2 2s2 2p4 (6 valence electrons)
    # O-2 should be: 1s2 2s2 2p6 (8 valence electrons)
    oxygen_anion = chemistry.BracketAtom("O", charge=-2)
    assert oxygen_anion.get_total_electrons_in_subshell("p") == 6, "Oxygen anion should have 6 electrons in p subshells"
    assert oxygen_anion.electrons_in_valency == 8, "Oxygen anion should have 8 valence electrons"

    # Test removing electrons (positive charge)
    # Carbon normally: 1s2 2s2 2p2 (4 valence electrons)
    # C+2 should be: 1s2 2s2 (2 valence electrons)
    carbon_cation = chemistry.BracketAtom("C", charge=2)
    assert carbon_cation.get_total_electrons_in_subshell("p") == 0, "Carbon cation should have 0 electrons in p subshells"
    assert carbon_cation.electrons_in_valency == 2, "Carbon cation should have 2 valence electrons"

    # Test with hydrogens (effectively adding electrons to valency)
    # Nitrogen normally: 1s2 2s2 2p3 (5 valence electrons)
    # N with 3 hydrogens (NH3) should be stable, effectively 8 valence electrons
    nitrogen_with_hydrogens = chemistry.BracketAtom("N", hidrogens=3)
    assert nitrogen_with_hydrogens.electrons_in_valency == 8, "Nitrogen with 3 hydrogens should have 8 valence electrons"

    # Test a more complex scenario: Sulfur (1s2 2s2 2p6 3s2 3p4) with -2 charge
    # Should become 1s2 2s2 2p6 3s2 3p6
    sulfur_anion = chemistry.BracketAtom("S", charge=-2)
    assert sulfur_anion.get_total_electrons_in_subshell("p") == 6, "Sulfur anion should have 6 electrons in 3p subshell"
    assert sulfur_anion.electrons_in_valency == 8, "Sulfur anion should have 8 valence electrons"

    # Test a scenario where electrons are removed from inner shells if outermost are empty
    # Lithium (1s2 2s1) with +1 charge should be 1s2
    lithium_cation = chemistry.BracketAtom("Li", charge=1)
    assert lithium_cation.get_total_electrons_in_subshell("s") == 2, "Lithium cation should have 2 electrons in 1s subshell"
    assert lithium_cation.electrons_in_valency == 2, "Lithium cation should have 2 valence electrons (1s2)"

    # Test a scenario where electrons are added to new layers if existing are full
    # Neon (1s2 2s2 2p6) with -1 charge (hypothetical)
    # Should become 1s2 2s2 2p6 3s1
    neon_anion = chemistry.BracketAtom("Ne", charge=-1)
    assert neon_anion.get_total_electrons_in_subshell("s") == 3, "Neon anion should have 3 electrons in s subshells (2s2 + 3s1)"
    assert neon_anion.electrons_in_valency == 1, "Neon anion should have 1 valence electron (in 3s)"
    assert neon_anion.valency_layer == 3, "Neon anion valency layer should be 3"
