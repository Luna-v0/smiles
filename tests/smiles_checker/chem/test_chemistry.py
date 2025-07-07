from smiles_checker.chem import Atom, chemistry
import pytest

def test_atom():
    assert chemistry.Atom("C") == chemistry.Atom(
        "C"
    ), "Atoms with the same symbol should be equal"
    assert chemistry.Atom("C") != chemistry.Atom(
        "N"
    ), "Atoms with different symbols should not be equal"
    assert chemistry.Atom("C") != chemistry.BracketAtom(
        "C", hidrogens=0
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
    assert (
        oxygen_anion.get_total_electrons_in_subshell("p") == 6
    ), "Oxygen anion should have 6 electrons in p subshells"
    assert (
        oxygen_anion.electrons_in_valency == 8
    ), "Oxygen anion should have 8 valence electrons"

    # Test removing electrons (positive charge)
    # Carbon normally: 1s2 2s2 2p2 (4 valence electrons)
    # C+2 should be: 1s2 2s2 (2 valence electrons)
    carbon_cation = chemistry.BracketAtom("C", charge=2)
    assert (
        carbon_cation.get_total_electrons_in_subshell("p") == 0
    ), "Carbon cation should have 0 electrons in p subshells"
    assert (
        carbon_cation.electrons_in_valency == 2
    ), "Carbon cation should have 2 valence electrons"

    # Test with hydrogens (effectively adding electrons to valency)
    # Nitrogen normally: 1s2 2s2 2p3 (5 valence electrons)
    # N with 3 hydrogens (NH3) should be stable, effectively 8 valence electrons
    nitrogen_with_hydrogens = chemistry.BracketAtom("N", hidrogens=3)
    assert (
        nitrogen_with_hydrogens.compute_valency() == True
    ), "Nitrogen with 3 hydrogens should be stable"

    # Test a more complex scenario: Sulfur (1s2 2s2 2p6 3s2 3p4) with -2 charge
    # Should become 1s2 2s2 2p6 3s2 3p6
    sulfur_anion = chemistry.BracketAtom("S", charge=-2)
    assert (
        sulfur_anion.get_electrons_in_specific_subshell(3, "p") == 6
    ), "Sulfur anion should have 6 electrons in 3p subshell"
    assert (
        sulfur_anion.electrons_in_valency == 8
    ), "Sulfur anion should have 8 valence electrons"

    # Test a scenario where electrons are removed from inner shells if outermost are empty
    # Lithium (1s2 2s1) with +1 charge should be 1s2
    lithium_cation = chemistry.BracketAtom("Li", charge=1)
    assert (
        lithium_cation.get_total_electrons_in_subshell("s") == 2
    ), "Lithium cation should have 2 electrons in 1s subshell"
    assert (
        lithium_cation.electrons_in_valency == 2
    ), "Lithium cation should have 2 valence electrons (1s2)"

    # Test a scenario where electrons are added to new layers if existing are full
    # Neon (1s2 2s2 2p6) with -1 charge (hypothetical)
    # Should become 1s2 2s2 2p6 3s1
    neon_anion = chemistry.BracketAtom("Ne", charge=-1)
    assert (
            neon_anion.get_electrons_in_specific_subshell(3, "s") == 1
        ), "Neon anion should have 1 electron in 3s subshell"
    assert (
        neon_anion.electrons_in_valency == 1
    ), "Neon anion should have 1 valence electron (in 3s)"
    assert neon_anion.valency_layer == 3, "Neon anion valency layer should be 3"

def test_number_of_electrons_per_bond():
    assert chemistry.number_of_electrons_per_bond("=") == 2
    assert chemistry.number_of_electrons_per_bond("#") == 3
    assert chemistry.number_of_electrons_per_bond("$") == 4
    assert chemistry.number_of_electrons_per_bond("/") == 1
    assert chemistry.number_of_electrons_per_bond("\\") == 1
    assert chemistry.number_of_electrons_per_bond("-") == 1
    assert chemistry.number_of_electrons_per_bond(".") == 0
    with pytest.raises(Exception, match="Invalid Bond invalid_bond"):
        chemistry.number_of_electrons_per_bond("invalid_bond")

def test_atom_invalid_symbol():
    with pytest.raises(Exception, match="Invalid Symbol X"):
        chemistry.Atom("X")

def test_bracket_atom_invalid_symbol():
    with pytest.raises(Exception, match="Invalid Symbol Xx"):
        chemistry.BracketAtom("Xx")

def test_atom_empty_electron_configuration():
    atom = Atom("H", electron_configuration="")
    assert atom.valency_layer == 0
    assert atom.electrons_in_valency == 0
    assert atom.layers == ()
    assert atom.electrons_by_layers == ()

def test_atom_no_match_in_electron_configuration():
    atom = Atom("H", electron_configuration="invalid_config")
    assert atom.valency_layer == 0
    assert atom.electrons_in_valency == 0
    assert atom.layers == ()
    assert atom.electrons_by_layers == ()

def test_bracket_atom_electron_addition_new_shell():
    # Test adding electrons to create a new shell (e.g., H- to H)
    # H normally: 1s1
    # H- should be: 1s2
    hydrogen_anion = chemistry.BracketAtom("H", charge=-1)
    assert hydrogen_anion.get_electrons_in_specific_subshell(1, "s") == 2
    assert hydrogen_anion.electrons_in_valency == 2
    assert hydrogen_anion.valency_layer == 1

def test_bracket_atom_electron_removal_empty_shell():
    # Test removing all electrons from an atom
    # H normally: 1s1
    # H+ should be: no electrons
    hydrogen_cation = chemistry.BracketAtom("H", charge=1)
    assert hydrogen_cation.electrons_in_valency == 0
    assert hydrogen_cation.valency_layer == 0
    assert hydrogen_cation.layers == ()
    assert hydrogen_cation.electrons_by_layers == ()

def test_next_subshell_index_error():
    # This test is designed to hit the IndexError in _next_subshell
    # by passing 'f' and then trying to get the next subshell.
    # This will force the 'return 's'' branch.
    atom = chemistry.BracketAtom("H")
    assert atom._next_subshell('f') == 's'

def test_max_electrons_in_subshell_invalid():
    atom = chemistry.BracketAtom("H")
    assert atom._max_electrons_in_subshell('x') == 0

def test_atom_lowercase_symbol():
    atom = chemistry.Atom("c")
    assert atom.symbol == "C"

def test_bracket_atom_lowercase_symbol():
    atom = chemistry.BracketAtom("c")
    assert atom.symbol == "C"
    assert atom.aromatic == True

def test_validate_valency_bracket():
    assert chemistry.validate_valency_bracket(isotope=None, symbol="C", chiral=None, hcount=4, charge=0, map=None) == True
    assert chemistry.validate_valency_bracket(isotope=None, symbol="O", chiral=None, hcount=0, charge=0, map=None) == False

def test_get_electrons_in_specific_subshell_not_found():
    atom = Atom("C")
    assert atom.get_electrons_in_specific_subshell(1, "p") == 0

def test_bracket_atom_electron_addition_empty_atom():
    atom = chemistry.BracketAtom("H", charge=-1)
    assert atom.electrons_in_valency == 2

def test_bracket_atom_electron_addition_multiple_shells():
    atom = chemistry.BracketAtom("Ne", charge=-1)
    assert atom.get_electrons_in_specific_subshell(3, "s") == 1

def test_validate_valency_bracket_no_hcount_no_charge():
    assert chemistry.validate_valency_bracket(isotope=None, symbol="C", chiral=None, hcount=None, charge=None, map=None) == False

def test_validate_valency_bracket_with_charge():
    assert chemistry.validate_valency_bracket(isotope=None, symbol="N", chiral=None, hcount=0, charge=1, map=None) == False

def test_validate_valency_bracket_with_hcount():
    assert chemistry.validate_valency_bracket(isotope=None, symbol="O", chiral=None, hcount=2, charge=0, map=None) == True

def test_atom_next_subshell_invalid_input():
    atom = Atom("H")
    assert atom._next_subshell('z') == 's'