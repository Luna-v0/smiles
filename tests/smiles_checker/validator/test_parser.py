import pytest

from smiles_checker.chem.atomic import Atom, BracketAtom
from smiles_checker.chem.chemistry import chemistry as chem
from smiles_checker.chem.chemistry import chemistry as chem
from smiles_checker.chem.chemistry import chemistry as chem
from smiles_checker.chem.chemistry import chemistry as chem
from smiles_checker.exceptions import ParserException
from smiles_checker.validator.parser_manager import ParserManager


@pytest.fixture(scope="function")
def parser_manager():
    """
    Fixture to create a ParserManager instance for each test.
    """
    pm = ParserManager()
    pm._reset()  # Ensure the parser manager is reset before each test
    return pm


def test_validate_branch(parser_manager: ParserManager):
    """
    Test the validate_branch method of the ParserManager.
    """
    # Test with no open ring numbers
    assert (
        parser_manager.validate_branch() is True
    ), "Should return True when no open ring numbers"


def test_internal_bracket(parser_manager: ParserManager):
    """ """
    assert parser_manager.internal_bracket(symbol="C") == chem.BracketAtom(symbol="C"), "Internal bracket should return BracketAtom('C')"
    assert (
        len(parser_manager.current_chain) == 1
    ), "Current chain should have one atom after internal bracket"


def test_listify(parser_manager: ParserManager):
    assert parser_manager.listify(base_element="C") == [
        "C"
    ], "Listify should return a list with one element 'C'"
    assert parser_manager.listify(base_element="C", recursion=["O"]) == [
        "C",
        "O",
    ], "Listify should combine base_element with a list recursion"
    assert parser_manager.listify(base_element="C", recursion="O") == [
        "C",
        "O",
    ], "Listify should combine base_element with a non-list recursion"


def test_atom(parser_manager: ParserManager):
    """
    Test the atom parser with valid and invalid inputs.
    """

    assert parser_manager.atom("C").symbol == Atom("C").symbol, "Parsing 'C' should return Atom('C')"
    assert parser_manager.atom("Na").symbol == Atom("Na").symbol, "Parsing 'Na' should return Atom('Na')"

    ## TODO missing tests with bracket atom

    with pytest.raises(ParserException) as exc_info:
        parser_manager.atom("X")
    assert exc_info.value.message == "Invalid Atom Symbol: X" and exc_info.value.rule == "Atom"

    with pytest.raises(ParserException) as exc_info:
        parser_manager.atom("Xx")
    assert exc_info.value.message == "Invalid Atom Symbol: Xx" and exc_info.value.rule == "Atom"


def test_fifteen(parser_manager: ParserManager):
    """
    Test the parser manager with a simple SMILES string.
    """
    assert parser_manager.fifteen("1", "5") == 15, "Parsing 1 and 5 should return 15"
    # assert error is raised when exceeding 15
    with pytest.raises(ParserException) as exc_info:
        parser_manager.fifteen("1", "6")
        assert exc_info.value.message == "Cannot exceed 15" and exc_info.value.rule == "fifteen"


def test_chiral(parser_manager: ParserManager):
    """
    Test the chiral parser.
    """
    assert (
        parser_manager.chiral(chiral1="@", chiral2=None) == "clockwise"
    ), "Chiral without second symbol should be clockwise"
    assert (
        parser_manager.chiral(chiral1="@", chiral2="@") == "counterclockwise"
    ), "Chiral with second symbol should be counterclockwise"


def test_charge(parser_manager: ParserManager):
    """
    Test the charge parser.
    """
    assert parser_manager.charge("+", None) == 1, "Charge '+' should return 1"
    assert parser_manager.charge("-", None) == -1, "Charge '-' should return -1"
    assert (
        parser_manager.charge("+", "+") == 2
    ), "Charge '+' with second '+' should return 2"
    assert (
        parser_manager.charge("-", "-") == -2
    ), "Charge '-' with second '-' should return -2"
    assert (
        parser_manager.charge("-", 3) == -3
    ), "Charge '-' with integer should return negative integer"
    assert (
        parser_manager.charge("+", 3) == 3
    ), "Charge '+' with integer should return positive integer"

    with pytest.raises(ParserException) as exc_info:
        parser_manager.charge("+", "-")
        assert exc_info.value.message == "Charge mismatch" and exc_info.value.rule == "charge"


def test_hcount(parser_manager: ParserManager):
    """
    Test the hydrogen count parser.
    """
    assert parser_manager.hcount(_="H", digit="1") == 1, "Hydrogen count '1' should return 1"
    assert parser_manager.hcount(_="H", digit="2") == 2, "Hydrogen count '2' should return 2"
    assert parser_manager.hcount(_="H", digit="3") == 3, "Hydrogen count '3' should return 3"

    with pytest.raises(ParserException) as exc_info:
        parser_manager.hcount(_="H", digit="A")
    assert exc_info.value.message == "Invalid hydrogen count" and exc_info.value.rule == "hcount"


def test_ring_number(parser_manager: ParserManager):
    """
    Test the ring number parser.
    """
    assert parser_manager.ring_number(ring_number_or_symbol="1", ring_number1=None, ring_number2=None) == 1, "Ring number 1 opened"
    assert parser_manager.ring_number(ring_number_or_symbol="2", ring_number1=None, ring_number2=None) == 2, "Ring number 2 opened"
    assert parser_manager.ring_number(ring_number_or_symbol="3", ring_number1=None, ring_number2=None) == 3, "Ring number 3 opened"

    assert parser_manager.ring_number(ring_number_or_symbol="1", ring_number1=None, ring_number2=None) == 1, "Ring number 1 closed"
    assert parser_manager.ring_number(ring_number_or_symbol="2", ring_number1=None, ring_number2=None) == 2, "Ring number 2 closed"
    assert parser_manager.ring_number(ring_number_or_symbol="3", ring_number1=None, ring_number2=None) == 3, "Ring number 3 closed"

    assert parser_manager.ring_number(ring_number_or_symbol="4", ring_number1=None, ring_number2=None) == 4, "Ring number 4 opened"
    assert parser_manager.ring_number(ring_number_or_symbol="4", ring_number1=None, ring_number2=None) == 4, "Ring number 4 closed"

    assert parser_manager.ring_number(ring_number_or_symbol="%", ring_number1="5", ring_number2=None) == 5, "Ring number 5 opened"
    assert parser_manager.ring_number(ring_number_or_symbol="%", ring_number1="5", ring_number2=None) == 5, "Ring number 5 closed"

    with pytest.raises(ParserException) as exc_info:
        parser_manager.ring_number(ring_number_or_symbol="A", ring_number1=None, ring_number2=None)
    assert exc_info.value.message == "Ring number must be a digit or a number with a leading digit"

    with pytest.raises(ParserException) as exc_info:
        parser_manager.ring_number(ring_number_or_symbol="1", ring_number1=None, ring_number2=None)
    assert exc_info.value.message == "Ring number already closed"

    with pytest.raises(ParserException) as exc_info:
        parser_manager.ring_number(ring_number_or_symbol="%", ring_number1=None, ring_number2=None)
    assert exc_info.value.message == "Ring number cannot be just '%'"

    with pytest.raises(ParserException) as exc_info:
        parser_manager.ring_number(ring_number_or_symbol="1", ring_number1="2", ring_number2="3")
    assert exc_info.value.message == "Ring number cannot have more than one digit after the first" and exc_info.value.rule == "ring_number"
