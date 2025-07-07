import pytest

from smiles_checker.chem import Atom, BracketAtom
from smiles_checker.validator.parser_manager import ParserException, ParserManager


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
    assert parser_manager.internal_bracket(None, "C") == BracketAtom(
        "C"
    ), "Internal bracket should return BracketAtom('C', 'O')"
    assert (
        len(parser_manager.current_chain) == 1
    ), "Current chain should have one atom after internal bracket"


def test_listify(parser_manager: ParserManager):
    assert parser_manager.listify("C") == [
        "C"
    ], "Listify should return a list with one element 'C'"
    assert parser_manager.listify(["C", "O"]) == [
        "O",
        "C",
    ], "Listify should return a list with elements 'O' and 'C' in that order"
    assert parser_manager.listify(["C", "O", "N"]) == [
        "N",
        "O",
        "C",
    ], "Listify should return a list with elements 'N', 'O', and 'C' in that order"


def test_atom(parser_manager: ParserManager):
    """
    Test the atom parser with valid and invalid inputs.
    """

    assert parser_manager.atom("C") == Atom("C"), "Parsing 'C' should return Atom('C')"
    assert parser_manager.atom("Na") == Atom(
        "Na"
    ), "Parsing 'Na' should return Atom('Na')"

    ## TODO missing tests with bracket atom

    with pytest.raises(ParserException):
        parser_manager.atom("X")

    with pytest.raises(ParserException):
        parser_manager.atom("au")


def test_fifteen(parser_manager: ParserManager):
    """
    Test the parser manager with a simple SMILES string.
    """
    assert parser_manager.fifteen("1", "5") == 15, "Parsing 1 and 5 should return 15"
    # assert error is raised when exceeding 15
    with pytest.raises(ParserException) as exc_info:
        parser_manager.fifteen("1", "6")
        assert (
            ParserException(rule="fifteen", parameter="1 6", message="Cannot exceed 15")
            == exc_info.value
        ), "Should raise ParserException when exceeding 15"


def test_chiral(parser_manager: ParserManager):
    """
    Test the chiral parser.
    """
    assert (
        parser_manager.chiral("@", None) == "clockwise"
    ), "Chiral without second symbol should be clockwise"
    assert (
        parser_manager.chiral("@@", None) == "counterclockwise"
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
        assert (
            ParserException(rule="charge", parameter="+ -", message="Charge mismatch")
            == exc_info.value
        ), "Should raise ParserException for charge mismatch"


def test_hcount(parser_manager: ParserManager):
    """
    Test the hydrogen count parser.
    """
    assert parser_manager.hcount("1", None) == 1, "Hydrogen count '1' should return 1"
    assert parser_manager.hcount("2", None) == 2, "Hydrogen count '2' should return 2"
    assert parser_manager.hcount("3", None) == 3, "Hydrogen count '3' should return 3"

    with pytest.raises(ParserException) as exc_info:
        parser_manager.hcount("A", None)
        assert (
            ParserException(
                rule="hcount", parameter="A None", message="Invalid hydrogen count"
            )
            == exc_info.value
        ), "Should raise ParserException for invalid hydrogen count"


def test_ring_number(parser_manager: ParserManager):
    """
    Test the ring number parser.
    """
    assert parser_manager.ring_number("1", None) == 1, "Ring number 1 opened"
    assert parser_manager.ring_number("2", None) == 2, "Ring number 2 opened"
    assert parser_manager.ring_number("3", None) == 3, "Ring number 3 opened"

    assert parser_manager.ring_number("1", None) == 1, "Ring number 1 closed"
    assert parser_manager.ring_number("2", None) == 2, "Ring number 2 closed"
    assert parser_manager.ring_number("3", None) == 3, "Ring number 3 closed"

    assert parser_manager.ring_number("4", None) == 4, "Ring number 4 opened"
    assert parser_manager.ring_number("4", None) == 4, "Ring number 4 closed"

    assert parser_manager.ring_number("%", "5") == 5, "Ring number 5 opened"
    assert parser_manager.ring_number("%", "5") == 5, "Ring number 5 closed"

    with pytest.raises(ParserException) as exc_info:
        parser_manager.ring_number("A", None)
        assert (
            ParserException(
                rule="ring_number", parameter="A None", message="Invalid ring number"
            )
            == exc_info.value
        ), "Should raise ParserException for invalid ring number"

        parser_manager.ring_number("1", None)
        assert (
            ParserException(
                rule="ring_number",
                parameter="1 None",
                message="Ring number already closed",
            )
            == exc_info.value
        ), "Should raise ParserException for already closed ring number"

        parser_manager.ring_number("%", None)
        assert (
            ParserException(
                rule="ring_number", parameter="% None", message="Invalid ring number"
            )
            == exc_info.value
        ), "Should raise ParserException for invalid ring number symbol"

        parser_manager.ring_number("1", "2", "3")
        assert (
            ParserException(
                rule="ring_number",
                parameter="1 2 3",
                message="Ring number cannot have more than one digit after the first",
            )
            == exc_info.value
        ), "Should raise ParserException for ring number with more than one digit after the first"
