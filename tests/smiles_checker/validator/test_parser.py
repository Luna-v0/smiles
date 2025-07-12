import pytest
from smiles_checker.validator.parser_manager import ParserManager, ParserException
from smiles_checker.chem.atomic import Atom, BracketAtom
from smiles_checker.chem.chemistry import chemistry as chem


@pytest.fixture
def parser_manager():
    return ParserManager()

def test_validate_unclosed_rings(parser_manager: ParserManager):
    parser_manager.current_open_rnum = [1]
    with pytest.raises(ParserException) as exc_info:
        parser_manager.validate()
    assert exc_info.value.message == "Unclosed ring numbers"

def test_validate_empty_chain(parser_manager: ParserManager):
    assert parser_manager.validate() is True

def test_internal_bracket(parser_manager: ParserManager):
    """ """
    assert parser_manager.internal_bracket(symbol="C") == chem.BracketAtom(
        symbol="C"
    ), "Internal bracket should return BracketAtom('C')"

def test_atom(parser_manager: ParserManager):
    """
    Test the atom parser with valid and invalid inputs.
    """

    assert (
        parser_manager.atom("C").symbol == Atom("C").symbol
    ), "Parsing 'C' should return Atom('C')"
    assert (
        parser_manager.atom("Na").symbol == Atom("Na").symbol
    ), "Parsing 'Na' should return Atom('Na')"

    with pytest.raises(ParserException) as exc_info:
        parser_manager.atom("Xx")
    assert exc_info.value.message == "Invalid Symbol Xx"

def test_chiral(parser_manager: ParserManager):
    """
    Test the chiral parser.
    """
    assert (
        parser_manager.chiral(chiral1="@", chiral2=None) == "clockwise"
    ), "Chiral without second symbol should be clockwise"
    assert (
        parser_manager.chiral(chiral1="@", chiral2="@") == "counter-clockwise"
    ), "Chiral with second symbol should be counterclockwise"

def test_hcount(parser_manager: ParserManager):
    """
    Test the hydrogen count parser.
    """
    assert (
        parser_manager.hcount(_="H", digit="1") == 1
    ), "Hydrogen count '1' should return 1"
    assert (
        parser_manager.hcount(_="H", digit="2") == 2
    ), "Hydrogen count '2' should return 2"
    assert (
        parser_manager.hcount(_="H", digit="3") == 3
    ), "Hydrogen count '3' should return 3"

    with pytest.raises(ParserException) as exc_info:
        parser_manager.hcount(_="H", digit="A")
    assert exc_info.value.message == "Invalid hydrogen count"

def test_ring_number(parser_manager: ParserManager):
    """
    Test the ring number parser.
    """
    parser_manager.atom("C")  # Set a last_atom for ring opening
    assert (
        parser_manager.ring_number(
            ring_number_or_symbol="1", ring_number1=None, ring_number2=None
        )
        == 1
    ), "Ring number 1 opened"
    parser_manager.atom("C")  # Set a last_atom for ring opening
    assert (
        parser_manager.ring_number(
            ring_number_or_symbol="2", ring_number1=None, ring_number2=None
        )
        == 2
    ), "Ring number 2 opened"
    parser_manager.atom("C")  # Set a last_atom for ring opening
    assert (
        parser_manager.ring_number(
            ring_number_or_symbol="3", ring_number1=None, ring_number2=None
        )
        == 3
    ), "Ring number 3 opened"

    parser_manager.atom("C")  # Set a last_atom for ring opening
    assert (
        parser_manager.ring_number(
            ring_number_or_symbol="1", ring_number1=None, ring_number2=None
        )
        == 1
    ), "Ring number 1 closed"
    parser_manager.atom("C")  # Set a last_atom for ring opening
    assert (
        parser_manager.ring_number(
            ring_number_or_symbol="2", ring_number1=None, ring_number2=None
        )
        == 2
    ), "Ring number 2 closed"
    parser_manager.atom("C")  # Set a last_atom for ring opening
    assert (
        parser_manager.ring_number(
            ring_number_or_symbol="3", ring_number1=None, ring_number2=None
        )
        == 3
    ), "Ring number 3 closed"

    parser_manager.atom("C")  # Set a last_atom for ring opening
    assert (
        parser_manager.ring_number(
            ring_number_or_symbol="4", ring_number1=None, ring_number2=None
        )
        == 4
    ), "Ring number 4 opened"
    parser_manager.atom("C")  # Set a last_atom for ring opening
    assert (
        parser_manager.ring_number(
            ring_number_or_symbol="4", ring_number1=None, ring_number2=None
        )
        == 4
    ), "Ring number 4 closed"

    with pytest.raises(ParserException) as exc_info:
        parser_manager.ring_number(
            ring_number_or_symbol="%", ring_number1="5", ring_number2=None
        )
    assert exc_info.value.message == "'%' must be followed by two digits for ring numbers > 9."

    with pytest.raises(ParserException) as exc_info:
        parser_manager.ring_number(
            ring_number_or_symbol="A", ring_number1=None, ring_number2=None
        )
    assert exc_info.value.message == "Invalid character for ring number."