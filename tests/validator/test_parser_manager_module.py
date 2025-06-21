import pytest
from validator.parser_manager import ParserManager, ParserException
from chem import chem


def test_listify():
    pm = ParserManager()
    assert pm.listify(1, None) == 1
    assert pm.listify(1, [2, 3]) == [1, 2, 3]
    assert pm.listify(1, 4) == [1, 4]


def test_hcount_and_int():
    pm = ParserManager()
    assert pm.hcount(None, None) == 1
    assert pm.hcount(None, "3") == 3
    assert pm.int(["1", "2"]) == 12


def test_charge_and_chiral_and_fifteen():
    pm = ParserManager()
    assert pm.charge("+", None) == 1
    assert pm.charge("-", None) == -1
    assert pm.charge("+", "-") == -2
    assert pm.charge("+", "+") == 2
    assert pm.charge("-", 3) == -3
    assert pm.chiral("@", None) == "clockwise"
    assert pm.chiral("@", "@")=="counterclockwise"
    assert pm.fifteen("1", None) == 1
    assert pm.fifteen("1", "2") == 12
    with pytest.raises(ParserException):
        pm.fifteen("1", "7")


def test_ring_number_cycle():
    pm = ParserManager()
    assert pm.ring_number("1", None, None) == 1
    assert pm.current_open_rnum == [1]
    # closing the ring
    assert pm.ring_number("1", None, None) == 1
    assert pm.current_closed_rnum == [1]
    with pytest.raises(ParserException):
        pm.ring_number("1", None, None)


def test_atom_and_chain_and_internal_bracket():
    pm = ParserManager()
    atom = pm.atom("C")
    assert atom.aromatic is True
    br = pm.internal_bracket(None, "C")
    assert pm.current_chain[-1] == br
    assert isinstance(br, type(chem.BracketAtom("C")))
    # chain without bond just returns passed element
    assert pm.chain(None, atom, None, None) == atom
    assert pm.chain(None, None, 2, None) == 2
    assert pm.chain(None, None, None, "dot") == "dot"
    with pytest.raises(Exception):
        pm.chain(":", "C", None, None)
    result = pm.chain("-", atom, None, None)
    assert result[1] == 1
