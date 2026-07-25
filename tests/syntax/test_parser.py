import pytest
from syntax.parser_manager import ParserManager, ParserException
from chem.atomic import Atom, BracketAtom
from chem.chemistry import chemistry as chem


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
    result = parser_manager.internal_bracket(symbol="C")
    expected = chem.BracketAtom(symbol="C")
    # BracketAtoms are compared by atom_id, so instances won't be equal
    # But they should have the same symbol
    assert result.symbol == expected.symbol, "Internal bracket should return BracketAtom with symbol 'C'"

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


# ============================================================================
# Tests for branched molecule parsing (CC(C)C bug fix)
# ============================================================================

class TestBranchedMolecules:
    """Tests for branched molecule graph connectivity.
    
    These tests verify that:
    1. Atoms inside branches connect to the branch start atom
    2. Atoms after branches connect to the branch start atom (restored by end_branch)
    
    Note: atom() doesn't create bonds between consecutive atoms - that's done
    by chains(). These tests focus on branch-specific bonding behavior.
    """
    
    @pytest.fixture
    def pm(self):
        """Fresh ParserManager for each test."""
        return ParserManager()
    
    def _count_edges(self, pm: ParserManager) -> int:
        """Count total edges in the graph (each edge is counted once)."""
        graph = pm.graph_builder.get_graph()
        return sum(len(neighbors) for neighbors in graph.adjacency_list.values()) // 2
    
    def _get_neighbor_ids(self, pm: ParserManager, atom) -> list:
        """Get list of neighbor atom IDs for a given atom."""
        graph = pm.graph_builder.get_graph()
        return [n.atom_id for n, _ in graph.adjacency_list.get(atom, [])]
    
    def test_branch_atom_connects_to_branch_start(self, pm: ParserManager):
        """Test that first atom inside a branch connects to branch_start_atom."""
        c1 = pm.atom("C")
        pm.start_branch()
        c2 = pm.atom("C")  # Should auto-connect to c1
        pm.end_branch(None)
        
        # c2 should be connected to c1
        assert self._count_edges(pm) == 1, "Branch atom should connect to branch start"
        assert c1.atom_id in self._get_neighbor_ids(pm, c2)
    
    def test_atom_after_branch_connects_to_branch_start(self, pm: ParserManager):
        """Test that atom after branch connects to branch_start_atom (via _connect_next_atom)."""
        c1 = pm.atom("C")
        pm.start_branch()
        c2 = pm.atom("C")  # Inside branch
        pm.end_branch(None)
        c3 = pm.atom("C")  # After branch - should connect to c1
        
        assert self._count_edges(pm) == 2
        c1_neighbors = self._get_neighbor_ids(pm, c1)
        assert c2.atom_id in c1_neighbors, "Branch atom should connect to c1"
        assert c3.atom_id in c1_neighbors, "Post-branch atom should connect to c1"
    
    def test_multiple_branches(self, pm: ParserManager):
        """Test C(C)(C)C pattern - multiple branches from same atom."""
        c1 = pm.atom("C")
        pm.start_branch()
        c2 = pm.atom("C")  # First branch
        pm.end_branch(None)
        pm.start_branch()
        c3 = pm.atom("C")  # Second branch
        pm.end_branch(None)
        c4 = pm.atom("C")  # After both branches
        
        assert self._count_edges(pm) == 3, "C(C)(C)C should have 3 edges"
        
        # C1 should be connected to C2, C3, and C4
        c1_neighbors = self._get_neighbor_ids(pm, c1)
        assert c2.atom_id in c1_neighbors
        assert c3.atom_id in c1_neighbors
        assert c4.atom_id in c1_neighbors
    
    def test_branch_with_explicit_bond(self, pm: ParserManager):
        """Test branch with explicit bond type (=O pattern)."""
        c1 = pm.atom("C")
        pm.start_branch()
        pm.save_branch_bond("=")
        o1 = pm.atom("O")  # Double-bonded O
        pm.end_branch(None)
        o2 = pm.atom("O")  # After branch
        
        assert self._count_edges(pm) == 2
        
        # Verify bond types
        graph = pm.graph_builder.get_graph()
        c1_neighbors = graph.adjacency_list.get(c1, [])
        for neighbor, bond_type in c1_neighbors:
            if neighbor == o1:
                assert bond_type == "=", "O inside branch should have double bond"
            elif neighbor == o2:
                assert bond_type == "-", "O after branch should have single bond"
    
    def test_nested_branch(self, pm: ParserManager):
        """Test nested branches C(C(C)C)C."""
        c1 = pm.atom("C")
        pm.start_branch()
        c2 = pm.atom("C")  # In first branch
        pm.start_branch()
        c3 = pm.atom("C")  # Nested branch
        pm.end_branch(None)
        c4 = pm.atom("C")  # After nested, still in first branch
        pm.end_branch(None)
        c5 = pm.atom("C")  # After all branches
        
        assert self._count_edges(pm) == 4
        
        # Verify connectivity:
        # c1 -> c2 (first branch)
        # c2 -> c3 (nested branch)
        # c2 -> c4 (after nested)
        # c1 -> c5 (after all branches)
        assert c2.atom_id in self._get_neighbor_ids(pm, c1)
        assert c3.atom_id in self._get_neighbor_ids(pm, c2)
        assert c4.atom_id in self._get_neighbor_ids(pm, c2)
        assert c5.atom_id in self._get_neighbor_ids(pm, c1)