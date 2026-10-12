"""Tests for the MolecularGraph container.

The graph is now a plain adjacency-list holder (valence/aromaticity live in the
chemistry backend), so only its edge/cycle bookkeeping is tested here.  The old
octet/Hückel graph tests were removed with that code.
"""

from chem.atomic import Atom
from chem.structure import Graph


def test_add_edge():
    graph = Graph()
    atom1 = Atom("C")
    atom2 = Atom("N")
    graph.add_edge(atom1, atom2, bond_type="-")
    assert atom2 in [a for a, b in graph.adjacency_list[atom1]]
    assert atom1 in [a for a, b in graph.adjacency_list[atom2]]


def test_add_cycle():
    graph = Graph()
    atom1 = Atom("C")
    atom2 = Atom("N")
    atom3 = Atom("O")
    cycle = [atom1, atom2, atom3]
    graph.add_cycle(cycle)
    assert cycle in graph.cycles
