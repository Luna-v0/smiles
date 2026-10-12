"""Molecular graph structure for representing SMILES molecules."""

from typing import Dict, List, Tuple

from chem.atomic import Atom


class MolecularGraph:
    """
    Graph representation of a molecule.

    This is a plain adjacency-list container.  Valence and aromaticity are
    validated by the chemistry backend (``chem.adapter``) from this graph, so
    no octet/Hückel perception lives here.

    Attributes:
        adjacency_list: Dictionary mapping atoms to list of (neighbor, bond_type) tuples.
        cycles: List of cycles, where each cycle is a list of atoms in order.
    """

    def __init__(self):
        """Initialize an empty molecular graph."""
        self.adjacency_list: Dict[Atom, List[Tuple[Atom, str]]] = {}
        self.cycles: List[List[Atom]] = []

    def add_edge(self, atom1: Atom, atom2: Atom, bond: str = "-", bond_type: str = None):
        """
        Add an edge between two atoms.

        Args:
            atom1: First atom.
            atom2: Second atom.
            bond: Bond type (legacy parameter name).
            bond_type: Bond type (preferred parameter name).
        """
        # Support both parameter names
        if bond_type is None:
            bond_type = bond

        # Default to single bond if not specified
        if bond_type is None:
            bond_type = "-"

        # Initialize adjacency lists if needed
        if atom1 not in self.adjacency_list:
            self.adjacency_list[atom1] = []
        if atom2 not in self.adjacency_list:
            self.adjacency_list[atom2] = []

        # Add edge (undirected)
        self.adjacency_list[atom1].append((atom2, bond_type))
        self.adjacency_list[atom2].append((atom1, bond_type))

    def update_edge(self, atom1: Atom, atom2: Atom, bond_type: str):
        """
        Update the bond type of an existing edge in both directions.

        Args:
            atom1: First atom.
            atom2: Second atom.
            bond_type: New bond type to set.
        """
        for a, b in ((atom1, atom2), (atom2, atom1)):
            neighbors = self.adjacency_list.get(a, [])
            for i, (neighbor, _) in enumerate(neighbors):
                if neighbor is b:
                    neighbors[i] = (b, bond_type)
                    break

    def remove_edge(self, atom1: Atom, atom2: Atom):
        """
        Remove the edge between two atoms in both directions (no-op if absent).

        Args:
            atom1: First atom.
            atom2: Second atom.
        """
        for a, b in ((atom1, atom2), (atom2, atom1)):
            if a in self.adjacency_list:
                self.adjacency_list[a] = [
                    (n, bt) for (n, bt) in self.adjacency_list[a] if n is not b
                ]

    def add_cycle(self, cycle: List[Atom]):
        """
        Add a cycle to the graph.

        Args:
            cycle: List of atoms forming a cycle in order.
        """
        if cycle not in self.cycles:
            self.cycles.append(cycle)


# Alias for backward compatibility with tests
Graph = MolecularGraph
