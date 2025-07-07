import os
import re
from dataclasses import dataclass, field
from json import load
from pathlib import Path
from typing import List, Optional

from .atomic import Atom, BracketAtom


@dataclass
class Graph:
    """
    Represents a molecular structure as a graph, where atoms are nodes and bonds are edges.

    Uses an adjacency list representation for the graph and tracks cycles.

    Attributes:
        adjacency_list (dict[Atom, List[tuple[Atom, str]]]): Adjacency list where keys are Atom objects
                                                            and values are lists of (neighbor_atom, bond_type) tuples.
        cycles (List[List[Atom]]): A list of cycles detected in the graph, where each cycle is a list of Atom objects.
    """

    adjacency_list: dict[Atom, List[tuple[Atom, str]]] = field(default_factory=dict)
    cycles: List[List[Atom]] = field(default_factory=list)

    def add_cycle(self, cycle: List[Atom]):
        """
        Adds a cycle (a list of connected atoms forming a ring) to the graph's cycle list.

        Args:
            cycle (List[Atom]): A list of Atom objects representing the cycle.
        """
        self.cycles.append(cycle)

    def add_edge(self, atom1: Atom, atom2: Atom, bond_type: str = "-"):
        """
        Adds an edge (bond) between two atoms in the graph.

        Args:
            atom1 (Atom): The first atom.
            atom2 (Atom): The second atom.
            bond_type (str): The type of bond (e.g., "-" for single, "=" for double, "#" for triple, ":" for aromatic).
        """
        if atom1 not in self.adjacency_list:
            self.adjacency_list[atom1] = []
        if atom2 not in self.adjacency_list:
            self.adjacency_list[atom2] = []
        self.adjacency_list[atom1].append((atom2, bond_type))
        self.adjacency_list[atom2].append((atom1, bond_type))

    def get_acyclic_subgraphs(self) -> List[List[Atom]]:
        """
        Identifies and returns all acyclic subgraphs within the molecular graph.

        This method uses a depth-first search (DFS) to traverse the graph and identify
        connected components that do not contain cycles.

        Returns:
            List[List[Atom]]: A list of lists, where each inner list represents an acyclic subgraph
                              and contains the Atom objects belonging to that subgraph.
        """
        visited = set()
        acyclic_subgraphs = []

        def dfs(atom: Atom, current_subgraph: List[Atom]):
            visited.add(atom)
            current_subgraph.append(atom)

            for neighbor, _ in self.adjacency_list.get(atom, []):
                if neighbor not in visited:
                    dfs(neighbor, current_subgraph)

        for atom in self.adjacency_list:
            if atom not in visited:
                current_subgraph = []
                dfs(atom, current_subgraph)
                acyclic_subgraphs.append(current_subgraph)

        return acyclic_subgraphs

    def check_valency_for_aba(self) -> bool:
        """
        Checks the valency of all acyclic BracketAtom instances within the graph.

        Iterates through all atoms in acyclic subgraphs and verifies their valency
        using the compute_valency method of BracketAtom.

        Returns:
            bool: True if all acyclic BracketAtom instances have valid valency, False otherwise.
        """
        acyclic_subgraphs = self.get_acyclic_subgraphs()
        bracket_atoms: List[BracketAtom] = [
            atom
            for subgraph in acyclic_subgraphs
            for atom in subgraph
            if isinstance(atom, BracketAtom)
        ]

        for atom in bracket_atoms:
            bonds = len(self.adjacency_list.get(atom, []))
            if not atom.compute_valency(bonds=bonds):
                return False
        return True

    def huckel(self) -> bool:
        """
        Applies Huckel's rule (4n + 2 pi electrons) to determine if cyclic systems are aromatic.

        This method calculates the total number of pi electrons within each detected cycle
        based on bond types and checks if the 4n+2 rule is satisfied.

        Returns:
            bool: True if all cycles satisfy Huckel's rule, False otherwise.
        """

        for atoms_in_cycle in self.cycles:
            pi_electrons = 0
            for i, atom in enumerate(atoms_in_cycle):
                next_atom = atoms_in_cycle[(i + 1) % len(atoms_in_cycle)]
                bond_type = None
                for neighbor, b_type in self.adjacency_list.get(atom, []):
                    if neighbor == next_atom:
                        bond_type = b_type
                        break

                if bond_type == ":":  # Aromatic bond
                    pi_electrons += 1
                elif bond_type == "=":  # Double bond
                    pi_electrons += 2
                elif bond_type == "#":  # Triple bond
                    pi_electrons += 4

            if not (pi_electrons - 2) % 4 == 0:
                return False

        return True


class Chemistry:
    """
    Class for handling Chemistry Domain Logic.

    Attributes:
        organic_atoms: A fixed list of all possible organic atoms
        pt_symbols: All the periodic table symbols
        look_up_table: A look up table for all atoms
        mol_graph: A graph representation of the molecule
    """

    organic_atoms: set[str]
    pt_symbols: List[str]
    look_up_table: dict[str, Atom]
    mol_graph: Graph

    def number_of_electrons_per_bond(self, bond: str) -> int:
        """
        Get the number of electrons per bond.

        Args:
            bond: The bond to get the number of electrons from.

        Returns:
            The number of electrons per bond.
        """
        bonds = {"=": 2, "#": 3, "$": 4, "/": 1, "\\": 1, "-": 1, ".": 0}

        if bond in bonds:
            return bonds[bond]

        raise Exception(f"Invalid Bond {bond}")

    def __init__(self, periodic_table_path=None):
        if periodic_table_path is None:
            periodic_table_path = os.path.join(
                os.path.dirname(__file__), "..", "periodic-table-lookup.json"
            )

        upper_organic_atoms = {"N", "O", "P", "H", "S", "F", "Cl", "Br", "I", "C", "B"}

        with open(periodic_table_path) as JSON:  # loads the periodic table json
            look_up_table_json = dict(load(JSON))

        # set of all organic atoms
        self.organic_atoms = {
            atom.lower() for atom in upper_organic_atoms
        } | upper_organic_atoms

        # set of all periodic table symbols
        self.pt_symbols = [
            look_up_table_json[x]["symbol"] for x in look_up_table_json["order"]
        ]

        look_up_table_json.pop("order")

        # creates a look up table for all atoms
        self.look_up_table = {
            x["symbol"]: Atom(
                symbol=x["symbol"], electron_configuration=x["electron_configuration"]
            )
            for x in look_up_table_json.values()
        }

        mol_graph = Graph()

    def Atom(self, symbol: str, aromatic: bool = False) -> Atom:
        """
        Create an Atom with the given symbol, charge and hidrogens.

        Args:
            symbol: The symbol of the atom
            chiral: The chiral of the atom
            hcount: The amount of hydrogens in the atom
            charge: The charge of the atom
            map: The map of the atom
        """
        if symbol.title() != symbol:
            symbol = symbol.title()
        if symbol not in self.pt_symbols:
            raise Exception(f"Invalid Symbol {symbol}")
        base_atom = self.look_up_table[symbol]
        return Atom(
            symbol,
            base_atom.electron_configuration,
            aromatic=aromatic,
        )

    def BracketAtom(self, symbol: str, **kwargs) -> BracketAtom:
        """
        Create a Bracket Atom with the given symbol, charge and hidrogens.

        Args:
            symbol: The symbol of the atom
            chiral: The chiral of the atom
            hcount: The amount of hydrogens in the atom
            charge: The charge of the atom
            map: The map of the atom
        """
        if symbol not in self.pt_symbols:
            raise Exception(f"Invalid Symbol {symbol}")
        base_atom = self.look_up_table[symbol]
        if symbol.title() != symbol:  # if the symbol is not lowercase
            symbol = symbol.title()
            kwargs["aromatic"] = True

        hidrogens = kwargs.pop("hidrogens", None)
        charge = kwargs.pop("charge", None)

        return BracketAtom(
            symbol,
            (
                base_atom.electron_configuration
                if base_atom.electron_configuration
                else ""
            )
            .replace("[He] ", "")
            .replace("[Ne] ", "")
            .replace("[Ar] ", "")
            .replace("[Kr] ", "")
            .replace("[Xe] ", "")
            .replace("[Rn] ", ""),
            hidrogens=hidrogens,
            charge=charge,
            **kwargs,
        )

    def validate_valency_bracket(
        self,
        isotope: Optional[int],
        symbol: str,
        chiral: Optional[int],
        hcount: Optional[int],
        charge: Optional[int],
        map: Optional[int],
    ) -> bool:
        """
        Validate the valency of a single bracket atom (if it is not a part of a aromatic ring)
        Args:
            isotope: The isotope of the atom
            symbol: The symbol of the atom
            chiral: The chiral of the atom
            hcount: The amount of hydrogens in the atom
            charge: The charge of the atom
            map: The map of the atom
        Returns:
            If the valency of the atom is valid
        """
        return self.BracketAtom(
            symbol=symbol,
            charge=charge if charge is not None else 0,
            hidrogens=hcount if hcount is not None else 0,
        ).compute_valency()

    def validate_aromacity(self) -> bool:
        """
        Validate the aromaticity of a molecule.
        This method checks if the molecule is aromatic by checking if it satisfies Huckel's rule.
        Returns:
            If the molecule is aromatic
        """
        return self.mol_graph.huckel() and self.mol_graph.check_valency_for_aba()


chemistry = Chemistry()
