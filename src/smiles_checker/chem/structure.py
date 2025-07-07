import os
from dataclasses import dataclass, field
from json import load
from pathlib import Path
from typing import List, Optional, Tuple

from .atomic import Atom, BracketAtom
from .chemistry import _parse_electron_configuration


@dataclass
class Graph:
    """
    Class for handling molecular structures as graphs.
    This class is used to represent molecules as graphs, where atoms are nodes
    and bonds are edges.

    It uses the adjacency list representation for the graph.
    It also needs to keep track of the cycles in the graph, which needs to be done
    using a auxiliary data structure.
    """

    adjacency_list: dict[Atom, List[Atom]] = field(default_factory=dict)
    cycles: List[List[Atom]] = field(default_factory=list)

    def add_cycle(self, cycle: List[Atom]):
        """
        Add a cycle to the graph.
        A cycle is a list of atoms that are connected in a circular manner.
        """
        self.cycles.append(cycle)

    def add_edge(self, atom1: Atom, atom2: Atom):
        """
        Add an edge between two atoms in the graph.
        """
        if atom1 not in self.adjacency_list:
            self.adjacency_list[atom1] = []
        if atom2 not in self.adjacency_list:
            self.adjacency_list[atom2] = []
        self.adjacency_list[atom1].append(atom2)
        self.adjacency_list[atom2].append(atom1)

    def get_acyclic_subgraphs(self) -> List[List[Atom]]:
        """
        Get all acyclic subgraphs in the graph.
        This method uses a depth-first search to find all acyclic subgraphs in the graph.
        """
        visited = set()
        acyclic_subgraphs = []

        def dfs(atom: Atom, current_subgraph: List[Atom]):
            visited.add(atom)
            current_subgraph.append(atom)

            for neighbor in self.adjacency_list.get(atom, []):
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
        Checks the valency for acyclic bracket atoms (ABA).
        """
        acyclic_subgraphs = self.get_acyclic_subgraphs()
        bracket_atoms: List[BracketAtom] = [
            atom
            for subgraph in acyclic_subgraphs
            for atom in subgraph
            if isinstance(atom, BracketAtom)
        ]

        for atom in bracket_atoms:
            if not atom.compute_valency():
                return False
        return True

    def huckel(self) -> bool:
        """
        Computes the Huckel's 4n + 2 rule for aromaticity.
        """

        for atoms in self.cycles:
            pi_electrons = sum(
                atom.get_total_electrons_in_subshell("p") for atom in atoms
            )
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
            print(symbol)
        if symbol not in self.pt_symbols:
            raise Exception(f"Invalid Symbol {symbol}")
        base_atom = self.look_up_table[symbol]
        return Atom(
            symbol=symbol,
            electron_configuration=base_atom.electron_configuration,
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
        return BracketAtom(
            symbol=symbol,
            electron_configuration=base_atom.electron_configuration,
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
