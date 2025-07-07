import os
from dataclasses import dataclass, field
from json import load
from pathlib import Path
from typing import List, Optional
import re

from .atomic import Atom, BracketAtom
from .structure import Graph


class Chemistry:
    """
    Handles chemical domain logic, including atom and molecule properties and validation.

    Attributes:
        organic_atoms (set[str]): A set of common organic atom symbols (case-sensitive and lowercase).
        pt_symbols (list[str]): A list of all valid periodic table symbols.
        look_up_table (dict[str, Atom]): A lookup table for Atom objects by their symbol.
        mol_graph (Graph): A graph representation of the molecule being processed.
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

    def __init__(self, periodic_table_path: Optional[str] = None):
        """
        Initializes the Chemistry handler, loading periodic table data.

        Args:
            periodic_table_path (Optional[str]): Path to the periodic table JSON file.
                                                 Defaults to '../periodic-table-lookup.json'.
        """
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
                symbol=x["symbol"], electron_configuration=x.get("electron_configuration", "").replace("[He] ", "").replace("[Ne] ", "").replace("[Ar] ", "").replace("[Kr] ", "").replace("[Xe] ", "").replace("[Rn] ", "")
            )
            for x in look_up_table_json.values()
        }

        self.mol_graph = Graph()

    def Atom(self, symbol: str, aromatic: bool = False) -> Atom:
        """
        Creates an Atom object from a given symbol.

        Args:
            symbol (str): The chemical symbol of the atom (e.g., "C", "N", "O").
            aromatic (bool): Indicates if the atom is considered aromatic.

        Returns:
            Atom: An Atom object initialized with properties from the periodic table.

        Raises:
            Exception: If the provided symbol is not a valid periodic table element.
        """
        processed_symbol = symbol.title()
        if processed_symbol not in self.pt_symbols:
            raise Exception(f"Invalid Symbol {symbol}")
        base_atom = self.look_up_table[processed_symbol]
        return Atom(
            processed_symbol,
            (base_atom.electron_configuration if base_atom.electron_configuration else "").replace("[He] ", "").replace("[Ne] ", "").replace("[Ar] ", "").replace("[Kr] ", "").replace("[Xe] ", "").replace("[Rn] ", ""),
            aromatic=aromatic,
        )

    def BracketAtom(self, symbol: str, **kwargs) -> BracketAtom:
        """
        Creates a BracketAtom object with additional properties.

        Args:
            symbol (str): The chemical symbol of the atom.
            **kwargs: Additional properties like hidrogens, charge, isotope, chiral, mol_map.

        Returns:
            BracketAtom: A BracketAtom object initialized with specified properties.

        Raises:
            Exception: If the provided symbol is not a valid periodic table element.
        """
        processed_symbol = symbol.title()
        if processed_symbol not in self.pt_symbols:
            raise Exception(f"Invalid Symbol {symbol}")
        base_atom = self.look_up_table[processed_symbol]
        if symbol.title() != symbol:  # if the symbol is not lowercase
            kwargs["aromatic"] = True

        hidrogens = kwargs.pop("hidrogens", None)
        charge = kwargs.pop("charge", None)

        return BracketAtom(
            processed_symbol,
            (base_atom.electron_configuration if base_atom.electron_configuration else "").replace("[He] ", "").replace("[Ne] ", "").replace("[Ar] ", "").replace("[Kr] ", "").replace("[Xe] ", "").replace("[Rn] ", ""),
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
        Validates the valency of a single bracket atom.

        Args:
            isotope (Optional[int]): The isotope of the atom.
            symbol (str): The symbol of the atom.
            chiral (Optional[int]): The chiral specification of the atom.
            hcount (Optional[int]): The amount of hydrogens in the atom.
            charge (Optional[int]): The charge of the atom.
            map (Optional[int]): The map of the atom.

        Returns:
            bool: True if the valency of the atom is valid, False otherwise.
        """
        return self.BracketAtom(
            symbol=symbol,
            charge=charge if charge is not None else 0,
            hidrogens=hcount if hcount is not None else 0,
        ).compute_valency()

    def validate_aromacity(self) -> bool:
        """
        Validates the aromaticity of the molecule represented by the internal graph.

        This method checks if the molecule satisfies Huckel's rule (4n + 2 pi electrons in a cyclic, planar, fully conjugated system).

        Returns:
            bool: True if the molecule is aromatic, False otherwise.
        """
        return self.mol_graph.huckel() and self.mol_graph.check_valency_for_aba()


chemistry = Chemistry()
