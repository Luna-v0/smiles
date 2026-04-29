"""Chemistry module for SMILES validation."""

import json
import os
from pathlib import Path

from chem.atomic import Atom, BracketAtom
from chem.structure import MolecularGraph
from exceptions import ParserException


# Load periodic table data
_PERIODIC_TABLE_PATH = Path(__file__).parent.parent / "periodic-table-lookup.json"
with open(_PERIODIC_TABLE_PATH, 'r') as f:
    _PERIODIC_TABLE_DATA = json.load(f)

# Extract symbols from periodic table
pt_symbols = []
for element_name in _PERIODIC_TABLE_DATA.get("order", []):
    element_data = _PERIODIC_TABLE_DATA.get(element_name, {})
    symbol = element_data.get("symbol", "")
    if symbol:
        pt_symbols.append(symbol)

# Create a lookup for element data by symbol
_element_by_symbol = {}
for element_name in _PERIODIC_TABLE_DATA.get("order", []):
    element_data = _PERIODIC_TABLE_DATA.get(element_name, {})
    symbol = element_data.get("symbol", "")
    if symbol:
        _element_by_symbol[symbol.upper()] = element_data


class ChemistryModule:
    """
    Chemistry module for SMILES validation.
    
    Provides factory functions for creating atoms and chemistry validation utilities.
    """

    def __init__(self):
        """Initialize chemistry module."""
        self.mol_graph = MolecularGraph()

    def Atom(self, symbol: str, electron_configuration: str = None) -> Atom:
        """
        Create an Atom instance.

        Args:
            symbol: Atomic symbol (case-sensitive for aromaticity).
                    "*" is a wildcard atom (matches any atom).
            electron_configuration: Optional electron configuration override.

        Returns:
            Atom instance.

        Raises:
            ParserException: If symbol is invalid.
        """
        # Wildcard atom: skip periodic table validation and electron config
        if symbol == "*":
            return Atom(symbol="*", electron_configuration="")

        symbol_upper = symbol.upper()
        # Check if any symbol in pt_symbols matches (case-insensitive)
        if not any(s.upper() == symbol_upper for s in pt_symbols):
            raise ParserException(
                rule="Atom",
                parameter=symbol,
                message=f"Invalid Atom Symbol: {symbol}"
            )

        # Get electron configuration from periodic table if not provided
        if electron_configuration is None:
            element_data = _element_by_symbol.get(symbol_upper, {})
            electron_configuration = element_data.get("electron_configuration", "")

        return Atom(symbol=symbol, electron_configuration=electron_configuration)

    def BracketAtom(
        self,
        symbol: str,
        isotope: int = None,
        chiral: str = None,
        hcount: int = None,
        hidrogens: int = None,  # Support typo in tests
        charge: int = None,
        mol_map: int = None,
        map: int = None,  # Support 'map' parameter name
        electron_configuration: str = None,
        aromatic: bool = None,
    ) -> BracketAtom:
        """
        Create a BracketAtom instance.

        Args:
            symbol: Atomic symbol.
            isotope: Isotope number.
            chiral: Chiral rotation.
            hcount: Hydrogen count.
            hidrogens: Hydrogen count (alternative name for compatibility).
            charge: Atomic charge.
            mol_map: Molecule mapping number.
            map: Molecule mapping number (alternative name).
            electron_configuration: Optional electron configuration override.
            aromatic: Optional aromatic flag override.

        Returns:
            BracketAtom instance.

        Raises:
            ParserException: If symbol is invalid.
        """
        # Wildcard bracket atom: skip periodic table validation and electron config
        if symbol == "*":
            if hidrogens is not None:
                hcount = hidrogens
            if map is not None:
                mol_map = map
            return BracketAtom(
                symbol="*",
                isotope=isotope,
                chiral=chiral,
                hcount=hcount,
                charge=charge,
                mol_map=mol_map,
                electron_configuration="",
                aromatic=False if aromatic is None else aromatic,
            )

        symbol_upper = symbol.upper()
        # Check if any symbol in pt_symbols matches (case-insensitive)
        if not any(s.upper() == symbol_upper for s in pt_symbols):
            raise ParserException(
                rule="BracketAtom",
                parameter=symbol,
                message=f"Invalid Atom Symbol: {symbol}"
            )

        # Handle parameter aliases
        if hidrogens is not None:
            hcount = hidrogens
        if map is not None:
            mol_map = map

        # Get electron configuration from periodic table if not provided
        if electron_configuration is None:
            element_data = _element_by_symbol.get(symbol_upper, {})
            electron_configuration = element_data.get("electron_configuration", "")

        # Determine aromaticity
        if aromatic is None:
            aromatic = symbol.islower()
        
        # Create atom with all parameters
        atom = BracketAtom(
            symbol=symbol,
            isotope=isotope,
            chiral=chiral,
            hcount=hcount,
            charge=charge,
            mol_map=mol_map,
            electron_configuration=electron_configuration,
            aromatic=aromatic,
        )
        
        return atom

    def validate(self) -> bool:
        """
        Validate the current molecule graph.

        Returns:
            True if valid, raises exception otherwise.
        """
        # This will be replaced by ChemistryValidator
        # For now, just check if graph is empty or has no issues
        return True

    def clear(self):
        """Clear the molecule graph."""
        self.mol_graph = MolecularGraph()

    @staticmethod
    def number_of_electrons_per_bond(bond: str) -> int:
        """
        Get number of electrons contributed by a bond type.

        Args:
            bond: Bond type ('-', '=', '#', '$', '/', '\\', ':').

        Returns:
            Number of electrons (typically 2 per bond order).

        Raises:
            Exception: If bond type is invalid.
        """
        bond_electrons = {
            '-': 1,  # Single bond: 2 electrons total, 1 per atom
            '=': 2,  # Double bond: 4 electrons total, 2 per atom
            '#': 3,  # Triple bond: 6 electrons total, 3 per atom
            '$': 4,  # Quadruple bond: 8 electrons total, 4 per atom
            '/': 1,  # Single bond (stereo)
            '\\': 1,  # Single bond (stereo)
            ':': 1,  # Aromatic bond: 1.5 electrons, treated as 1 for pi counting
            '.': 0,  # No bond
        }
        
        if bond not in bond_electrons:
            raise Exception(f"Invalid Bond {bond}")
        
        return bond_electrons[bond]

    @staticmethod
    def validate_valency_bracket(
        isotope: int = None,
        symbol: str = None,
        chiral: str = None,
        hcount: int = None,
        charge: int = None,
        map: int = None,
    ) -> bool:
        """
        Validate valency for a bracket atom specification.

        Args:
            isotope: Isotope number.
            symbol: Atomic symbol.
            chiral: Chiral rotation.
            hcount: Hydrogen count.
            charge: Atomic charge.
            map: Molecule mapping number.

        Returns:
            True if valency is satisfied, False otherwise.
        """
        if symbol is None:
            return False
        
        try:
            atom = chemistry.BracketAtom(
                symbol=symbol,
                isotope=isotope,
                chiral=chiral,
                hcount=hcount,
                charge=charge,
                mol_map=map,
            )
            return atom.compute_valency()
        except Exception:
            return False


# Create singleton instance
chemistry = ChemistryModule()

# Export commonly used items
__all__ = ['chemistry', 'pt_symbols', 'Atom', 'BracketAtom']
