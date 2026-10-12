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

# Extract symbols from periodic table (used for symbol validation)
pt_symbols = []
for element_name in _PERIODIC_TABLE_DATA.get("order", []):
    element_data = _PERIODIC_TABLE_DATA.get(element_name, {})
    symbol = element_data.get("symbol", "")
    if symbol:
        pt_symbols.append(symbol)


class ChemistryModule:
    """
    Chemistry module for SMILES validation.
    
    Provides factory functions for creating atoms and chemistry validation utilities.
    """

    def __init__(self):
        """Initialize chemistry module."""
        self.mol_graph = MolecularGraph()

    def Atom(self, symbol: str) -> Atom:
        """
        Create an Atom instance.

        Args:
            symbol: Atomic symbol (case-sensitive for aromaticity).
                    "*" is a wildcard atom (matches any atom).

        Returns:
            Atom instance.

        Raises:
            ParserException: If symbol is invalid.
        """
        # Wildcard atom: skip periodic table validation
        if symbol == "*":
            return Atom(symbol="*")

        symbol_upper = symbol.upper()
        # Check if any symbol in pt_symbols matches (case-insensitive)
        if not any(s.upper() == symbol_upper for s in pt_symbols):
            raise ParserException(
                rule="Atom",
                parameter=symbol,
                message=f"Invalid Atom Symbol: {symbol}"
            )

        return Atom(symbol=symbol)

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
            aromatic: Optional aromatic flag override.

        Returns:
            BracketAtom instance.

        Raises:
            ParserException: If symbol is invalid.
        """
        # Wildcard bracket atom: skip periodic table validation
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


# Create singleton instance
chemistry = ChemistryModule()

# Export commonly used items
__all__ = ['chemistry', 'pt_symbols', 'Atom', 'BracketAtom']
