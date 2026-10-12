"""Atomic classes for representing atoms in SMILES molecules.

These are plain data holders for the molecular graph.  Valence/aromaticity are
validated by the chemistry backend (see ``chem.adapter``), so no electron-count
or octet machinery lives here anymore.
"""

from dataclasses import dataclass, field
from typing import Optional

# Global counter for unique atom IDs
_atom_id_counter = 0


def _get_next_atom_id():
    """Get next unique atom ID."""
    global _atom_id_counter
    _atom_id_counter += 1
    return _atom_id_counter


@dataclass
class Atom:
    """
    Represents a basic atom in a SMILES string.

    Attributes:
        symbol: The atomic symbol (always uppercase).
        aromatic: Whether the atom is aromatic (lowercase symbol).
        atom_id: Unique per-instance id used for graph identity.
    """

    symbol: str
    aromatic: bool = False
    atom_id: int = field(default_factory=_get_next_atom_id, init=False)

    def __post_init__(self):
        """Normalize symbol and determine aromaticity."""
        original_symbol = self.symbol
        self.symbol = self.symbol.upper()
        # Aromatic if the original symbol was lowercase.
        self.aromatic = original_symbol.islower()

    def __eq__(self, other):
        """Check equality based on atom_id for graph purposes."""
        if not isinstance(other, Atom):
            return False
        return self.atom_id == other.atom_id

    def __hash__(self):
        """Make Atom hashable using atom_id."""
        return hash(self.atom_id)

    def __repr__(self):
        """String representation."""
        return f"Atom(symbol='{self.symbol}', aromatic={self.aromatic})"


@dataclass
class BracketAtom(Atom):
    """
    Represents a bracket atom with additional properties.

    Attributes:
        isotope: Isotope number if specified.
        chiral: Chiral rotation designator.
        hcount: Hydrogen count (``hidrogens`` alias accepted for compatibility).
        charge: Atomic charge.
        mol_map: Molecule mapping number.
    """

    isotope: Optional[int] = None
    chiral: Optional[str] = None
    hcount: Optional[int] = None
    charge: Optional[int] = None
    mol_map: Optional[int] = None

    def __post_init__(self):
        """Normalize symbol/aromaticity, honouring an explicit aromatic flag."""
        # Handle 'hidrogens' alias (used by some callers/tests).
        if hasattr(self, 'hidrogens'):
            self.hcount = self.hidrogens
            delattr(self, 'hidrogens')

        # Store explicit aromatic flag before super() recomputes it from case.
        explicit_aromatic = getattr(self, 'aromatic', None)
        super().__post_init__()
        if explicit_aromatic is not None:
            self.aromatic = explicit_aromatic

    def __eq__(self, other):
        """Check equality based on atom_id for graph purposes."""
        if not isinstance(other, BracketAtom):
            return False
        return self.atom_id == other.atom_id

    def __hash__(self):
        """Make BracketAtom hashable using atom_id."""
        return hash(self.atom_id)

    def __repr__(self):
        """String representation."""
        attrs = []
        if self.isotope is not None:
            attrs.append(f"isotope={self.isotope}")
        attrs.append(f"symbol='{self.symbol}'")
        if self.chiral:
            attrs.append(f"chiral='{self.chiral}'")
        if self.hcount is not None:
            attrs.append(f"hcount={self.hcount}")
        if self.charge is not None:
            attrs.append(f"charge={self.charge}")
        if self.mol_map is not None:
            attrs.append(f"mol_map={self.mol_map}")
        return f"BracketAtom({', '.join(attrs)})"
