"""Atomic classes for representing atoms in SMILES molecules."""

import re
from dataclasses import dataclass, field
from typing import Optional, Tuple

from exceptions import ParserException

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
        electron_configuration: The electron configuration string from periodic table.
        valency_layer: The outermost shell number containing valence electrons.
        electrons_in_valency: Number of electrons in the valence shell.
        layers: Tuple of shell numbers that contain electrons.
        electrons_by_layers: Tuple of electron counts per shell.
    """

    symbol: str
    electron_configuration: str = ""
    aromatic: bool = False
    atom_id: int = field(default_factory=_get_next_atom_id, init=False)

    def __post_init__(self):
        """Normalize symbol and determine aromaticity."""
        # Store original to check if lowercase
        original_symbol = self.symbol
        self.symbol = self.symbol.upper()
        
        # Aromatic if original was lowercase
        self.aromatic = original_symbol.islower()
        
        # Parse electron configuration
        self._parse_electron_configuration()

    def _parse_electron_configuration(self):
        """Parse electron configuration to extract layer and electron information."""
        if not self.electron_configuration:
            self.valency_layer = 0
            self.electrons_in_valency = 0
            self.layers = ()
            self.electrons_by_layers = ()
            return

        # Pattern: "1s2 2s2 2p4" etc.
        pattern = r'(\d+)([spdf])(\d+)'
        matches = re.findall(pattern, self.electron_configuration)
        
        if not matches:
            self.valency_layer = 0
            self.electrons_in_valency = 0
            self.layers = ()
            self.electrons_by_layers = ()
            return

        # Group electrons by shell number
        shell_electrons = {}
        for shell_str, subshell, count_str in matches:
            shell_num = int(shell_str)
            count = int(count_str)
            shell_electrons[shell_num] = shell_electrons.get(shell_num, 0) + count

        if not shell_electrons:
            self.valency_layer = 0
            self.electrons_in_valency = 0
            self.layers = ()
            self.electrons_by_layers = ()
            return

        # Find valency layer (outermost shell)
        self.valency_layer = max(shell_electrons.keys())
        self.electrons_in_valency = shell_electrons[self.valency_layer]
        
        # Create tuples sorted by shell number
        sorted_shells = sorted(shell_electrons.keys())
        self.layers = tuple(sorted_shells)
        self.electrons_by_layers = tuple(shell_electrons[s] for s in sorted_shells)

    def get_electrons_in_specific_subshell(self, shell: int, subshell: str) -> int:
        """
        Get number of electrons in a specific subshell.

        Args:
            shell: Shell number (e.g., 1, 2, 3).
            subshell: Subshell type ('s', 'p', 'd', 'f').

        Returns:
            Number of electrons in the specified subshell.
        """
        if not self.electron_configuration:
            return 0

        pattern = rf'{shell}{re.escape(subshell)}(\d+)'
        match = re.search(pattern, self.electron_configuration)
        if match:
            return int(match.group(1))
        return 0

    def get_total_electrons_in_subshell(self, subshell: str) -> int:
        """
        Get total electrons across all shells for a subshell type.

        Args:
            subshell: Subshell type ('s', 'p', 'd', 'f').

        Returns:
            Total number of electrons in all shells of this subshell type.
        """
        if not self.electron_configuration:
            return 0

        pattern = rf'\d+{re.escape(subshell)}(\d+)'
        matches = re.findall(pattern, self.electron_configuration)
        return sum(int(count) for count in matches)

    def _next_subshell(self, subshell: str) -> str:
        """
        Get the next subshell in order (s -> p -> d -> f -> s).

        Args:
            subshell: Current subshell type.

        Returns:
            Next subshell type.
        """
        subshell_order = ['s', 'p', 'd', 'f']
        try:
            idx = subshell_order.index(subshell)
            if idx + 1 < len(subshell_order):
                return subshell_order[idx + 1]
            return 's'  # Wrap around
        except (ValueError, IndexError):
            return 's'

    def _max_electrons_in_subshell(self, subshell: str) -> int:
        """
        Get maximum electrons that can fit in a subshell.

        Args:
            subshell: Subshell type.

        Returns:
            Maximum electrons (s=2, p=6, d=10, f=14).
        """
        max_electrons = {'s': 2, 'p': 6, 'd': 10, 'f': 14}
        return max_electrons.get(subshell, 0)

    def __eq__(self, other):
        """Check equality based on atom_id for graph purposes."""
        if not isinstance(other, Atom):
            return False
        # Use atom_id for graph identity (each atom instance is unique)
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
        chiral: Chiral rotation ('clockwise' or 'counterclockwise').
        hcount: Hydrogen count (hidrogens in tests).
        charge: Atomic charge.
        mol_map: Molecule mapping number.
    """

    isotope: Optional[int] = None
    chiral: Optional[str] = None
    hcount: Optional[int] = None
    charge: Optional[int] = None
    mol_map: Optional[int] = None

    def __post_init__(self):
        """Initialize atom and adjust electron configuration for charge/hydrogens."""
        # Handle hidrogens parameter (typo in tests)
        if hasattr(self, 'hidrogens'):
            self.hcount = self.hidrogens
            delattr(self, 'hidrogens')
        
        # Store original aromatic flag before calling super
        explicit_aromatic = getattr(self, 'aromatic', None)
        
        super().__post_init__()
        
        # Restore explicit aromatic flag if it was set
        if explicit_aromatic is not None:
            self.aromatic = explicit_aromatic
        
        # Adjust electron configuration based on charge and hydrogens
        if self.charge is not None or self.hcount is not None:
            self._adjust_electrons_for_charge_and_hydrogens()

    def _adjust_electrons_for_charge_and_hydrogens(self):
        """Adjust electron configuration based on charge and hydrogen count."""
        if not self.electron_configuration:
            return

        # Calculate total electrons to add/remove
        # Negative charge = add electrons, positive = remove
        # Hydrogens contribute 1 electron each (bonding electron)
        electron_change = 0
        if self.charge is not None:
            electron_change -= self.charge  # Negative charge adds electrons
        if self.hcount is not None:
            electron_change += self.hcount  # Each H adds one electron

        if electron_change == 0:
            return

        # Parse current configuration
        pattern = r'(\d+)([spdf])(\d+)'
        matches = re.findall(pattern, self.electron_configuration)
        
        if not matches:
            return

        # Convert to list of (shell, subshell, count) tuples
        config_list = [(int(s), sub, int(c)) for s, sub, c in matches]
        
        # Sort by shell then subshell order
        subshell_order = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
        config_list.sort(key=lambda x: (x[0], subshell_order.get(x[1], 4)))

        # Apply electron changes
        if electron_change > 0:
            # Add electrons
            self._add_electrons(config_list, electron_change)
        else:
            # Remove electrons
            self._remove_electrons(config_list, abs(electron_change))

        # Rebuild electron configuration string
        self.electron_configuration = ' '.join(
            f'{shell}{subshell}{count}' 
            for shell, subshell, count in config_list
            if count > 0
        )
        
        # Re-parse to update valency info
        self._parse_electron_configuration()

    def _add_electrons(self, config_list: list, num_electrons: int):
        """Add electrons to configuration, creating new shells if needed."""
        # Electron filling order: 1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f, 5d, 6p, 7s, 5f, 6d, 7p
        # After p subshell, go to next shell's s (not same shell's d)
        
        for _ in range(num_electrons):
            added = False
            # Try to add to existing subshells first
            for i, (shell, subshell, count) in enumerate(config_list):
                max_electrons = self._max_electrons_in_subshell(subshell)
                if count < max_electrons:
                    config_list[i] = (shell, subshell, count + 1)
                    added = True
                    break
            
            if not added:
                # Need to add to a new subshell
                # Find the last shell and subshell
                if config_list:
                    last_shell, last_subshell, last_count = config_list[-1]
                    max_for_last = self._max_electrons_in_subshell(last_subshell)
                    
                    # If last subshell is full, move to next
                    if last_count >= max_for_last:
                        if last_subshell == 'p':
                            # After p, go to next shell's s
                            new_shell = last_shell + 1
                            config_list.append((new_shell, 's', 1))
                        else:
                            next_subshell = self._next_subshell(last_subshell)
                            if next_subshell == 's':
                                # New shell
                                new_shell = last_shell + 1
                                config_list.append((new_shell, 's', 1))
                            else:
                                # Same shell, next subshell
                                config_list.append((last_shell, next_subshell, 1))
                    else:
                        # Last subshell not full, add to it
                        config_list[-1] = (last_shell, last_subshell, last_count + 1)
                else:
                    # Empty config, start with 1s
                    config_list.append((1, 's', 1))

    def _remove_electrons(self, config_list: list, num_electrons: int):
        """Remove electrons from configuration, starting from outermost."""
        # Sort in reverse order (outermost first)
        subshell_order = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
        config_list.sort(key=lambda x: (-x[0], -subshell_order.get(x[1], 4)))
        
        for _ in range(num_electrons):
            removed = False
            for i, (shell, subshell, count) in enumerate(config_list):
                if count > 0:
                    config_list[i] = (shell, subshell, count - 1)
                    removed = True
                    break
            
            if not removed:
                # No more electrons to remove
                break
        
        # Remove empty entries and re-sort
        config_list[:] = [(s, sub, c) for s, sub, c in config_list if c > 0]
        config_list.sort(key=lambda x: (x[0], subshell_order.get(x[1], 4)))

    def compute_valency(self) -> bool:
        """
        Compute if the atom's valency is satisfied.

        Returns:
            True if valency is satisfied (typically 8 electrons for most atoms,
            or 2 for H/He), False otherwise.
        """
        # Hydrogen and helium: 2 electrons
        if self.symbol in ['H', 'He']:
            return self.electrons_in_valency == 2
        
        # Lithium and beryllium: can be satisfied with 2 electrons (like He)
        # when they have +1 or +2 charge respectively
        if self.symbol.upper() in ['LI', 'BE']:
            if self.electrons_in_valency == 2:
                return True
        
        # Most other elements: 8 electrons (octet rule)
        return self.electrons_in_valency == 8

    def __eq__(self, other):
        """Check equality based on atom_id for graph purposes."""
        if not isinstance(other, BracketAtom):
            return False
        # Use atom_id for graph identity (each atom instance is unique)
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
