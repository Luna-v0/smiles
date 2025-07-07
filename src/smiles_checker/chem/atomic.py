from dataclasses import dataclass, field
from typing import List, Optional
import re


def _parse_electron_configuration(electron_config: str) -> frozenset[tuple[tuple[int, str], int]]:
    """
    Parses an electron configuration string into a frozenset of (principal_quantum_number, subshell_type), electron_count.

    Args:
        electron_config: The electron configuration string (e.g., "1s2 2s2 2p6").

    Returns:
        A frozenset of tuples, where each tuple contains (principal_quantum_number, subshell_type) and electron_count.
    """
    parsed_config = {}
    # Remove noble gas notations like [He], [Ne], etc.
    cleaned_config = re.sub(r'\[[A-Za-z]+\]\s*', '', electron_config).strip()
    
    if not cleaned_config: # Handle cases where only noble gas notation was present
        return {}
    orbitals = cleaned_config.split()
    for orbital_str in orbitals:
        match = re.match(r'(\d+)([spdf])(\d+)', orbital_str)
        if match:
            n = int(match.group(1))
            subshell = match.group(2)
            electrons = int(match.group(3))
            parsed_config[(n, subshell)] = electrons
    return parsed_config


@dataclass(frozen=True)
class Atom:
    """
    Represents a chemical atom with its electron configuration and related properties.

    Attributes:
        symbol (str): The chemical symbol of the atom.
        valency_layer (int): The principal quantum number of the outermost electron shell.
        electrons_in_valency (int): The number of electrons in the valency_layer.
        layers (tuple[int, ...]): A sorted tuple of all occupied principal quantum numbers.
        electrons_by_layers (tuple[int, ...]): A tuple of electron counts for each corresponding layer.
        electron_configuration (str): The raw electron configuration string.
        aromatic (bool): Indicates if the atom is part of an aromatic system.
    """

    symbol: str
    valency_layer: int = field(init=False)
    electrons_in_valency: int = field(init=False)
    layers: tuple[int, ...] = field(init=False)
    electrons_by_layers: tuple[int, ...] = field(init=False)
    electron_configuration: str = field(default_factory=str)
    aromatic: bool = field(default=False)
    _subshell_electrons: frozenset[tuple[tuple[int, str], int]] = field(init=False)

    def _calculate_electron_properties(self, subshell_electrons_data: frozenset[tuple[tuple[int, str], int]]):
        """
        Calculates and sets valency_layer, electrons_in_valency, layers, and electrons_by_layers.

        Args:
            subshell_electrons_data: A frozenset of (subshell, electron_count) tuples.
        """
        if not subshell_electrons_data:
            object.__setattr__(self, 'valency_layer', 0)
            object.__setattr__(self, 'electrons_in_valency', 0)
            object.__setattr__(self, 'layers', ())
            object.__setattr__(self, 'electrons_by_layers', ())
            return

        object.__setattr__(self, 'valency_layer', max([n for (n, _), _ in subshell_electrons_data]))
        object.__setattr__(self, 'electrons_in_valency', sum(
            electrons
            for (n, _), electrons in subshell_electrons_data
            if n == self.valency_layer
        ))

        _layers = {n for (n, _), _ in subshell_electrons_data}

        object.__setattr__(self, 'layers', tuple(sorted(list(_layers), reverse=True)))
        object.__setattr__(self, 'electrons_by_layers', tuple([
            sum(electrons for (n, _), electrons in subshell_electrons_data if n == layer_n)
            for layer_n in self.layers
        ]))

    

    def __post_init__(self):
        """
        Builds valency_layer and electrons_in_valency.
        """
        
        initial_subshell_electrons = frozenset(_parse_electron_configuration(self.electron_configuration).items())
        object.__setattr__(self, '_subshell_electrons', initial_subshell_electrons)
        
        self._calculate_electron_properties(self._subshell_electrons)

    def get_total_electrons_in_subshell(self, subshell_type: str) -> int:
        """
        Returns the total number of electrons in a given subshell type (s, p, d, f) across all principal quantum numbers.
        """
        return sum(electrons for (n, s), electrons in self._subshell_electrons if s == subshell_type)

    def get_electrons_in_specific_subshell(self, principal_quantum_number: int, subshell_type: str) -> int:
        """
        Returns the number of electrons in a specific subshell, e.g., (3, 'p').
        """
        for (n, s), e in self._subshell_electrons:
            if n == principal_quantum_number and s == subshell_type:
                return e
        return 0

    def _next_subshell(self, subshell_type: str) -> str:
        """
        Determines the next subshell in the standard order (s, p, d, f).

        Args:
            subshell_type: The current subshell type (e.g., 's', 'p').

        Returns:
            The next subshell type, or 's' if at the end of the sequence or invalid input.
        """
        subshell_order = ['s', 'p', 'd', 'f']
        try:
            index = subshell_order.index(subshell_type)
            return subshell_order[index + 1]
        except (ValueError, IndexError):
            return 's' # Start over or handle error

    def _max_electrons_in_subshell(self, subshell_type: str) -> int:
        """
        Returns the maximum number of electrons a given subshell type can hold.

        Args:
            subshell_type: The subshell type (s, p, d, f).

        Returns:
            The maximum electron capacity of the subshell.
        """
        if subshell_type == 's': return 2
        if subshell_type == 'p': return 6
        if subshell_type == 'd': return 10
        if subshell_type == 'f': return 14
        return 0 # Should not happen


@dataclass(frozen=True)
class BracketAtom(Atom):
    """
    Represents an atom within square brackets in SMILES, allowing for explicit specification of properties.

    Inherits from Atom and adds properties like isotope, chiral, hydrogen count, charge, and mapping.

    Attributes:
        hidrogens (Optional[int]): Number of implicit hydrogens attached to the atom.
        charge (Optional[int]): The charge of the atom.
        isotope (Optional[int]): The isotope number of the atom.
        chiral (Optional[int]): Chiral specification (e.g., @ or @@).
        mol_map (Optional[int]): Atom mapping number.
    """

    hidrogens: Optional[int] = field(default=None)
    charge: Optional[int] = field(default=None)
    isotope: Optional[int] = field(default=None)
    chiral: Optional[int] = field(default=None)
    mol_map: Optional[int] = field(default=None)

    def __post_init__(self):
        super().__post_init__()
        
        if self.charge is not None:
            charge_effect = -self.charge  # Positive charge means removing electrons, negative means adding

            # Create a mutable dictionary from the frozenset for manipulation
            subshell_electrons = {k: v for k, v in self._subshell_electrons}

            if charge_effect < 0: # Remove electrons
                electrons_to_remove = -charge_effect
                # Sort subshells from outermost to innermost (higher n, then higher subshell type)
                sorted_subshells = sorted([k for k, v in subshell_electrons.items()], key=lambda x: (x[0], {'s': 0, 'p': 1, 'd': 2, 'f': 3}[x[1]]), reverse=True)

                for n, s in sorted_subshells:
                    if electrons_to_remove == 0:
                        break
                    available_electrons = subshell_electrons[(n, s)]
                    removed = min(electrons_to_remove, available_electrons)
                    subshell_electrons[(n, s)] -= removed
                    if subshell_electrons[(n, s)] == 0:
                        del subshell_electrons[(n, s)]
                    electrons_to_remove -= removed

            elif charge_effect > 0: # Add electrons
                electrons_to_add = charge_effect
                # Sort subshells from innermost to outermost to fill them in order
                sorted_subshells = sorted([k for k, v in subshell_electrons.items()], key=lambda x: (x[0], {'s': 0, 'p': 1, 'd': 2, 'f': 3}[x[1]]))
                
                # Fill existing subshells first
                for n, s in sorted_subshells:
                    if electrons_to_add == 0:
                        break
                    capacity = self._max_electrons_in_subshell(s)
                    available_space = capacity - subshell_electrons.get((n, s), 0)
                    add = min(electrons_to_add, available_space)
                    subshell_electrons[(n, s)] = subshell_electrons.get((n, s), 0) + add
                    electrons_to_add -= add

                # If electrons still remain, add new subshells
                if electrons_to_add > 0:
                    if sorted_subshells:
                        last_n, last_s = sorted_subshells[-1]
                    else: # Case where atom had no electrons to begin with
                        last_n, last_s = 1, 's' 
                        subshell_electrons[(1,'s')] = 0

                    current_n, current_s = last_n, last_s

                    while electrons_to_add > 0:
                        current_n += 1
                        current_s = 's'
                        capacity = self._max_electrons_in_subshell(current_s)
                        add = min(electrons_to_add, capacity)
                        subshell_electrons[(current_n, current_s)] = add
                        electrons_to_add -= add

            object.__setattr__(self, '_subshell_electrons', frozenset(subshell_electrons.items()))

            self._calculate_electron_properties(self._subshell_electrons)

        object.__setattr__(self, 'solo_valency', self.compute_valency())
        
        

    

    def _octate_rule(self, electrons: int) -> bool:
        """
        Check if the given number of electrons satisfies the octet rule (8 valence electrons) or duet rule (2 valence electrons for H/He).
        """
        if self.valency_layer == 1:
            return electrons == 2
        return electrons == 8

    def compute_valency(self, bonds: int = 0) -> bool:
        """
        Check if the valency of the current atom is stable, considering hydrogens, charge, and bonds.
        """
        # Calculate effective valence electrons: initial valence electrons + hydrogens (as bonds)
        effective_valence_electrons = self.electrons_in_valency + (self.hidrogens if self.hidrogens is not None else 0) + bonds

        # Check against octet/duet rule
        return self._octate_rule(effective_valence_electrons)

    def __eq__(self, other):
        if not isinstance(other, BracketAtom):
            return NotImplemented
        return (
            super().__eq__(other)
            and self.hidrogens == other.hidrogens
            and self.charge == other.charge
            and self.isotope == other.isotope
            and self.chiral == other.chiral
            and self.mol_map == other.mol_map
        )