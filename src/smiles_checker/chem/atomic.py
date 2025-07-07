from dataclasses import dataclass, field
from typing import List, Optional
import re


def _parse_electron_configuration(electron_config: str) -> frozenset[tuple[tuple[int, str], int]]:
    """
    Parses a single electron configuration string (e.g., "1s2 2s2 2p6")
    into a dictionary mapping (principal_quantum_number, subshell_type) to electron count.
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
    Class for handling only the periodic table atom properties, more properties from
    the periodic table lookup json could be added later.

    Attribute:
        symbol: The symbol of the atom
        valency_layer: The last electron layer number
        electrons_in_valency: Amount of electrons in the valency layer
        layers: Indexes of all layers
        electrons_by_layer: Amount of electrons in each layer
        electron_configuration: The electron configuration of the atom
    """

    symbol: str
    valency_layer: int = field(init=False)
    electrons_in_valency: int = field(init=False)
    layers: tuple[int, ...] = field(init=False)
    electrons_by_layers: tuple[int, ...] = field(init=False)
    electron_configuration: str = field(default_factory=str)
    aromatic: bool = field(default=False)
    _subshell_electrons: frozenset[tuple[tuple[int, str], int]] = field(init=False)

    def __post_init__(self):
        """
        Builds valency_layer and electrons_in_valency.
        """
        
        object.__setattr__(self, '_subshell_electrons', frozenset(_parse_electron_configuration(self.electron_configuration).items()))
        

        if not self._subshell_electrons:  # Handle empty electron configuration
            object.__setattr__(self, 'valency_layer', 0)
            object.__setattr__(self, 'electrons_in_valency', 0)
            object.__setattr__(self, 'layers', ())
            object.__setattr__(self, 'electrons_by_layers', ())
            return

        object.__setattr__(self, 'valency_layer', max([n for (n, _), _ in self._subshell_electrons]))
        object.__setattr__(self, 'electrons_in_valency', sum(
            electrons
            for (n, _), electrons in self._subshell_electrons
            if n == self.valency_layer
        ))

        _layers = {n for (n, _), _ in self._subshell_electrons}

        object.__setattr__(self, 'layers', tuple(_layers))
        object.__setattr__(self, 'electrons_by_layers', tuple([
            sum(electrons for (n, _), electrons in self._subshell_electrons if n == layer_n)
            for layer_n in self.layers
        ]))

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


@dataclass(frozen=True)
class BracketAtom(Atom):
    """

    Class for handling atoms with different properties than the default values of periodic table atoms.

    Attributes:
        isotope: The isotope of the atom
        symbol: The symbol of the atom
        chiral: The chiral of the atom
        hidrogens: The amount of hydrogens in the atom
        charge: The charge of the atom
        mol_map: The map of the atom
    """

    hidrogens: Optional[int] = field(default=None)
    charge: Optional[int] = field(default=None)
    isotope: Optional[int] = field(default=None)
    chiral: Optional[int] = field(default=None)
    mol_map: Optional[int] = field(default=None)

    def __post_init__(self):
        super().__post_init__()
        
        # Adjust electrons based on charge
        if self.charge is not None:
            charge_effect = -self.charge  # Positive charge means removing electrons, negative means adding

            # Create a mutable copy of the subshell electrons
            subshell_electrons = {k: v for k, v in self._subshell_electrons}

            if charge_effect < 0: # Remove electrons
                electrons_to_remove = -charge_effect
                # Sort subshells from outermost to innermost
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

            # Recalculate valency layer and electrons
            if self._subshell_electrons:
                object.__setattr__(self, 'valency_layer', max(n for (n, s), e in self._subshell_electrons))
                object.__setattr__(self, 'electrons_in_valency', sum(e for (n, s), e in self._subshell_electrons if n == self.valency_layer))
            else:
                object.__setattr__(self, 'valency_layer', 0)
                object.__setattr__(self, 'electrons_in_valency', 0)
                object.__setattr__(self, 'layers', ())
                object.__setattr__(self, 'electrons_by_layers', ())

        object.__setattr__(self, 'solo_valency', self.compute_valency())
        
        

    def _next_subshell(self, subshell_type: str) -> str:
        subshell_order = ['s', 'p', 'd', 'f']
        try:
            index = subshell_order.index(subshell_type)
            return subshell_order[index + 1]
        except (ValueError, IndexError):
            return 's' # Start over or handle error

    def _max_electrons_in_subshell(self, subshell_type: str) -> int:
        if subshell_type == 's': return 2
        if subshell_type == 'p': return 6
        if subshell_type == 'd': return 10
        if subshell_type == 'f': return 14
        return 0 # Should not happen

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