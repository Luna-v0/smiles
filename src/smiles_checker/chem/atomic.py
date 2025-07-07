from dataclasses import dataclass, field
from typing import List, Optional
import re
from .chemistry import _parse_electron_configuration


@dataclass
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
    layers: List[int] = field(init=False)
    electrons_by_layers: List[int] = field(init=False)
    electron_configuration: List[str] = field(default_factory=list)
    aromatic: bool = field(default=False)
    _subshell_electrons: dict[tuple[int, str], int] = field(init=False)

    def __post_init__(self):
        """
        Builds valency_layer and electrons_in_valency, also corrects the electron_configuration to list of strings
        """
        if isinstance(self.electron_configuration, str):
            self.electron_configuration = self.electron_configuration.split(" ")

        self._subshell_electrons = _parse_electron_configuration(self.electron_configuration)

        self.valency_layer = max([n for n, _ in self._subshell_electrons.keys()])
        self.electrons_in_valency = sum(
            electrons
            for (n, _), electrons in self._subshell_electrons.items()
            if n == self.valency_layer
        )

        _layers = {n for n, _ in self._subshell_electrons.keys()}

        self.layers = list(_layers)
        self.layers.sort(reverse=True)

        self.electrons_by_layers = [
            sum(electrons for (n, _), electrons in self._subshell_electrons.items() if n == layer_n)
            for layer_n in self.layers
        ]

    def get_total_electrons_in_subshell(self, subshell_type: str) -> int:
        """
        Returns the total number of electrons in a given subshell type (s, p, d, f) across all principal quantum numbers.
        """
        return sum(electrons for (_, s), electrons in self._subshell_electrons.items() if s == subshell_type)

    def __eq__(self, other):
        if not isinstance(other, Atom):
            return NotImplemented
        return self.symbol == other.symbol

    def __hash__(self):
        return hash(self.symbol)


@dataclass
class BracketAtom:
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
        if self.hidrogens is None:
            self.hidrogens = 0
        if self.charge is None:
            self.charge = 0

        # Adjust electrons based on charge and hydrogens
        # This is a simplified model and might need refinement for complex cases
        total_charge_effect = self.charge - self.hidrogens # Positive charge means fewer electrons, negative means more

        # Distribute charge effect across subshells, prioritizing outermost
        # This is a very basic model and doesn't account for orbital filling rules
        for (n, subshell), electrons in reversed(list(self._subshell_electrons.items())):
            if total_charge_effect == 0:
                break
            
            if total_charge_effect > 0: # Remove electrons
                removed = min(electrons, total_charge_effect)
                self._subshell_electrons[(n, subshell)] -= removed
                total_charge_effect -= removed
            else: # Add electrons
                added = min(self._max_electrons_in_subshell(subshell) - electrons, -total_charge_effect)
                self._subshell_electrons[(n, subshell)] += added
                total_charge_effect += added

        self.solo_valency = self.compute_valency()

    def _max_electrons_in_subshell(self, subshell_type: str) -> int:
        if subshell_type == 's': return 2
        if subshell_type == 'p': return 6
        if subshell_type == 'd': return 10
        if subshell_type == 'f': return 14
        return 0 # Should not happen

    def _octate_rule(self, layer: int) -> bool:
        """
        Check if the atom is a noble gas or if it is a Helium

        Args:
            layer: The layer to check if it is a noble gas or Helium
        Returns:
            If the atom is a noble gas or Helium
        """
        return self.electrons_by_layers[layer] == 8 or (
            layer == len(self.electrons_by_layers) - 1
            and self.electrons_by_layers[layer] == 2
        )

    def compute_valency(self) -> bool:
        """
        Check if the valency of the current atom would be stable given more charge and hidrogens

        Returns:
            If that kept the Atom with a stable valency
        """
        if self.hidrogens is None:
            self.hidrogens = 0
        if self.charge is None:
            self.charge = 0

        acc = self.hidrogens - self.charge

        if sum(self.electrons_by_layers) < -acc:
            return False

        if acc < 0:
            return self._handle_negative_acc(acc)

        return self._handle_positive_acc(acc)

    def _handle_negative_acc(self, acc: int) -> bool:
        """
        Handle the case where the accumulated charge is negative (removing electrons).

        Args:
            acc: The accumulated charge to be handled.

        Returns:
            If the atom remains stable after removing electrons.
        """
        for x in range(len(self.electrons_by_layers)):
            electron = self.electrons_by_layers[x]

            if electron > -acc:
                return self._octate_rule(x)

            self.electrons_by_layers[x] = 0

            acc += electron

        # filter the layers that are not 0
        self.electrons_by_layers = [x for x in self.electrons_by_layers if x != 0]
        self.valency_layer = len(self.electrons_by_layers) - 1
        self.electrons_in_valency = self.electrons_by_layers[0]

        # if it just keeps removing electrons more than the maximum
        return False

    def _handle_positive_acc(self, acc: int) -> bool:
        """
        Handle the case where the accumulated charge is positive (adding electrons).

        Args:
            acc: The accumulated charge to be handled.

        Returns:
            If the atom remains stable after adding electrons and the amount of electrons.
        """
        max_valency_per_layer = [2, 8, 18, 32, 32, 18, 8, 2]

        # first lets try adding electrons to the valency layer
        if (
            self.electrons_in_valency + acc
            <= max_valency_per_layer[self.valency_layer - 1]
        ):
            return (
                self.electrons_in_valency + acc
                == max_valency_per_layer[self.valency_layer - 1]
            )

        acc -= max_valency_per_layer[self.valency_layer - 1] - self.electrons_in_valency
        self.electrons_in_valency = max_valency_per_layer[self.valency_layer - 1]

        # now we check if there is any layer left to add electrons
        for x in range(len(self.electrons_by_layers)):
            electron = self.electrons_by_layers[x]
            max_electrons_current_layer = max_valency_per_layer[
                len(self.electrons_by_layers) - 1
            ]
            if electron + acc > max_electrons_current_layer:
                # update the valency layer
                self.valency_layer = len(self.electrons_by_layers) - 1
                self.electrons_in_valency = electron + acc
                self.electrons_by_layers.insert(electron + acc, 0)
                return self._octate_rule(x)

            acc -= electron
            self.electrons_by_layers[x] = max_valency_per_layer[
                len(self.electrons_by_layers) - 1
            ]
            self.electrons_by_layers.insert(self.electrons_by_layers[x], 0)

        # if it just keeps adding electrons more than the maximum
        for x in range(self.valency_layer, len(max_valency_per_layer)):
            electron = self.electrons_by_layers[x]
            if electron + acc > max_valency_per_layer[x]:
                # update the valency layer
                self.valency_layer = len(self.electrons_by_layers) - 1
                self.electrons_in_valency = self.electrons_by_layers[0]

                return self._octate_rule(x)

            acc -= electron
            self.electrons_by_layers[x] = max_valency_per_layer[x]
            self.electrons_by_layers.insert(max_valency_per_layer[x], 0)

        return False







