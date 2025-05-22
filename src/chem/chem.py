from dataclasses import dataclass, field
from json import load
from typing import List, Optional


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
        electrons_by_layer: Amount of electons in each layer
        electron_configuration: The electron configuration of the atom
    """

    symbol: str
    valency_layer: int = field(init=False)
    electrons_in_valency: int = field(init=False)
    layers: List[int] = field(init=False)
    electrons_by_layers: List[int] = field(init=False)
    electron_configuration: List[str] = field(default_factory=list)
    aromatic: bool = field(default=False)

    def __post_init__(self):
        """
        Builds valency_layer and electrons_in_valency, also corrects the electron_configuration to list of strings
        """
        if isinstance(self.electron_configuration, str):
            self.electron_configuration = self.electron_configuration.split(" ")

        self.valency_layer = max([int(x[0]) for x in self.electron_configuration])
        self.electrons_in_valency = sum(
            [
                int(x[2])
                for x in self.electron_configuration
                if int(x[0]) == self.valency_layer
            ]
        )

        _layers = {int(x[0]) for x in self.electron_configuration}

        self.layers = list(_layers)
        self.layers.sort(reverse=True)

        self.electrons_by_layers = [
            sum([int(x[2]) for x in self.electron_configuration if x[0] == layer_n])
            for layer_n in self.layers
        ]

    def __eq__(self, other):
        if not isinstance(other, Atom):
            return NotImplemented
        return self.symbol == other.symbol

    def __hash__(self):
        return hash(self.symbol)


@dataclass
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
        """
        Create a Bracket Atom with the given symbol, charge and hidrogens.

        Args:
            isotope: The isotope of the atom
            symbol: The symbol of the atom
            chiral: The chiral of the atom
            hcount: The amount of hydrogens in the atom
            charge: The charge of the atom
            mol_map: The map of the atom
        """
        super().__post_init__()
        self.solo_valency = self.compute_valency()

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
                self.electrons_in_valency = electon + acc
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


class Chem:
    """
    Class for handling Chemistry Domain Logic.

    Attributes:
        organic_atoms: A fixed list of all possible organic atoms
        pt_symbols: All the periodic table symbols
        look_up_table: A look up table for all atoms
    """

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

    def __init__(self, periodic_table_path="src/periodic-table-lookup.json"):

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


chem = Chem()
