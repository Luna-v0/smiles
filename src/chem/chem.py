from typing import List, Optional
from json import load
from dataclasses import dataclass, field


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
    

    def __post_init__(self):
        """
        Builds valency_layer and electrons_in_valency, also corrects the electron_configuration to list of strings
        """
        if isinstance(self.electron_configuration, str):
            self.electron_configuration = self.electron_configuration.split(
                ' ')

        self.valency_layer = max([int(x[0])
                                 for x in self.electron_configuration])
        self.electrons_in_valency = sum(
            [int(x[2]) for x in self.electron_configuration if int(x[0]) == self.valency_layer])
        
        
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
    
    def check_valency(self, charge: int = 0, hidrogens: int = 0) -> bool:
        """
        Check if the valency of the current atom would be stable given more charge and hidrogens
        
        Args:
            charge: The amount of charge added to the Atom.
            hidrogens: The amount of Hidrogens added to the Atom.
        
        Returns:
            If that kept the Atom with a stable valency
        """
        
        max_valency_per_layer = [2, 8, 18, 32, 32, 18, 8, 2]
        
        if hidrogens is None: hidrogens = 0
        if charge is None: charge = 0
        
        acc = hidrogens-charge
        
        if sum(self.electrons_by_layers) < -acc:
            return False
        
        if acc < 0:
            for x in range(len(self.electrons_by_layers)):
                electron = self.electrons_by_layers[x]
                
                if electron > -acc:
                    # the octate rule or if it is a noble gas like Helium
                    return electron + acc == 8 or (x == len(self.electrons_by_layers)-1 and electron == 2)

                acc += electron

            # if it just keep removing electrons more than the maximum
            return False

        # first lets try adding electrons to the valency layer
        
        pass
        
        

class Chem:
    """
    Class for handling Chemistry Domain Logic.

    Attributes:
        organic_atoms: A fixed list of all possible organic atoms
        pt_symbols: All the periodic table symbols
        look_up_table: A look up table for all atoms
    """

    def __init__(self, periodic_table_path='../periodic-table-lookup.json'):

        upper_organic_atoms = {"N", "O", "P",
                               "S", "F", "Cl", "Br", "I", "C", "B"}

        with open(periodic_table_path) as JSON:  # loads the periodic table json
            look_up_table_json = dict(load(JSON))

        # set of all organic atoms
        self.organic_atoms = {atom.lower()
                              for atom in upper_organic_atoms} + upper_organic_atoms

        # set of all periodic table symbols
        self.pt_symbols = [look_up_table_json[x]['symbol']
                           for x in look_up_table_json['order']]

        look_up_table_json.pop("order")

        # creates a look up table for all atoms
        self.look_up_table = {x["symbol"]: Atom(
            symbol=x["symbol"],
            electron_configuration=x["electron_configuration"])
            for x in look_up_table_json.values()}



    def check_valency(self,elect_config: list[str], charge=0, hcount=0):
        
        if hcount is None:
            hcount = 0
        if charge is None:
            charge = 0

        max_valency_per_layer = [2, 8, 18, 32, 32, 18, 8, 2]

        acc = hcount-charge

        unique_elect_config = set([x[0] for x in elect_config])

        sorted_elect_config = list(unique_elect_config)
        sorted_elect_config.sort(reverse=True)

        electron_layers = [
            sum([int(x[2]) for x in elect_config if x[0] == layer_n]) for layer_n in sorted_elect_config
        ]  # decrescent list by layer of electrons

        if sum(electron_layers) < -acc:
            return False

        print(electron_layers)

        if acc < 0:
            for x in range(len(electron_layers)):
                electron = electron_layers[x]

                if electron > -acc:
                    return electron + acc == 8 or (x == len(electron_layers)-1 and electron == 2)

                acc += electron

            return False

        for _ in range(len(max_valency_per_layer)):
            electron = electron_layers[0]
            if electron+acc > max_valency_per_layer[len(electron_layers)-1]:
                return electron + acc == 8 or (x == len(electron_layers)-1 and electron == 2)

            acc -= electron
            electron_layers[0] = max_valency_per_layer[len(electron_layers)-1]
            electron_layers.insert(0, 0)

        return False

    def validate_valency_bracket(isotope: Optional[int], symbol: str, chiral: Optional[int],
                                 hcount: Optional[int], charge: Optional[int],
                                 map: Optional[int]) -> bool:

        elect_config, valency = get_electric_config(symbol)

        if not hcount and not charge:
            return (valency == 8 or symbol == 'He')

        return check_valency(elect_config, charge, hcount)
