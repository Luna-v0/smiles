from typing import Optional
from json import load

upper_organic_atoms = ["N", "O", "P", "S", "F", "Cl", "Br", "I", "C", "B"] 
organic_atoms = [atom.lower() for atom in upper_organic_atoms ] + upper_organic_atoms


def load_table(file:str) -> dict:
    with open(file) as JSON:
        return dict(load(JSON))

look_up_table_json = load_table('periodic-table-lookup.json')

tp_symbols = [look_up_table_json[x]['symbol'] for x in look_up_table_json['order']]

look_up_table_json.pop("order")

look_up_table = {x["symbol"]:x for x in look_up_table_json.values()}


def get_electric_config(symbol:str):
    
    mol_prop = look_up_table[symbol]
    
    electric_config:str = mol_prop['electron_configuration']
    
    splited_elect_config = electric_config.split(' ')
    
    valency_layer = max([int(x[0]) for x in splited_elect_config])
    # if 1s2 2s2 2p6 3s1 (Sodium) for instance the valency layer is 3
    
    elect_in_valency_layer = sum([int(x[2]) for x in splited_elect_config if int(x[0])==valency_layer])
    
    return splited_elect_config, elect_in_valency_layer


def check_valency(elect_config:list[str],charge=0,hcount=0): # elect_config example ['1s2','2s2','2p2']
    if hcount is None:
        hcount = 0
    if charge is None:
        charge = 0 
    
    max_valency_per_layer = [2,8,18,32,32,18,8,2]
    
    acc = hcount-charge
    
    unique_elect_config = set([x[0] for x in elect_config])
    
    sorted_elect_config = list(unique_elect_config)
    sorted_elect_config.sort(reverse=True)
    
    electron_layers = [
        sum([int(x[2]) for x in elect_config if x[0]==layer_n]) for layer_n in sorted_elect_config
    ] # decrescent list by layer of electrons
    
    if sum(electron_layers) < -acc: return False
    
    print(electron_layers)
    
    if acc < 0:
        for x in range(len(electron_layers)):
            electron = electron_layers[x]

            if electron > -acc:
                return electron + acc == 8 or (x == len(electron_layers)-1 and electron == 2)

            acc +=electron

        return False
    
    for _ in range(len(max_valency_per_layer)):
        electron = electron_layers[0]
        if electron+acc > max_valency_per_layer[len(electron_layers)-1]:
            return electron + acc == 8 or (x == len(electron_layers)-1 and electron == 2)

        acc -=electron
        electron_layers[0] = max_valency_per_layer[len(electron_layers)-1]
        electron_layers.insert(0,0)
    
    return False


def validate_valency_bracket(isotope: Optional[int], symbol: str, chiral: Optional[int],
                         hcount: Optional[int], charge: Optional[int],
                         map: Optional[int]) -> bool:

    elect_config, valency = get_electric_config(symbol)
    
    if not hcount and not charge:
        return (valency == 8 or symbol == 'He')
    
    return check_valency(elect_config,charge,hcount)