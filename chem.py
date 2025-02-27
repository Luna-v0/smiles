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

def helium_special_case(elect_config:list[str],charge:int) -> bool:
    
    electrons = sum([int(x[2]) for x in elect_config]) - charge
    
    return True if electrons == 2 else False

def check_valency(elect_config:list[str],charge=0,hcount=0): # elect_config example ['1s2','2s2','2p2']
    if hcount is None:
        hcount = 0
    if charge is None:
        charge = 0
    
    max_valency_per_layer = [2,8,18,32,32,18,8,2]
    
    acc = hcount-charge
    
    unique_elect_config = set([x[0] for x in elect_config])
    
    electron_layers = [
        sum([int(x[2]) for x in electron_layers if elect_config==layer_n]) for layer_n in list(unique_elect_config).sort(reverse=True)
    ] # decrescent list of electrons in each layer
    
    if acc > 0:
        return electron_layers[0] + acc == 8 or (len(electron_layers) == 1 and electron_layers[0]==2)
    
    for x in range(len(electron_layers)):
        electron = electron_layers[x]
        
        if electron > -acc:
            return electron + acc == 8 or (x == len(electron_layers)-1 and electron == 2)
        
        acc +=electron_layers
    
    


def validate_valency_bracket(isotope: Optional[int], symbol: str, chiral: Optional[int],
                         hcount: Optional[int], charge: Optional[int],
                         map: Optional[int]) -> bool:

    elect_config, valency = get_electric_config(symbol)
    
    print(valency,hcount,charge, not hcount)
    
    
    if not hcount and not charge:
        return (valency == 8 or symbol == 'He')
    
    if not charge and hcount and (valency+hcount==8 or (symbol=='H' and hcount==1)): # this or can be removed depending of the standard
        return True

    if not hcount and charge and (valency-charge==8 or valency-charge==0 or helium_special_case(elect_config,charge)):
        return True
        
    
    if hcount and charge and (valency-charge+hcount == 8 or valency-charge+hcount == 0 or helium_special_case(elect_config,hcount-charge)):
        return True
        
        
    if hcount and charge and valency-charge<=0:
        return False
    
    return False