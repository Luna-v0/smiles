from sly import Lexer
import re
from src.chem.chem import tp_symbols as tp

def generate_regex_from_list(elem_list: list[str]) -> str:
    re_elem = []
    for elem in elem_list:
        re_elem.append(re.escape(elem))

    return '|'.join(re_elem)


def generate_lower(elem_list: list[str]) -> list[str]:
    re_elem = []
    for elem in elem_list:

        re_elem.append(elem.lower())

    return sorted(elem_list,reverse=True) + sorted(re_elem,reverse=True)


atoms = list(set(tp) - set("H"))
bonds = ["=", "#", "$", "/", "\\"]


class SmilesLex(Lexer):
    literals = {".", "@", "-", "+", ":", "%", "H", ")", "(", "]", "[", 'H'}

    tokens = {'semi_bond', 'digit', 'semi_symbol'}

    semi_symbol = rf'{generate_regex_from_list(generate_lower(atoms))}'
    semi_bond = rf'{generate_regex_from_list(bonds)}'
    digit = r'\d'

