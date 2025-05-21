import re

from sly import Lexer

from src.chem.chem import chem

pt = chem.pt_symbols


def generate_regex_from_list(elem_list: list[str]) -> str:
    """
    Generate a regex from a list of elements.

    Args:
        elem_list: A list of elements to generate a regex from.
    Returns:
        A regex string that matches any of the elements in the list.
    """
    re_elem = []
    for elem in elem_list:
        re_elem.append(re.escape(elem))

    return "|".join(re_elem)


def generate_lower(elem_list: list[str]) -> list[str]:
    """
    Generate a list of lower case elements from a list of elements and return both the lower and upper case elements.

    Args:
        elem_list: A list of elements to generate a lower case list from.
    Returns:
        A list with all elements.
    """
    re_elem = []
    for elem in elem_list:

        re_elem.append(elem.lower())

    return sorted(elem_list, reverse=True) + sorted(re_elem, reverse=True)


atoms = list(set(pt) - set("H"))
bonds = ["=", "#", "$", "/", "\\"]


class SmilesLex(Lexer):
    """
    Tokenizer for SMILES strings.

    Attributes:
        tokens: A set of all tokens
        literals: A set of all literals
        semi_symbol: A regex for semi symbols
        semi_bond: A regex for semi bonds
        digit: A regex for digits
    """

    literals = {".", "@", "-", "+", ":", "%", "H", ")", "(", "]", "[", "H"}

    tokens = {"semi_bond", "digit", "semi_symbol"}

    semi_symbol = rf"{generate_regex_from_list(generate_lower(atoms))}"
    semi_bond = rf"{generate_regex_from_list(bonds)}"
    digit = r"\d"
