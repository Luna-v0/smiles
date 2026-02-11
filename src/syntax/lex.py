import re

from sly import Lexer

from chem.chemistry import pt_symbols

# Aromatic symbols as defined in the SMILES grammar
AROMATIC_SYMBOLS = ['b', 'c', 'n', 'o', 'p', 's', 'se', 'as']

# Get periodic table symbols from JSON (loaded via chemistry module)
pt = pt_symbols


def generate_regex_from_list(elem_list: list[str]) -> str:
    """
    Generate a regex from a list of elements.

    Args:
        elem_list: A list of elements to generate a regex from.
    Returns:
        A regex string that matches any of the elements in the list.
    """
    return "|".join(re.escape(e) for e in elem_list)


def generate_all_symbols(elem_list: list[str]) -> list[str]:
    """
    Generate all symbol variants (both upper and lower case) from a list of elements.
    Symbols are sorted by length (longest first) to ensure proper regex matching.
    
    Args:
        elem_list: A list of elements to generate symbol variants from.
    Returns:
        A sorted list with all elements and their lowercase variants, sorted by length (descending).
    """
    result = []
    for elem in elem_list:
        result.append(elem)
        result.append(elem.lower())
    return sorted(set(result), key=len, reverse=True)


# Pre-compute all symbols (periodic table symbols + aromatic symbols)
# Sort by length (longest first) to ensure longer matches are tried first in regex
# Use only uppercase periodic table symbols (lowercase comes from AROMATIC_SYMBOLS)
pt_uppercase = [s for s in pt if s != "H"]  # Exclude H, it's a literal
_ALL_SYMBOLS = sorted(set(pt_uppercase + AROMATIC_SYMBOLS), key=len, reverse=True)
_SEMI_SYMBOL_PATTERN = generate_regex_from_list(_ALL_SYMBOLS)


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

    literals = {".", "@", "-", "+", ":", "%", "H", ")", "(", "]", "["}

    tokens = {"semi_bond", "digit", "semi_symbol"}

    semi_symbol = _SEMI_SYMBOL_PATTERN
    semi_bond = r'[=#$/\\]'

    @_(r'\d')
    def digit(self, t):
        t.value = int(t.value)
        return t
