import copy
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


class BracketLex(Lexer):
    """
    Tokenizer state for the contents of a bracket atom (``[...]``).

    Entered when :class:`SmilesLex` sees ``[`` and exited on ``]``.  Inside
    brackets the full periodic table is legal, ``:`` is the atom-class
    separator (never a bond) and digits belong to isotope / hcount / charge /
    class fields rather than ring closures.

    Attributes:
        tokens: A set of all tokens.
        literals: A set of all literals.
        semi_symbol: A regex for element symbols.
    """

    literals = {"@", "-", "+", ":", "H", "*"}

    tokens = {"digit", "semi_symbol"}

    semi_symbol = _SEMI_SYMBOL_PATTERN

    @_(r'\d')
    def digit(self, t):
        t.value = int(t.value)
        return t

    @_(r'\]')
    def rbracket(self, t):
        """Leave the bracket state and hand ``]`` to the grammar."""
        t.type = ']'
        self.begin(SmilesLex)
        return t

    def tokenize(self, text, lineno=1, index=0):
        """Recover from a prior parse that aborted mid-bracket.

        SLY switches lexer state by reassigning ``self.__class__``, so an
        error inside ``[...]`` strands the instance in this state.  A fresh
        ``tokenize`` call always restarts from the outer state.
        """
        self.begin(SmilesLex)
        return self.tokenize(text, lineno, index)


class SmilesLex(Lexer):
    """
    Tokenizer for SMILES strings (outside bracket atoms).

    ``[`` switches to :class:`BracketLex` until the matching ``]``; bracket
    contents therefore never share token rules with the main chain, where
    only the organic subset is legal bare and digits are ring closures.

    Attributes:
        tokens: A set of all tokens
        literals: A set of all literals
        semi_symbol: A regex for semi symbols
        semi_bond: A regex for semi bonds
        digit: A regex for digits
    """

    literals = {".", "@", "-", "+", ":", "%", "H", ")", "(", "]", "[", "*"}

    tokens = {"semi_bond", "digit", "semi_symbol"}

    semi_symbol = _SEMI_SYMBOL_PATTERN
    semi_bond = r'[=#$/\\]'

    @_(r'\[')
    def lbracket(self, t):
        """Enter the bracket state and hand ``[`` to the grammar."""
        t.type = '['
        self.begin(BracketLex)
        return t

    # Two-character element symbols that are legal as *bare* (non-bracket) atoms.
    # Outside brackets the organic subset only allows Cl and Br (plus the
    # aromatic se/as); every other two-letter periodic symbol (e.g. Cn, Sc, Os,
    # Na) is only legal inside ``[...]``.  The base regex greedily matches the
    # longest periodic-table symbol, which mis-reads e.g. ``Cn1ccnc1`` as
    # Copernicium instead of C + aromatic n, so :meth:`tokenize` splits such
    # tokens back into single atoms when they occur outside brackets.
    _BARE_TWO_CHAR = {"Cl", "Br", "se", "as"}

    @_(r'\d')
    def digit(self, t):
        t.value = int(t.value)
        return t

    def tokenize(self, text, lineno=1, index=0):
        """Tokenize, splitting greedily-merged bare two-letter atoms.

        Wraps the SLY tokenizer so that a two-character ``semi_symbol`` found
        outside brackets — and not one of the genuinely two-character bare
        atoms — is emitted as two single-character atom tokens.  Bracket
        contents (where the full periodic table is valid) are left untouched.
        """
        bracket_depth = 0
        for token in super().tokenize(text, lineno, index):
            if token.type == "[":
                bracket_depth += 1
            elif token.type == "]":
                bracket_depth -= 1

            if (
                bracket_depth == 0
                and token.type == "semi_symbol"
                and len(token.value) == 2
                and token.value not in self._BARE_TWO_CHAR
            ):
                first = copy.copy(token)
                first.value = token.value[0]
                first.end = token.index + 1
                second = copy.copy(token)
                second.value = token.value[1]
                second.index = token.index + 1
                yield first
                yield second
            else:
                yield token
