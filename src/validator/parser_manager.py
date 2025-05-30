import functools
from dataclasses import dataclass, field
from typing import List, Optional, Union

from chem import Atom, BracketAtom, chem


def fill_none(func):
    """
    Decorator to fill None values in the function arguments
    """

    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        needed = func.__code__.co_argcount - 1  # minus self
        padded = (list(args) + [None] * needed)[:needed]
        return func(self, *padded, **kwargs)

    return wrapper


@dataclass
class ParserException(Exception):
    """
    Exception for parser errors.

    Args:
        rule: The rule that caused the error.
        parameter: The parameter that caused the error.
        message: The error message.
    """

    rule: str
    parameter: str
    message: str


class ParserManager:
    """
    Parser manager to parse the SMILES strings.

    Attributes:
        current_open_rnum: The current open ring numbers.
        current_closed_rnum: The current closed ring numbers.
        current_chain: The current chain.
    """

    current_open_rnum = list()
    current_closed_rnum = list()
    current_chain: list[Atom] = list()

    def __exit__(self):
        self._reset()

    def _reset(self):
        """
        Resets the parser manager to its initial state.
        """
        self.current_open_rnum = list()
        self.current_closed_rnum = list()
        self.current_chain = list()

    def chain(self, bond=None, atom=None, rnum=None, dot_proxy=None):

        if bond is None:

            if atom is not None:
                return atom

            if rnum is not None:
                return rnum

            return dot_proxy

        if bond == ":" and atom and type(atom) == str and atom[0].isupper():
            raise Exception(
                f"Aromatic bond cannot be use with Uppercase and collon {atom}"
            )

        # TODO: need to check if the atom is not bracketed too

        return [atom, chem.number_of_electrons_per_bond(bond)]

    def inner_branch(self, bond_dot=None, line=None, inner_branch=None):
        """ """
        if bond_dot == ".":
            self.validate_branch()
            self._reset()

        pass

    def validate_branch(self) -> bool:
        """
        Validates based on the current state of the parser manager
        """

        if len(self.current_open_rnum) != 0:
            raise ParserException(
                rule="validate_branch",
                parameter=f"{self.current_open_rnum}",
                message="Unclosed ring numbers",
            )

        starting_aromacity = self.current_chain[0].aromatic
        for atom in self.current_chain[1:]:
            if atom.aromatic != starting_aromacity:
                raise ParserException(
                    rule="validate_branch",
                    parameter=f"{atom}",
                    message="Aromaticity mismatch",
                )

        return True

    @fill_none
    def internal_bracket(self, istope, symbol, chiral, hcount, charge, mol_map):
        """
        Parses the internal bracket and checks for valency.
        """

        br_atom = chem.BracketAtom(
            isotope=istope,
            symbol=symbol,
            chiral=chiral,
            hidrogens=hcount,
            charge=charge,
            mol_map=mol_map,
        )

        self.current_chain.append(br_atom)

        return br_atom

    @fill_none
    def listify(self, base_element, recursion):
        """
        Generic rule for dealing with rules in the following format:

        x -> y x
        x -> y
        Args:
            base_element: Base element.
            recursion: The chain element.
        Returns:
            The parsed atom or chain branch.
        """
        if recursion is None:
            return base_element

        if type(recursion) == list:
            return [base_element] + recursion

        return [base_element, recursion]

    def atom(self, symbol_or_bracket: str):
        """
        Parses the atom symbol or bracket and from this point on always returns the parser manager
        Args:
            symbol_or_bracket: The atom symbol or bracket atom.
        Returns:
            The atom
        """
        if type(symbol_or_bracket) != str:
            return symbol_or_bracket

        if len(symbol_or_bracket) == 1 or symbol_or_bracket in chem.organic_atoms:
            return chem.Atom(symbol_or_bracket, aromatic=True)

        elem1, elem2 = symbol_or_bracket

        if elem1 in chem.organic_atoms and elem2 in chem.organic_atoms:
            return [elem1, elem2]

        raise Exception(
            f"Inorganic Atom Outside Bracket {symbol_or_bracket} not allowed"
        )

    @fill_none
    def ring_number(
        self,
        ring_number_or_symbol: str,
        ring_number1: Optional[str],
        ring_number2: Optional[str],
    ) -> int:
        """
        Parses the ring numbers provided.
        Args:
            ring_number_or_symbol: A number or the % symbol.
            ring_number1: The first digit, if any.
            ring_number2: The second digit, if any.
        Returns:
            The parsed ring number as an integer.
        """
        if ring_number_or_symbol == "%":
            assert (
                ring_number1 is not None and ring_number2 is not None
            ), "Ring number must be a string"
            rnum = self.int([ring_number1, ring_number2])
        else:
            rnum = int(ring_number_or_symbol)

        if rnum in self.current_open_rnum:
            self.current_open_rnum.remove(rnum)
            self.current_closed_rnum.append(rnum)
        elif rnum in self.current_closed_rnum:
            raise ParserException(
                rule="ring_number",
                parameter=f"{rnum}",
                message="Ring number already closed",
            )
        else:
            self.current_open_rnum.append(rnum)

        return rnum

    def int(self, digits: List[str]) -> int:
        """
        Parses the provided digits to an integer.
        Args:
            digits: The digits to be parsed.
        Returns:
            The parsed integer.
        """
        return int("".join(digits))

    @fill_none
    def hcount(self, _, digit: Optional[str]) -> int:
        """
        Parses the hydrogen count.
        Args:
            digit: The digit to be parsed.
        Returns:
            The parsed hydrogen count.
        """
        return int(digit) if digit else 1

    @fill_none
    def charge(self, charge1: str, charge2: Union[str, None, int]) -> int:
        """
        Parsers the charge string to an integer.
        Args:
            charge1: either "+" or "-".
            charge2: either "+", "-", None or an integer.
        Returns:
            The parsed charge as an integer.
        """
        if charge2 is None:
            return 1 if charge1 == "+" else -1

        if charge2 == "-":
            return -2

        if charge2 == "+":
            return 2

        if charge1 == "-":
            return charge2 * -1

        return charge2

    @fill_none
    def chiral(self, chiral1: str, chiral2: Optional[str]) -> str:
        """
        Fixes the current chiral rotation

        Args:
            chiral1: The first chiral symbol.
            chiral2: The second chiral symbol, if any.
        Returns:
            The current chiral rotation.
        """
        return "counterclockwise" if chiral2 else "clockwise"

    @fill_none
    def fifteen(self, digit1: str, digit2: Optional[str]) -> int:
        """
        Fixes fifteen as maximum value for valency
        Args:
            digit1: The first digit to be parsed.
            digit2: The second digit to be parsed, if any.
        Returns:
            The parsed rules.
        """
        if digit2:
            x = int(digit1 + digit2)

            if x > 15:
                raise ParserException(
                    rule="fifteen",
                    parameter=f"{digit1} {digit2}",
                    message="Cannot exceed 15",
                )
            return x

        return int(digit1)


parser_manager = ParserManager()
