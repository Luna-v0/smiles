from typing import List, Optional, Union
from dataclasses import dataclass, field
import functools

import chem

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
    """
    current_chiral = None
    current_open_rnum = list()
    current_closed_rnum = list()
  
    @fill_none
    def listify(self,base_element, recursion):
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
        if recursion is None: return base_element

        if type(recursion) == list:
            return [base_element] + recursion 

        return [base_element, recursion]


    def atom(self, symbol_or_bracket:str):
        """
        Parses the atom symbol or bracket and from this point on always returns the parser manager
        Args:
            symbol_or_bracket: The atom symbol or bracket atom.
        Returns:
            The atom
        """
        
        if type(symbol_or_bracket) != str \
            or len(symbol_or_bracket) == 1 \
            or symbol_or_bracket in chem.organic_atoms:
                return symbol_or_bracket
        
        elem1, elem2 = symbol_or_bracket
        
        if elem1 in chem.organic_atoms and elem2 in chem.organic_atoms:
            return [elem1,elem2]
        
        raise Exception(f"Inorganic Atom Outside Bracket {rules[0]}")

    
    @fill_none
    def ring_number(self, ring_number_or_symbol:str, ring_number1:Optional[str], ring_number2: Optional[str]) -> int:
        """
        Parses the ring numbers provided.
        Args:
            ring_number_or_symbol: A number or the % symbol.
            ring_number1: The first digit, if any.
            ring_number2: The second digit, if any.
        Returns:
            The parsed ring number as an integer.
        """
        if ring_number_or_symbol == '%':
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
                message="Ring number already closed")
        else:
            self.current_open_rnum.append(rnum)
        
        return rnum
        

    def int(self, digits:List[str]) -> int:
        """
        Parses the provided digits to an integer.
        Args:
            digits: The digits to be parsed.
        Returns:
            The parsed integer.
        """
        return int(''.join(digits))
    
    @fill_none
    def hcount(self, _, digit:Optional[str]) -> int:
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

        if charge2 == '-':
            return -2

        if charge2 == '+':
            return 2

        if charge1 == '-': 
            return charge2 * -1

        return charge2
        
    @fill_none
    def chiral(self, chiral1: str, chiral2: Optional[str]) -> bool:
        """
        Fixes the current chiral rotation

        True for clockwise and False for counterclockwise
        Args:
            chiral1: The first chiral symbol.
            chiral2: The second chiral symbol, if any.
        Returns:
            The current chiral rotation.
        """
        self.current_chiral = chiral2 is None

        return self.current_chiral

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
                    message="Cannot exceed 15")
            return x

        return int(digit1)


parser_manager = ParserManager()
