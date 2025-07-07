from dataclasses import field
from typing import List, Optional, Union

from smiles_checker.chem.atomic import Atom, BracketAtom
from smiles_checker.chem.chemistry import chemistry as chem
from smiles_checker.exceptions import ParserException


class ParserManager:
    def __init__(self):
        self.current_open_rnum = []
        self.current_closed_rnum = []
        self.current_chain = []
        self.last_atom = None
        self.ring_atoms = {}

    def __enter__(self):
        """
        Initializes the parser manager.
        """
        self._reset()
        return self

    def __exit__(self):
        self._reset()

    def _reset(self):
        """
        Resets the parser manager to its initial state.
        """
        self.current_open_rnum = []
        self.current_closed_rnum = []
        self.current_chain = []
        self.last_atom = None
        self.ring_atoms = {}

    def chain(self, bond=None, atom=None, rnum=None):
        """
        Handles the 'chain' rule in SMILES parsing, connecting atoms and managing ring closures.

        Args:
            bond (str, optional): The type of bond (e.g., "-", "=", "#"). Defaults to None.
            atom (Atom, optional): The atom object. Defaults to None.
            rnum (int, optional): The ring number. Defaults to None.
        """
        if atom is not None:
            if self.last_atom:
                # Default to single bond if not specified
                actual_bond = bond if bond is not None else "-"
                chem.mol_graph.add_edge(self.last_atom, atom, actual_bond)
            self.last_atom = atom
            return atom
        elif rnum is not None:
            if rnum in self.current_open_rnum:
                ring_opening_atom = self.ring_atoms.pop(rnum)
                if self.last_atom:
                    # Determine bond type for ring closure
                    bond_type = bond
                    if bond is None:
                        bond_type = ":" if self.last_atom.aromatic and ring_opening_atom.aromatic else "-"
                    chem.mol_graph.add_edge(self.last_atom, ring_opening_atom, bond_type)
                self.current_open_rnum.remove(rnum)
                self.current_closed_rnum.append(rnum)
            else:
                if self.last_atom:
                    self.ring_atoms[rnum] = self.last_atom
                self.current_open_rnum.append(rnum)
            return rnum
        else:
            raise ParserException(
                rule="chain",
                parameter=f"bond={bond}, atom={atom}, rnum={rnum}",
                message="Invalid chain rule combination"
            )

    def dot_proxy(self, atom: Atom):
        """
        Handles the '.' (dot) notation in SMILES, indicating disconnected structures.

        Args:
            atom (Atom): The atom following the dot, which starts a new disconnected chain.
        """
        self.last_atom = None  # Disconnect the current chain
        self.current_chain.append(atom) # Add the new atom to the current chain
        return atom

    def inner_branch(self, bond_dot=None, line=None, inner_branch_content=None):
        """
        Handles the parsing of branches in SMILES strings.

        Args:
            bond_dot (str, optional): The bond or dot preceding the branch. Defaults to None.
            line (list): The content of the branch (list of atoms/ring closures).
            inner_branch_content (any, optional): Content of a nested inner branch. Defaults to None.
        """
        # Save the current state before processing the branch
        original_last_atom = self.last_atom
        original_current_chain = list(self.current_chain)

        if bond_dot == ".":
            self.last_atom = None  # Disconnect chain for dot notation

        # Process the main line of the branch
        if line:
            # The 'line' is already processed by the 'chain' method in yacc.py
            # and its effects on last_atom and mol_graph are already applied.
            # We just need to ensure the connection from the original_last_atom
            # to the first atom of the branch if a bond is specified.
            if original_last_atom and line[0] and bond_dot and bond_dot != ".":
                chem.mol_graph.add_edge(original_last_atom, line[0], bond_dot)

        # Process nested inner branches if they exist
        if inner_branch_content:
            # This part assumes inner_branch_content is already processed recursively
            # by the parser and its effects are handled. No explicit action needed here
            # unless there's a specific connection logic for nested branches.
            pass

        # Restore the last_atom to connect the branch back to the main chain
        self.last_atom = original_last_atom
        self.current_chain = original_current_chain

        return line # Return the processed branch content

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

        if not self.current_chain:
            return True

        # Validate the molecular graph structure
        chem.mol_graph.validate_graph()

        # Validate aromaticity
        if not chem.validate_aromacity():
            raise ParserException(
                rule="validate_branch",
                parameter="Aromaticity check failed",
                message="Invalid aromaticity in molecule",
            )

        return True

    def internal_bracket(self, istope=None, symbol=None, chiral=None, hcount=None, charge=None, mol_map=None):
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

    def listify(self, base_element=None, recursion=None):
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
            return [base_element]

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

        atom_obj = chem.Atom(symbol_or_bracket, aromatic=symbol_or_bracket.islower())
        self.current_chain.append(atom_obj)
        self.last_atom = atom_obj
        return atom_obj

    def ring_number(
        self,
        ring_number_or_symbol: str = None,
        ring_number1: Optional[str] = None,
        ring_number2: Optional[str] = None,
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
        rnum = -1

        if ring_number_or_symbol == "%":
            if ring_number1 is None:
                raise ParserException(
                    rule="ring_number",
                    parameter=f"{ring_number_or_symbol} {ring_number1} {ring_number2}",
                    message="Ring number cannot be just '%'"
                )

            digits = ring_number1 + (ring_number2 or "")
            if not digits.isdigit():
                raise ParserException(
                    rule="ring_number",
                    parameter=f"{ring_number_or_symbol} {ring_number1} {ring_number2}",
                    message="Ring number must be a digit or '%'"
                )
            rnum = int(digits)
        elif ring_number_or_symbol.isdigit():
            if ring_number2 is not None:
                raise ParserException(
                    rule="ring_number",
                    parameter=f"{ring_number_or_symbol} {ring_number1} {ring_number2}",
                    message="Ring number cannot have more than one digit after the first"
                )
            digits = ring_number_or_symbol + (ring_number1 or "")
            if not digits.isdigit():
                raise ParserException(
                    rule="ring_number",
                    parameter=f"{ring_number_or_symbol} {ring_number1} {ring_number2}",
                    message="Ring number must be a digit or a number with a leading digit"
                )
            rnum = int(digits)
        else:
            raise ParserException(
                rule="ring_number",
                parameter=f"{ring_number_or_symbol} {ring_number1} {ring_number2}",
                message="Ring number must be a digit or a number with a leading digit"
            )

        if rnum < 1:
            raise ParserException(
                rule="ring_number",
                parameter=f"{ring_number_or_symbol} {ring_number1} {ring_number2}",
                message="Ring number must be greater than 0"
            )
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

    def hcount(self, _, digit: Optional[str] = None) -> int:
        """
        Parses the hydrogen count.
        Args:
            digit: The digit to be parsed.
        Returns:
            The parsed hydrogen count.
        """
        if digit is None:
            return 1
        if not digit.isdigit():
            raise ParserException(
                rule="hcount",
                parameter=digit,
                message="Invalid hydrogen count",
            )
        return int(digit)

    def charge(self, charge1: str = None, charge2: Union[str, None, int] = None) -> int:
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

        if type(charge2) == str and charge2 != charge1:
            raise ParserException(
                rule="charge",
                parameter=f"{charge1} {charge2}",
                message="Charge mismatch",
            )

        if charge2 == "-":
            return -2

        if charge2 == "+":
            return 2

        if charge1 == "-":
            return charge2 * -1

        return charge2

    def chiral(self, chiral1: str = None, chiral2: Optional[str] = None) -> str:
        """
        Fixes the current chiral rotation

        Args:
            chiral1: The first chiral symbol.
            chiral2: The second chiral symbol, if any.
        Returns:
            The current chiral rotation.
        """
        return "counterclockwise" if chiral1 == "@" and chiral2 == "@" else "clockwise"

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
