import functools
from dataclasses import dataclass, field
from types import prepare_class
from typing import Dict, List, Optional, Union

from smiles_checker.chem.atomic import Atom, BracketAtom
from smiles_checker.chem.chemistry import chemistry as chem


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

    def __init__(self):
        self.open_cycles: Dict[int, List[(Atom, str)]] = {}
        self.closed_cycles: set[int] = set()

    def clear(self) -> None:
        """
        Clear the current open and closed cycles.
        """
        self.open_cycles = {}
        self.closed_cycles = set()

    def has_open_cycles(self) -> bool:
        """
        Check if there are any open cycles.

        Returns:
            bool: True if there are open cycles, False otherwise.
        """
        return len(self.open_cycles) > 0

    def add_atom_to_cycles(self, atom: Atom) -> None:
        for cycle in self.open_cycles.values():
            if atom in cycle:
                raise ParserException(
                    rule="add_atom_to_cycles",
                    parameter=str(atom),
                    message="Atom already in cycle",
                )

            cycle.append(atom)

    def fifteen(self, value=-1) -> int:
        """
        Function to parse the 'fifteen' rule. (Must be at most fifteen)

        fifteen -> digit | digit digit

        """
        if value >= 15:
            raise ParserException(
                rule="fifteen",
                parameter=str(value),
                message="Value must be at most fifteen.",
            )
        return value

    def chiral(self, rotation: str) -> Optional[bool]:
        """
        Function to parse the 'chiral' rule.

        chiral -> "@" | "@@"

        """
        return rotation

    def mol_map(self, value=-1) -> int:
        """
        Function to parse the 'mol_map' rule.

        mol_map -> ":" digit digit digit | ":" digit digit | ":" digit

        """
        return value

    def charge(self, charge: int = 0) -> int:
        """
        Function to parse the 'charge' rule.

        charge -> "-" fifteen | "+" fifteen | "-" | "+" | "-" "-" | "+" "+"

        """
        return charge

    def hcount(self, hcount: int = 1) -> int:
        """
        Function to parse the 'hcount' rule.

        hcount -> "H" | "H" digit

        """
        return hcount

    def isotope(self, **kwargs) -> int:
        """
        Function to parse the 'isotope' rule.
        isotope -> digit digit digit | digit digit | digit
        """
        return self.digit_matching(kwargs, "isotope")

    def rnum(self, **kwargs) -> int:
        """
        Function to parse the 'rnum' rule.
        rnum -> digit digit digit | digit digit | digit
        """
        cycle_num = self.digit_matching(kwargs, "rnum")

        if cycle_num in self.closed_cycles:
            raise ParserException(
                rule="rnum",
                parameter=str(cycle_num),
                message=f"Cycle number {cycle_num} already closed.",
            )

        if cycle_num in self.open_cycles:
            chem.mol_graph.add_cycle(self.open_cycles[cycle_num])
            del self.open_cycles[cycle_num]
            self.closed_cycles.add(cycle_num)

    def atom(self, **kwargs) -> Union[Atom, BracketAtom]:
        """
        Function to parse the 'atom' rule.
        atom -> symbol | bracket_atom
        """
        match kwargs:
            case {"symbol": str(symbol)}:
                return chem.Atom(symbol=symbol)
            case {"bracket_atom": bracket_atom}:
                return bracket_atom
            case _:
                raise ParserException(
                    rule="atom",
                    parameter=str(kwargs),
                    message="Invalid atom rule",
                )

    def semi_bond_rule(self, semi_bond: str) -> str:
        """
        Function to parse the 'semi_bond_rule' rule.
        semi_bond_rule -> " - " | " = "
        """
        ## TODO check if not missing something
        return semi_bond

    def bond(self, semi_bond_rule: str) -> str:
        """
        Function to parse the 'bond' rule.
        bond -> semi_bond_rule | "-"
        """
        return semi_bond_rule

    def bond_dot(self, bond: str) -> str:
        """
        Function to parse the 'bond_dot' rule.
        bond_dot -> bond | "."
        """
        return bond

    def inner_branch(self, **kwargs):
        match kwargs:
            case {"bond_dot": str(bond_dot), "line": line}:
                raise NotImplementedError("Not implemented inner branch with bond dot")
            case {"line": line}:
                return line
            case {
                "bond_dot": str(bond_dot),
                "line": line,
                "inner_branch": inner_branch,
            }:
                raise NotImplementedError(
                    "Not implemented inner branch with bond dot and inner branch"
                )
            case {"line": line, "inner_branch": inner_branch}:
                raise NotImplementedError(
                    "Not implemented inner branch with line and inner branch"
                )
            case _:
                raise ParserException(
                    rule="inner_branch",
                    parameter=str(kwargs),
                    message="Invalid inner branch rule",
                )

    def branch(self, inner_branch):
        return inner_branch

    def symbol(self, semi_symbol: str) -> str:
        """
        Function to parse the 'symbol' rule.
        symbol -> semi_symbol | "H"
        """
        return semi_symbol

    def dot_proxy(self, atom: str) -> str:
        """
        Function to parse the 'dot_proxy' rule.
        dot_proxy -> "." atom
        """
        if self.has_open_cycles():
            raise ParserException(
                rule="dot_proxy",
                parameter=str(atom),
                message="Cannot use dot proxy with open cycles.",
            )

        chem.validate()
        self.clear()
        chem.clear()
        return atom

    def chain(self, **kwargs) -> str:
        """
        Function to parse the 'chain' rule.
        chain -> dot_proxy | bond atom | bond rnum | atom | rnum
        """
        match kwargs:
            case {"dot_proxy": str(dot_proxy)}:
                return {"atom": self.dot_proxy(atom=dot_proxy)}
            case {"bond": str(bond), "atom": atom}:
                return {"bond": bond, "atom": atom}
            case {"bond": str(bond), "rnum": int(rnum)}:
                return {"bond": bond, "rnum": int(rnum)}
            case {"atom": atom}:
                return {"atom": atom}
            case {"rnum": int(rnum)}:
                return {"rnum": int(rnum)}
            case _:
                raise ParserException(
                    rule="chain",
                    parameter=str(kwargs),
                    message="Invalid chain rule",
                )

    def chains(self, **kwargs) -> List[str]:
        """
        Function to parse the 'chains' rule.
        chains -> chain | chain chains
        """
        if kwargs.get("chains") is None:
            chain = kwargs["chain"]
            if chain.get("rnum"):  # starts with a cycle number
                raise ParserException(
                    rule="chains",
                    parameter=str(chain),
                    message="Cannot start with a cycle number.",
                )
                return
            if chain.get("bond"):  # starts with a bond
                raise ParserException(
                    rule="chains",
                    parameter=str(chain),
                    message="Cannot start with a bond.",
                )
                return
            return {"chains": kwargs.get("chain")}

        match kwargs:
            case {"chains": {"rnum": rnum}, "chain": {"rnum": rnum2}}:
                # Double cycle numbers in sequence
                raise ParserException(
                    rule="chains",
                    parameter=str(kwargs),
                    message="Cannot have two cycle numbers in the same chain.",
                )
                return
            case {"chains": {"bond": bond1}, "chain": {"bond": bond2}}:
                # Double bonds in sequence
                raise ParserException(
                    rule="chains",
                    parameter=str(kwargs),
                    message="Cannot have two bonds in the same chain.",
                )
            case {"chains": {"bond": bond, "atom": atom1}, "chain": {"atom": atom2}}:
                chem.mol_graph.add_edge(atom1, atom2, bond=bond)
            case {"chains": {"atom": atom1}, "chain": {"atom": atom2}}:
                chem.mol_graph.add_edge(atom1, atom2)
                return {"chains": {"atom": atom2}}
            case _:
                raise ParserException(
                    rule="chains",
                    parameter=str(kwargs),
                    message="Invalid chains rule",
                )
                return

    def internal_bracket(self, internal_bracket) -> BracketAtom:
        """
        Function to parse the 'internal_bracket' rule.
        internal_bracket -> isotope? symbol chiral? hcount? charge? mol_map?
        """
        return chem.BracketAtom(**internal_bracket)

    def bracket_atom(self, internal_bracket) -> BracketAtom:
        return internal_bracket

    def chain_branch(self, **kwargs) -> Union[str, List[str]]:
        """
        Function to parse the 'chain_branch' rule.
        chain_branch -> chains | branch | chains chain_branch | branch chain_branch
        """
        if kwargs.get("chain_branch") is None:
            if "chains" in kwargs:
                return kwargs

            return {"branches_opened": [kwargs["branch"]]}

        match kwargs:
            case _:
                raise ParserException(
                    rule="chain_branch",
                    parameter=str(kwargs),
                    message="Invalid chain branch rule",
                )

    def line(self, **kwargs) -> Union[str, List[str]]:
        """
        Function to parse the 'line' rule.
        line -> atom | atom chain_branch | chains
        """
        if kwargs.get("chain_branch") is None:
            return kwargs["atom"]

        match kwargs:
            case {"atom": atom, "chain_branch": {"chains": {"atom": atom2}}}:
                chem.mol_graph.add_edge(atom, atom2)
                return {}

            case _:
                raise ParserException(
                    rule="line",
                    parameter=str(kwargs),
                    message="Invalid line rule",
                )


parser_manager = ParserManager()
