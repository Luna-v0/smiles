import functools
from dataclasses import dataclass, field
from types import prepare_class
from typing import Dict, List, Optional, Union

from chem.atomic import Atom, BracketAtom
from chem.chemistry import chemistry as chem
from chem.graph_builder import GraphBuilder
from exceptions import ParserException


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
        self.graph_builder = GraphBuilder()
        self.last_atom: Optional[Atom] = None
        # Old API compatibility
        self.current_open_rnum: List[int] = []
        # Branch context: pending bond to apply to next atom
        self.pending_branch_bond: Optional[str] = None
        self.branch_start_atom: Optional[Atom] = None
        self.branch_stack: List[Tuple[Atom, Optional[str], bool, bool]] = []  # Stack for nested branches (last_atom, pending_bond, first_branch_processed, connect_next)
        # Track if first atom in branch has been processed (for implicit single bonds)
        self._first_branch_atom_processed: bool = False
        # Track if next atom should connect to last_atom (for atoms following branches)
        self._connect_next_atom: bool = False

    def _get_last_atom_from_chains(self, chains_data):
        """Helper to extract the last atom from a chains structure."""
        if isinstance(chains_data, dict):
            if "atom" in chains_data:
                return chains_data["atom"]
            elif "chains" in chains_data:
                return self._get_last_atom_from_chains(chains_data["chains"])
        return None
    
    def clear(self) -> None:
        """
        Clear the current open and closed cycles.
        """
        self.open_cycles = {}
        self.closed_cycles = set()
        self.graph_builder.clear()
        self.last_atom = None
        self.pending_branch_bond = None
        self.branch_start_atom = None
        self.branch_stack = []
        self._first_branch_atom_processed = False
        self._connect_next_atom = False
    
    def validate(self) -> bool:
        """
        Validate that all rings are closed.
        
        Returns:
            True if valid, raises ParserException if unclosed rings exist.
        """
        # Support old API: current_open_rnum
        if hasattr(self, 'current_open_rnum') and self.current_open_rnum:
            raise ParserException(
                rule="validate",
                parameter="open_cycles",
                message="Unclosed ring numbers"
            )
        if self.has_open_cycles():
            raise ParserException(
                rule="validate",
                parameter="open_cycles",
                message="Unclosed ring numbers"
            )
        return True
    
    def ring_number(self, ring_number_or_symbol: str = None, ring_number1: str = None, ring_number2: str = None) -> int:
        """
        Parse ring number (old API compatibility).
        
        Args:
            ring_number_or_symbol: Ring number as string or "%" for two-digit.
            ring_number1: First digit if using "%".
            ring_number2: Second digit if using "%".
            
        Returns:
            Ring number as integer.
        """
        if ring_number_or_symbol == "%":
            if ring_number1 is None or ring_number2 is None:
                raise ParserException(
                    rule="ring_number",
                    parameter="%",
                    message="'%' must be followed by two digits for ring numbers > 9."
                )
            try:
                rnum = int(ring_number1) * 10 + int(ring_number2)
            except ValueError:
                raise ParserException(
                    rule="ring_number",
                    parameter=f"{ring_number1}{ring_number2}",
                    message="Invalid character for ring number."
                )
        else:
            try:
                rnum = int(ring_number_or_symbol)
            except (ValueError, TypeError):
                raise ParserException(
                    rule="ring_number",
                    parameter=str(ring_number_or_symbol),
                    message="Invalid character for ring number."
                )
        
        # Use rnum method
        return self.rnum(rnum)

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

    def chiral(self, rotation: str = "", chiral1: str = None, chiral2: str = None) -> Optional[str]:
        """
        Function to parse the 'chiral' rule.

        chiral -> "@" | "@@"
        
        Args:
            rotation: Rotation direction (new API).
            chiral1: First chiral symbol (old API compatibility).
            chiral2: Second chiral symbol (old API compatibility).

        """
        # Support old API
        if chiral1 is not None:
            if chiral2 is None or chiral2 == "":
                return "clockwise"
            else:
                return "counter-clockwise"
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

    def hcount(self, hcount: int = 1, _: str = None, digit: str = None) -> int:
        """
        Function to parse the 'hcount' rule.

        hcount -> "H" | "H" digit
        
        Args:
            hcount: Hydrogen count (new API).
            _: Hydrogen symbol "H" (old API compatibility).
            digit: Digit string (old API compatibility).

        """
        # Support old API: hcount(_="H", digit="1")
        if _ is not None and digit is not None:
            try:
                return int(digit)
            except ValueError:
                raise ParserException(
                    rule="hcount",
                    parameter=str(digit),
                    message="Invalid hydrogen count"
                )
        return hcount

    def isotope(self, value=-1) -> int:
        """
        Function to parse the 'isotope' rule.
        isotope -> digit digit digit | digit digit | digit
        """
        return value

    def rnum(self, ring_number=-1, bond_type: str = "-") -> int:
        """
        Function to parse the 'rnum' rule.
        rnum -> digit digit digit | digit digit | digit

        Ring numbers can be reused after being closed (valid in SMILES).

        Returns:
            The ring number.
        """
        cycle_num = ring_number

        if cycle_num in self.open_cycles:
            # Close the ring using graph builder
            if self.last_atom:
                # Use provided bond_type or opening bond
                opening_atom, opening_bond, _ = self.open_cycles[cycle_num]
                bond = bond_type if bond_type != "-" else opening_bond
                closing_atom = self.last_atom
                self.graph_builder.close_ring(cycle_num, closing_atom, bond)
                # Update last_atom to the closing atom (graph_builder.close_ring already does this)
                self.last_atom = self.graph_builder.last_atom
            del self.open_cycles[cycle_num]
            self.closed_cycles.add(cycle_num)
        else:
            # Open a new ring (or reopen a previously closed one)
            # Remove from closed_cycles if reusing a ring number
            if cycle_num in self.closed_cycles:
                self.closed_cycles.discard(cycle_num)

            if self.last_atom:
                self.graph_builder.open_ring(cycle_num, self.last_atom, bond_type)
                # Store reference for parser_manager tracking
                opening_atom, opening_bond, opening_index = self.graph_builder.open_cycles[cycle_num]
                self.open_cycles[cycle_num] = (opening_atom, opening_bond, opening_index)
            else:
                # No last_atom - this shouldn't happen in valid SMILES, but handle gracefully
                # Store the ring number for later when we have an atom
                pass

        return cycle_num

    def _apply_pending_branch_bond(self, new_atom: Atom):
        """Apply pending branch bond if one exists, or implicit single bond for first branch atom."""
        if self.branch_start_atom and not self._first_branch_atom_processed:
            # First atom in branch - connect to branch_start_atom
            bond_type = self.pending_branch_bond if self.pending_branch_bond else "-"
            self.graph_builder.add_bond(
                self.branch_start_atom, new_atom, bond_type
            )
            # Mark that we've processed the first atom in this branch
            self._first_branch_atom_processed = True
            # Clear the pending bond
            self.pending_branch_bond = None

    def atom(self, symbol=None, **kwargs) -> Union[Atom, BracketAtom]:
        """
        Function to parse the 'atom' rule.
        atom -> symbol | bracket_atom

        Args:
            symbol: Atomic symbol (for compatibility with old API).
            **kwargs: Other arguments (bracket_atom, etc.).
        """
        # Helper to finalize atom addition with proper bond connections
        def finalize_atom(atom):
            self.graph_builder.add_atom(atom)
            
            # Check if we need to connect to last_atom (after a branch)
            if self._connect_next_atom and self.last_atom:
                self.graph_builder.add_bond(self.last_atom, atom, "-")
                self._connect_next_atom = False
            else:
                # Normal branch bond handling
                self._apply_pending_branch_bond(atom)
            
            self.last_atom = atom
            return atom
        
        # Support old API: atom("C")
        if symbol is not None and not kwargs:
            try:
                atom = chem.Atom(symbol=symbol)
                return finalize_atom(atom)
            except ParserException as e:
                # Convert error message format for old API compatibility
                if "Invalid Atom Symbol" in e.message:
                    raise ParserException(
                        rule="atom",
                        parameter=symbol,
                        message=f"Invalid Symbol {symbol}"
                    )
                raise

        # New API with kwargs
        if "symbol" in kwargs:
            symbol = kwargs["symbol"]
            atom = chem.Atom(symbol=symbol)
            return finalize_atom(atom)

        match kwargs:
            case {"bracket_atom": bracket_atom}:
                return finalize_atom(bracket_atom)
            case _:
                raise ParserException(
                    rule="atom",
                    parameter=str(kwargs) if kwargs else str(symbol),
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
                # Branch starts with a bond like (=C1)
                # The bond connects the previous atom to the first atom in this branch
                # We need to return the bond info so branch() can use it
                return {"bond": bond_dot, "line": line}
            case {"line": line}:
                return line
            case {
                "bond_dot": str(bond_dot),
                "line": line,
                "inner_branch": inner_branch,
            }:
                # Multiple elements in branch with bond
                return {"bond": bond_dot, "line": line, "inner_branch": inner_branch}
            case {"line": line, "inner_branch": inner_branch}:
                # Multiple lines in branch (no initial bond)
                return {"line": line, "inner_branch": inner_branch}
            case _:
                raise ParserException(
                    rule="inner_branch",
                    parameter=str(kwargs),
                    message="Invalid inner branch rule",
                )

    def start_branch(self):
        """
        Called when '(' is encountered - saves state for branch processing.
        """
        # Save current state to the stack for nested branches (including both flags)
        self.branch_stack.append((self.last_atom, self.pending_branch_bond, self._first_branch_atom_processed, self._connect_next_atom))
        self.branch_start_atom = self.last_atom
        self.pending_branch_bond = None
        self._first_branch_atom_processed = False  # Reset for new branch
        self._connect_next_atom = False  # Reset - atoms inside branch shouldn't use this flag

    def save_branch_bond(self, bond_dot: str) -> str:
        """
        Called when a bond is encountered at the start of a branch.
        Saves the bond to be applied when the first atom in the branch is parsed.

        Returns:
            The bond string (for yacc rule processing).
        """
        self.pending_branch_bond = bond_dot
        return bond_dot

    def end_branch(self, inner_branch):
        """
        Called when ')' is encountered - restores state after branch processing.
        """
        # Restore last_atom to the branch start atom
        # This ensures atoms after the branch connect to where we branched from
        if self.branch_stack:
            saved_last_atom, saved_pending_bond, saved_first_atom_flag, saved_connect_next = self.branch_stack.pop()
            self.last_atom = saved_last_atom
            # Restore the flag state for nested branches
            self._first_branch_atom_processed = saved_first_atom_flag
            # Don't restore _connect_next_atom - we always set it True after branch closes
        elif self.branch_start_atom:
            self.last_atom = self.branch_start_atom

        self.branch_start_atom = None
        self.pending_branch_bond = None
        # Flag that next atom should connect to last_atom (for atoms following branches)
        self._connect_next_atom = True
        return inner_branch

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

        # Validation will be done by ChemistryValidator
        # Just clear state for next molecule
        graph = self.graph_builder.get_graph()
        self.clear()
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
                # rnum() already called from yacc rule, just return the data
                return {"bond": bond, "rnum": int(rnum)}
            case {"atom": atom}:
                return {"atom": atom}
            case {"rnum": int(rnum)}:
                # rnum() already called from yacc rule, just return the data
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
            # Only check for invalid starts if this is the very first chain
            # and there's no preceding atom context
            if chain.get("rnum") and self.last_atom is None:  # starts with a cycle number without atom
                raise ParserException(
                    rule="chains",
                    parameter=str(chain),
                    message="Cannot start with a cycle number.",
                )
                return
            if chain.get("bond") and "atom" not in chain and "rnum" not in chain:  # starts with a bond without atom or rnum
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
                self.graph_builder.add_bond(atom1, atom2, bond)
                self.last_atom = atom2
                return {"chains": {"atom": atom2}}
            case {"chains": {"atom": atom1}, "chain": {"atom": atom2}}:
                self.graph_builder.add_bond(atom1, atom2)
                self.last_atom = atom2
                return {"chains": {"atom": atom2}}
            case {"chains": {"rnum": rnum}, "chain": {"atom": atom}}:
                # Ring was closed, now adding another atom
                # The last_atom should be the atom where the ring closed
                if self.last_atom:
                    self.graph_builder.add_bond(self.last_atom, atom)
                    self.last_atom = atom
                return {"chains": {"atom": atom}}
            case {"chains": {"atom": atom1}, "chain": {"rnum": rnum}}:
                # Atom followed by ring number - ring opens at atom1
                # rnum() already handled opening/closing
                # Keep last_atom as atom1 so next atom connects to it
                self.last_atom = atom1
                return {"chains": {"atom": atom1}}
            case {"chains": {"atom": atom1}, "chain": {"bond": bond, "rnum": rnum}}:
                # Atom followed by bonded ring closure (e.g., c1...=1)
                # The rnum() was already called with the bond type from yacc
                # Just update last_atom
                self.last_atom = atom1
                return {"chains": {"atom": atom1}}
            case {"chains": {"chains": nested_chains}, "chain": chain_data}:
                # Handle nested chains structure
                # Use self.last_atom if available (most reliable), otherwise extract from nested structure
                last_atom_to_use = self.last_atom
                if not last_atom_to_use:
                    # Fallback: extract from nested structure
                    if isinstance(nested_chains, dict):
                        if "atom" in nested_chains:
                            last_atom_to_use = nested_chains["atom"]
                        elif "chains" in nested_chains:
                            # Recursively find last atom
                            last_atom_to_use = self._get_last_atom_from_chains(nested_chains["chains"])
                
                # Process chain_data
                if chain_data.get("atom"):
                    new_atom = chain_data["atom"]
                    if last_atom_to_use and last_atom_to_use != new_atom:
                        # Connect last atom to new atom (avoid self-loops)
                        self.graph_builder.add_bond(last_atom_to_use, new_atom)
                    self.last_atom = new_atom
                    return {"chains": {"atom": new_atom}}
                elif chain_data.get("rnum"):
                    # Ring number - rnum() already called from yacc rule, just update last_atom
                    # Don't call rnum() again to avoid double-closing
                    if last_atom_to_use:
                        self.last_atom = last_atom_to_use
                    return {"chains": nested_chains}
                else:
                    # Unknown chain_data type - try to extract atom and process
                    # This handles cases where chain_data might have unexpected structure
                    # Extract the actual atom from nested_chains to maintain connection
                    if last_atom_to_use:
                        self.last_atom = last_atom_to_use
                    # Return nested structure to preserve state
                    return {"chains": nested_chains}
            case {"chains": {"atom": atom1}, "chain": {"atom": atom2}}:
                self.graph_builder.add_bond(atom1, atom2)
                self.last_atom = atom2
                return {"chains": {"atom": atom2}}
            case _:
                # Debug: print what we got
                debug_kwargs = {}
                for k, v in kwargs.items():
                    try:
                        if isinstance(v, dict):
                            debug_kwargs[k] = {kk: str(vv)[:50] for kk, vv in v.items()}
                        else:
                            debug_kwargs[k] = str(v)[:100]
                    except:
                        debug_kwargs[k] = repr(v)[:100]
                raise ParserException(
                    rule="chains",
                    parameter=str(debug_kwargs),
                    message=f"Invalid chains rule: {debug_kwargs}",
                )
                return

    def internal_bracket(self, internal_bracket=None, **kwargs) -> BracketAtom:
        """
        Function to parse the 'internal_bracket' rule.
        internal_bracket -> isotope? symbol chiral? hcount? charge? mol_map?
        """
        # Support both dict and keyword arguments
        if internal_bracket is None:
            internal_bracket = kwargs
        elif isinstance(internal_bracket, dict):
            internal_bracket.update(kwargs)
        return chem.BracketAtom(**internal_bracket)

    def bracket_atom(self, internal_bracket) -> BracketAtom:
        return internal_bracket

    def chain_branch(self, **kwargs) -> Union[str, List[str]]:
        """
        Function to parse the 'chain_branch' rule.
        chain_branch -> chains | branch | chains chain_branch | branch chain_branch
        """
        # Check if this is a combined chains + branch structure
        if "chains" in kwargs and "branch" in kwargs:
            # We have both chains and branch - need to process atoms after branch
            # The branch may contain nested chains with atoms that need connection
            branch_data = kwargs.get("branch")
            
            # If branch_data has nested chains, those atoms need to be connected
            # after the branch closes (to last_atom which was restored by end_branch)
            if isinstance(branch_data, dict) and "chains" in branch_data:
                # Extract atoms from the nested chains inside branch_data
                nested_chains = branch_data.get("chains")
                self._process_nested_chains_after_branch(nested_chains)
            
            return kwargs
        
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
    
    def _process_nested_chains_after_branch(self, chains_data):
        """
        Process nested chains that appear after a branch.
        These contain atoms that should be connected to last_atom.
        """
        if not chains_data:
            return
        
        if isinstance(chains_data, dict):
            # Look for atoms in the chains structure
            if "atom" in chains_data:
                atom = chains_data["atom"]
                if self.last_atom and self.last_atom != atom:
                    self.graph_builder.add_bond(self.last_atom, atom)
                    self.last_atom = atom
            elif "chains" in chains_data:
                # Recursively process nested chains
                self._process_nested_chains_after_branch(chains_data["chains"])

    def line(self, **kwargs) -> Union[str, List[str]]:
        """
        Function to parse the 'line' rule.
        line -> atom | atom chain_branch | chains
        """
        if kwargs.get("chain_branch") is None:
            # Atom was already processed and last_atom should be set
            atom = kwargs.get("atom")
            if atom and not self.last_atom:
                # Ensure last_atom is set if it wasn't already
                self.last_atom = atom
            return atom if atom else kwargs

        match kwargs:
            case {"atom": atom, "chain_branch": {"chains": {"atom": atom2}}}:
                # Simple case: atom + single atom in chains
                self.graph_builder.add_bond(atom, atom2)
                self.last_atom = atom2
                return {}
            case {"atom": atom, "chain_branch": {"chains": chain_data, "branch": _}}:
                # Chain branch contains both chains and branch(es)
                # Extract first atom from chains to connect to initial atom
                first_atom_in_chains = None
                if isinstance(chain_data, dict):
                    if "atom" in chain_data:
                        first_atom_in_chains = chain_data["atom"]
                    elif "chains" in chain_data:
                        first_atom_in_chains = self._get_last_atom_from_chains(chain_data["chains"])
                
                if atom and first_atom_in_chains and not self.graph_builder.are_bonded(atom, first_atom_in_chains):
                    self.graph_builder.add_bond(atom, first_atom_in_chains)
                    self.last_atom = first_atom_in_chains
                elif atom:
                    self.last_atom = atom
                return {}
            case {"atom": atom, "chain_branch": {"branch": _, "chains": chain_data}}:
                # Same as above but with different key order
                first_atom_in_chains = None
                if isinstance(chain_data, dict):
                    if "atom" in chain_data:
                        first_atom_in_chains = chain_data["atom"]
                    elif "chains" in chain_data:
                        first_atom_in_chains = self._get_last_atom_from_chains(chain_data["chains"])
                
                if atom and first_atom_in_chains and not self.graph_builder.are_bonded(atom, first_atom_in_chains):
                    self.graph_builder.add_bond(atom, first_atom_in_chains)
                    self.last_atom = first_atom_in_chains
                elif atom:
                    self.last_atom = atom
                return {}
            case {"atom": atom, "chain_branch": {"branches_opened": _}}:
                # Only branches, no chains following - atom is already connected
                return {}
            case {"atom": atom, "chain_branch": {"chains": chain_data}}:
                # Handle chains that may contain rnum or other chain elements
                # Extract the first atom from chains to connect to the initial atom
                first_atom_in_chains = None
                if isinstance(chain_data, dict):
                    if "atom" in chain_data:
                        first_atom_in_chains = chain_data["atom"]
                    elif "chains" in chain_data:
                        # Recursively find first atom
                        first_atom_in_chains = self._get_last_atom_from_chains(chain_data["chains"])
                
                # Connect initial atom to first atom in chains if both exist
                if atom and first_atom_in_chains:
                    self.graph_builder.add_bond(atom, first_atom_in_chains)
                    self.last_atom = first_atom_in_chains
                elif atom:
                    # No atom in chains yet, set last_atom to initial atom
                    self.last_atom = atom
                return {}
            case _:
                raise ParserException(
                    rule="line",
                    parameter=str(kwargs),
                    message=f"Invalid line rule: {kwargs}",
                )


parser_manager = ParserManager()
