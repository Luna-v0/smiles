"""Graph builder for constructing molecular graphs during SMILES parsing."""

from typing import Dict, List, Optional, Tuple

from chem.atomic import Atom
from chem.structure import MolecularGraph
from exceptions import ParserException


class GraphBuilder:
    """
    Builds a molecular graph incrementally during SMILES parsing.

    Tracks atoms, bonds, and rings as the parser processes the SMILES string.
    """

    def __init__(self):
        """Initialize an empty graph builder."""
        self.graph = MolecularGraph()
        self.open_cycles: Dict[int, Tuple[Atom, str, int]] = {}  # ring_num -> (opening_atom, bond_type, atom_index)
        self.closed_cycles: set[int] = set()
        self.ring_open_history: Dict[int, int] = {}  # ring_num -> atom_index where it opened (even after closing)
        self.ring_close_history: Dict[int, int] = {}  # ring_num -> atom_index where it closed (for fused ring detection)
        self.last_atom: Optional[Atom] = None
        self.atom_sequence: List[Atom] = []  # Track sequence of atoms as we parse

    def add_atom(self, atom: Atom):
        """
        Add an atom to the graph.

        Args:
            atom: Atom to add.
        """
        if atom not in self.graph.adjacency_list:
            self.graph.adjacency_list[atom] = []
        self.last_atom = atom
        # Track atom sequence for cycle detection
        if atom not in self.atom_sequence:
            self.atom_sequence.append(atom)

    def add_bond(self, atom1: Atom, atom2: Atom, bond_type: str = None):
        """
        Add a bond between two atoms.

        Args:
            atom1: First atom.
            atom2: Second atom.
            bond_type: Type of bond ('-', '=', '#', ':', etc.). If None, automatically
                      determines based on atom aromaticity (aromatic atoms -> ':', else '-').
        """
        # If bond_type not specified, determine automatically
        if bond_type is None:
            # If both atoms are aromatic, use aromatic bond
            if getattr(atom1, 'aromatic', False) and getattr(atom2, 'aromatic', False):
                bond_type = ":"
            else:
                bond_type = "-"

        self.graph.add_edge(atom1, atom2, bond_type=bond_type)
        self.last_atom = atom2

    def are_bonded(self, atom1: Atom, atom2: Atom) -> bool:
        """
        Check if two atoms are already bonded.

        Args:
            atom1: First atom.
            atom2: Second atom.

        Returns:
            True if atoms are bonded, False otherwise.
        """
        return any(neighbor == atom2 for neighbor, _ in self.graph.adjacency_list.get(atom1, []))

    def open_ring(self, ring_number: int, atom: Atom, bond_type: str = "-"):
        """
        Open a ring at a given atom.

        Ring numbers can be reused after being closed (valid in SMILES).

        Args:
            ring_number: Ring number identifier.
            atom: Atom where ring opens.
            bond_type: Bond type for the ring closure.
        """
        # Remove from closed_cycles if reusing a ring number
        if ring_number in self.closed_cycles:
            self.closed_cycles.discard(ring_number)

        # Store opening atom, bond type, and current position in atom sequence
        # The atom_index should be the position where this atom appears in the sequence
        # Since atom is added via add_atom before open_ring is called, it should be in sequence
        # But use last_atom's position to be safe
        if self.last_atom and self.last_atom in self.atom_sequence:
            atom_index = self.atom_sequence.index(self.last_atom)
        elif atom in self.atom_sequence:
            atom_index = self.atom_sequence.index(atom)
        else:
            # Fallback: use current length (atom will be at this index)
            atom_index = len(self.atom_sequence)
        self.open_cycles[ring_number] = (atom, bond_type, atom_index)
        self.ring_open_history[ring_number] = atom_index  # Track even after closing
        self.last_atom = atom

    def close_ring(self, ring_number: int, atom: Atom, bond_type: str = "-"):
        """
        Close a ring at a given atom.

        Args:
            ring_number: Ring number identifier.
            atom: Atom where ring closes.
            bond_type: Bond type for the ring closure.
        """
        # Ring numbers can be reused, so we allow closing a ring number
        # that was previously closed (as long as it's currently open)
        
        if ring_number not in self.open_cycles:
            # Ring opens and closes at same position (self-loop)
            opening_atom = self.last_atom if self.last_atom else atom
            opening_index = len(self.atom_sequence) - 1 if self.atom_sequence else 0
            if opening_atom in self.atom_sequence:
                opening_index = self.atom_sequence.index(opening_atom)
            self.open_cycles[ring_number] = (opening_atom, "-", opening_index)
        
        # Get the opening atom, bond, and index
        opening_atom, opening_bond, opening_index = self.open_cycles[ring_number]
        
        # Recalculate opening_index to ensure it's correct
        if opening_atom in self.atom_sequence:
            opening_index = self.atom_sequence.index(opening_atom)
        
        # Build cycle by finding path from opening to closing atom through the graph
        # IMPORTANT: Find the path BEFORE adding the closing bond
        if opening_atom == atom:
            # Self-loop
            cycle = [opening_atom]
        else:
            # Get closing index
            if self.last_atom and self.last_atom in self.atom_sequence:
                closing_index = self.atom_sequence.index(self.last_atom)
            elif atom in self.atom_sequence:
                closing_index = self.atom_sequence.index(atom)
            else:
                closing_index = len(self.atom_sequence) - 1
            
            # Find rings that opened after this one (check both open and closed rings)
            next_ring_open_index = closing_index + 1
            # Check open cycles
            for other_ring_num, (other_atom, _, other_index) in self.open_cycles.items():
                if other_index > opening_index and other_index < next_ring_open_index:
                    next_ring_open_index = other_index
            # Check closed cycles (using ring_open_history)
            for other_ring_num, other_index in self.ring_open_history.items():
                if other_ring_num != ring_number and other_index > opening_index and other_index < next_ring_open_index:
                    next_ring_open_index = other_index
            
            # Use atoms up to where next ring opened (if any)
            # For fused rings, we include the atom where the next ring opens (shared atom)
            # but not atoms after that (they belong to the other ring)
            if next_ring_open_index <= closing_index:
                # Include only up to the shared atom (where next ring opens)
                end_index = next_ring_open_index
            else:
                end_index = closing_index
            
            allowed_atoms = set()
            if opening_index <= end_index:
                allowed_atoms = set(self.atom_sequence[opening_index:end_index + 1])
            else:
                allowed_atoms = set(self.atom_sequence[opening_index:] + self.atom_sequence[:end_index + 1])
            
            # For fused rings, also include atoms from other rings that closed between opening and closing
            # This allows paths to go through other rings (e.g., naphthalene ring 1 goes through ring 2)
            for other_ring_num, other_close_index in self.ring_close_history.items():
                if other_ring_num != ring_number and opening_index < other_close_index < closing_index:
                    # Include the closing atom of the other ring (shared atom in fused rings)
                    if other_close_index < len(self.atom_sequence):
                        allowed_atoms.add(self.atom_sequence[other_close_index])
            
            # Always include opening and closing atoms
            allowed_atoms.add(opening_atom)
            allowed_atoms.add(atom)

            # Find path using BFS
            from collections import deque
            queue = deque([(opening_atom, [opening_atom])])
            visited = {opening_atom}
            found_path = None

            all_paths = []
            while queue:
                current, path = queue.popleft()

                if current == atom:
                    # Found a path - collect all paths, prefer longer ones
                    all_paths.append(path)
                    continue

                if current in self.graph.adjacency_list:
                    for neighbor, _ in self.graph.adjacency_list[current]:
                        if neighbor not in allowed_atoms:
                            continue
                        if neighbor not in visited or neighbor == atom:
                            visited.add(neighbor)
                            queue.append((neighbor, path + [neighbor]))
            
            # Prefer paths with more atoms (actual cycles), not just direct connections
            # Also prefer paths that avoid atoms that belong exclusively to later rings
            if all_paths:
                # Filter out paths that are just direct connections (2 atoms)
                good_paths = [p for p in all_paths if len(p) > 2]
                if good_paths:
                    # Score paths: prefer shorter paths that stay close to opening_index
                    # But allow paths through shared atoms from other rings (fused rings)
                    def score_path(p):
                        length_penalty = len(p)
                        # Prefer paths that use atoms near opening_index (part of this ring)
                        distance_penalty = 0
                        for a in p:
                            if a in self.atom_sequence:
                                idx = self.atom_sequence.index(a)
                                # Prefer atoms between opening_index and end_index
                                if opening_index <= idx <= end_index:
                                    distance_penalty += abs(idx - opening_index) * 0.1
                                elif idx in [self.ring_close_history.get(rn) for rn in self.ring_close_history if rn != ring_number]:
                                    # Allow atoms where other rings closed (shared atoms in fused rings)
                                    distance_penalty += 5  # Small penalty, not heavy
                                else:
                                    # Heavy penalty for atoms outside the expected range
                                    distance_penalty += 100
                        return length_penalty + distance_penalty
                    found_path = min(good_paths, key=score_path)
                else:
                    # All paths are direct connections - use the sequence fallback instead
                    found_path = None
            
            if found_path and len(found_path) > 2:
                cycle = found_path
            else:
                # Fallback: use atom sequence slice, using end_index calculated above
                # end_index already accounts for rings that opened after this one
                if opening_index <= end_index:
                    cycle = self.atom_sequence[opening_index:end_index + 1]
                else:
                    cycle = self.atom_sequence[opening_index:] + self.atom_sequence[:end_index + 1]

                # For fused rings, include atoms from other rings that closed between opening and closing
                # These are the shared atoms that connect the rings
                atoms_to_insert = []
                for other_ring_num, other_close_index in self.ring_close_history.items():
                    if other_ring_num != ring_number and opening_index < other_close_index < closing_index:
                        if other_close_index < len(self.atom_sequence):
                            other_atom = self.atom_sequence[other_close_index]
                            if other_atom not in cycle:
                                atoms_to_insert.append((other_close_index, other_atom))

                # Insert the atoms at their proper positions in the cycle
                # Sort by index so we insert in the right order
                for insert_index, insert_atom in sorted(atoms_to_insert):
                    # Find the position in the cycle where this atom should go
                    # It should go after the last atom with index < insert_index
                    insert_pos = 0
                    for i, cycle_atom in enumerate(cycle):
                        if cycle_atom in self.atom_sequence:
                            cycle_atom_index = self.atom_sequence.index(cycle_atom)
                            if cycle_atom_index < insert_index:
                                insert_pos = i + 1
                    cycle.insert(insert_pos, insert_atom)

                # Ensure closing atom is included (it might be after end_index)
                if atom not in cycle and atom != opening_atom:
                    cycle.append(atom)
                
                # Ensure cycle starts with opening_atom and ends with closing atom
                if cycle and cycle[0] != opening_atom:
                    if opening_atom in cycle:
                        idx = cycle.index(opening_atom)
                        cycle = cycle[idx:] + cycle[:idx]
                    else:
                        cycle = [opening_atom] + cycle
                
                if cycle and cycle[-1] != atom and atom != opening_atom:
                    if atom in cycle:
                        cycle.remove(atom)
                    cycle.append(atom)
                
                # Remove duplicates while preserving order
                seen = set()
                unique_cycle = []
                for a in cycle:
                    if a not in seen:
                        seen.add(a)
                        unique_cycle.append(a)
                cycle = unique_cycle
        
        # Add the closing bond
        # If bond_type is "-" (default) and both atoms are aromatic, use aromatic bond
        if bond_type == "-" and getattr(opening_atom, 'aromatic', False) and getattr(atom, 'aromatic', False):
            bond_type = ":"
        self.graph.add_edge(opening_atom, atom, bond_type=bond_type)
        
        # Track where this ring closed
        if atom in self.atom_sequence:
            closing_index = self.atom_sequence.index(atom)
            self.ring_close_history[ring_number] = closing_index
        
        # Mark cycle in graph
        self.graph.add_cycle(cycle)
        
        # Move to closed cycles (but keep in ring_open_history)
        del self.open_cycles[ring_number]
        self.closed_cycles.add(ring_number)
        self.last_atom = atom

    def get_graph(self) -> MolecularGraph:
        """
        Get the constructed molecular graph.

        Returns:
            The molecular graph.
        """
        return self.graph

    def clear(self):
        """Clear the graph builder state."""
        self.graph = MolecularGraph()
        self.open_cycles = {}
        self.closed_cycles = set()
        self.ring_open_history = {}
        self.ring_close_history = {}
        self.last_atom = None
        self.atom_sequence = []

    def has_open_cycles(self) -> bool:
        """
        Check if there are any open cycles.

        Returns:
            True if there are open cycles, False otherwise.
        """
        return len(self.open_cycles) > 0

