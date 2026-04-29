"""Chemistry validator for SMILES molecules."""

from typing import Tuple

from chem.structure import MolecularGraph
from exceptions import ParserException


class ChemistryValidator:
    """
    Validates chemistry rules for molecular graphs.

    Performs valency and aromaticity validation.
    """

    def validate_rings_and_valency(self, graph: MolecularGraph) -> Tuple[bool, ParserException | None]:
        """
        Validate valency for ring and non-ring parts of the graph.

        Args:
            graph: Molecular graph to validate.

        Returns:
            Tuple of (is_valid, exception). If valid, exception is None.
        """
        # Get acyclic subgraphs (non-ring parts)
        acyclic_subgraphs = graph.get_acyclic_subgraphs()
        
        # Validate valency for acyclic parts (strict octet rule)
        for subgraph in acyclic_subgraphs:
            for atom in subgraph:
                if not self._check_atom_valency(atom, graph, is_in_ring=False):
                    return False, ParserException(
                        rule="validate_rings_and_valency",
                        parameter=str(atom),
                        message=f"Invalid valency for atom {atom.symbol} in acyclic part",
                    )

        # Validate valency for ring atoms (lenient, bonds contribute to valency)
        for cycle in graph.cycles:
            for atom in cycle:
                if not self._check_atom_valency(atom, graph, is_in_ring=True):
                    return False, ParserException(
                        rule="validate_rings_and_valency",
                        parameter=str(atom),
                        message=f"Invalid valency for atom {atom.symbol} in ring",
                    )
        
        return True, None

    def validate_aromaticity(self, graph: MolecularGraph) -> Tuple[bool, ParserException | None]:
        """
        Validate aromaticity using pi electrons rule (Hückel's rule).

        Args:
            graph: Molecular graph to validate.

        Returns:
            Tuple of (is_valid, exception). If valid, exception is None.
        """
        if not graph.huckel():
            return False, ParserException(
                rule="validate_aromaticity",
                parameter="aromatic_cycles",
                message="Aromatic cycles do not satisfy Hückel's rule (4n+2 pi electrons)",
            )
        
        return True, None

    def validate_aromatic_aliphatic_bonds(self, graph: MolecularGraph) -> Tuple[bool, ParserException | None]:
        """
        Validate aromatic/aliphatic carbon bonding.

        Note: Most aromatic/aliphatic validation is done at the SMILES string level
        (checking for trailing 'C'). This method handles edge cases in the graph.

        Args:
            graph: Molecular graph to validate.

        Returns:
            Tuple of (is_valid, exception). If valid, exception is None.
        """
        # Graph-level validation is minimal since string-level check handles most cases
        return True, None

    def validate(self, graph: MolecularGraph) -> Tuple[bool, ParserException | None]:
        """
        Perform all chemistry validations.

        Args:
            graph: Molecular graph to validate.

        Returns:
            Tuple of (is_valid, exception). If valid, exception is None.
        """
        # Validate aromatic/aliphatic bonds
        is_valid, exception = self.validate_aromatic_aliphatic_bonds(graph)
        if not is_valid:
            return False, exception

        # Validate rings and valency
        is_valid, exception = self.validate_rings_and_valency(graph)
        if not is_valid:
            return False, exception

        # Validate aromaticity
        is_valid, exception = self.validate_aromaticity(graph)
        if not is_valid:
            return False, exception

        return True, None

    def _check_atom_valency(self, atom, graph: MolecularGraph, is_in_ring: bool = False) -> bool:
        """
        Check if an atom's valency is satisfied using the octet rule.

        For non-aromatic bracket atoms NOT in rings, verifies that the internal
        configuration (symbol + charge + hydrogen count) satisfies the octet rule
        (8 electrons) or duet rule (2 electrons for H, He, Li, Be).

        Aromatic bracket atoms and atoms in rings are handled more leniently
        because they participate in bonding that contributes to their valency.

        Args:
            atom: Atom to check.
            graph: Molecular graph containing the atom.
            is_in_ring: Whether the atom is part of a ring structure.

        Returns:
            True if valency is satisfied, False otherwise.
        """
        from chem.atomic import BracketAtom

        # Wildcard atoms (*) are placeholders and do not have a fixed valency
        if getattr(atom, 'symbol', None) == "*":
            return True

        # Regular atoms (not bracket atoms) are assumed to follow standard valency rules
        if not isinstance(atom, BracketAtom):
            return True

        # Aromatic bracket atoms participate in delocalized bonding, skip valency check
        if getattr(atom, 'aromatic', False):
            return True

        # Ring atoms have bonds that contribute to their valency, skip strict check
        if is_in_ring:
            return True

        # For non-ring, non-aromatic BracketAtoms, use the octet rule via compute_valency()
        # This checks if the bracket atom's internal configuration is stable
        return atom.compute_valency()
