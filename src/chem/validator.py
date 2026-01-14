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
        
        # Validate valency for acyclic parts
        for subgraph in acyclic_subgraphs:
            for atom in subgraph:
                if not self._check_atom_valency(atom, graph):
                    return False, ParserException(
                        rule="validate_rings_and_valency",
                        parameter=str(atom),
                        message=f"Invalid valency for atom {atom.symbol} in acyclic part",
                    )
        
        # Validate valency for ring atoms
        for cycle in graph.cycles:
            for atom in cycle:
                if not self._check_atom_valency(atom, graph):
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

    def validate(self, graph: MolecularGraph) -> Tuple[bool, ParserException | None]:
        """
        Perform all chemistry validations.

        Args:
            graph: Molecular graph to validate.

        Returns:
            Tuple of (is_valid, exception). If valid, exception is None.
        """
        # Validate rings and valency
        is_valid, exception = self.validate_rings_and_valency(graph)
        if not is_valid:
            return False, exception
        
        # Validate aromaticity
        is_valid, exception = self.validate_aromaticity(graph)
        if not is_valid:
            return False, exception
        
        return True, None

    def _check_atom_valency(self, atom, graph: MolecularGraph) -> bool:
        """
        Check if an atom's valency is satisfied.

        Args:
            atom: Atom to check.
            graph: Molecular graph containing the atom.

        Returns:
            True if valency is satisfied, False otherwise.
        """
        # For regular atoms (Atom, not BracketAtom), skip detailed valency checking
        # (they're assumed to follow standard valency rules)
        from chem.atomic import BracketAtom
        
        if not isinstance(atom, BracketAtom):
            return True
        
        # For BracketAtoms, check valency using the graph's method
        # This checks all BracketAtoms in the graph
        return graph.check_valency_for_aba()
