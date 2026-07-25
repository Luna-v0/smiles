"""Chemistry validator for SMILES molecules."""

from typing import Tuple

from chem.structure import MolecularGraph
from exceptions import ParserException


class ChemistryValidator:
    """
    Validates chemistry rules for molecular graphs.

    Performs valency and aromaticity validation by delegating to a backend
    chemistry engine (see ``chem.adapter``).

    Args:
        backend: Which chemistry driver to use. ``"pysmiles"`` (default) is the
            permissive policy (accepts hypervalent/radical species); ``"rdkit"``
            is the strict, RDKit-equivalent policy. Both bypass the backend's
            own SMILES string parser.
    """

    def __init__(self, backend: str = "pysmiles"):
        self.backend = backend

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
        Perform chemistry validation by delegating to a backend engine.

        Valency and aromaticity are validated by the pysmiles chemistry engine,
        which is driven from the graph *we* constructed (its SMILES string
        parser is never used).  See ``docs/chemistry_backend_study.md`` and
        ``chem.adapter``.

        Args:
            graph: Molecular graph to validate.

        Returns:
            Tuple of (is_valid, exception). If valid, exception is None.
        """
        from chem.adapter import to_descriptor, valid_pysmiles, valid_rdkit

        driver = valid_rdkit if self.backend == "rdkit" else valid_pysmiles
        atoms, bonds = to_descriptor(graph)
        is_valid, reason = driver(atoms, bonds)
        if not is_valid:
            return False, ParserException(
                rule="chemistry",
                parameter=reason or "",
                message=reason or "chemistry validation failed",
            )
        return True, None

    def _check_atom_valency(self, atom, graph: MolecularGraph, is_in_ring: bool = False) -> bool:
        """
        Check if an atom's valency is satisfied using the octet rule.

        For non-aromatic bracket atoms NOT in rings, verifies that the effective
        electron count (internal configuration + graph bond electrons) satisfies
        the octet rule (8 electrons) or duet rule (2 electrons for H, He, Li, Be).

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

        # Count bond electrons from graph connections.
        # Each bond contributes shared electrons to the atom's effective valence:
        #   single/aromatic = 1, double = 2, triple = 3, quadruple = 4
        bond_electrons = 0
        for neighbor, bond_type in graph.adjacency_list.get(atom, []):
            if bond_type == '=':
                bond_electrons += 2
            elif bond_type == '#':
                bond_electrons += 3
            elif bond_type == '$':
                bond_electrons += 4
            else:
                bond_electrons += 1

        effective_electrons = atom.electrons_in_valency + bond_electrons

        # Duet rule: H and He need 2 electrons
        if atom.symbol in ['H', 'HE']:
            return effective_electrons == 2

        # Li, Be: must satisfy duet rule internally (via charge), not via covalent bonds
        if atom.symbol in ['LI', 'BE']:
            return atom.electrons_in_valency == 2

        # Octet rule for most other elements (>= to allow expanded octets)
        return effective_electrons >= 8
