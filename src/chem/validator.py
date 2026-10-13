"""Chemistry validator for SMILES molecules."""

from typing import TYPE_CHECKING, List, Tuple

from chem.structure import MolecularGraph
from exceptions import ParserException

if TYPE_CHECKING:
    from chem.adapter import AtomDescriptor, BondDescriptor


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
        from chem.adapter import to_descriptor

        return self.validate_descriptor(*to_descriptor(graph))

    def validate_descriptor(
        self, atoms: "List[AtomDescriptor]", bonds: "List[BondDescriptor]"
    ) -> Tuple[bool, ParserException | None]:
        """Validate an ``(atoms, bonds)`` descriptor from ``chem.adapter``.

        Args:
            atoms: Atom descriptors (see ``chem.adapter.to_descriptor``).
            bonds: Bond descriptors.

        Returns:
            Tuple of (is_valid, exception). If valid, exception is None.
        """
        from chem.adapter import valid_pysmiles, valid_rdkit

        driver = valid_rdkit if self.backend == "rdkit" else valid_pysmiles
        is_valid, reason = driver(atoms, bonds)
        if not is_valid:
            return False, ParserException(
                rule="chemistry",
                parameter=reason or "",
                message=reason or "chemistry validation failed",
            )
        return True, None
