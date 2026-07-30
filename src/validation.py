"""Single public entry point for SMILES validation.

This module is the one place where the validation pipeline is assembled:
lexing and parsing (``syntax``), ring-closure semantics (``ParserManager``)
and chemistry validation (``ChemistryValidator``).  Both ``src`` and
``syntax.yacc`` re-export :func:`validate_smiles` from here, so the two
historical entry points can no longer disagree.

Failure reporting is tiered: validation *verdicts* (invalid molecules) are
returned as values, while internal errors — a missing chemistry backend, a
non-string input, a bug in the pipeline — raise, so they can never be
mistaken for "invalid molecule".
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from sly.lex import LexError

from chem.structure import MolecularGraph
from chem.validator import ChemistryValidator
from exceptions import ParserException, RingSemanticsException
from syntax.parser_manager import parser_manager
from syntax.yacc import lexer as _lexer
from syntax.yacc import parser as _parser

Tier = Literal["lex", "grammar", "ring_semantics", "chemistry"]

#: Exception types that constitute a validation verdict ("this molecule is
#: invalid").  Anything else that escapes the pipeline is an internal error
#: and propagates to the caller.
_VALIDATION_ERRORS = (ParserException, LexError)


@dataclass
class ValidationResult:
    """Structured outcome of validating one SMILES string.

    Attributes:
        valid: Whether the input is a valid molecule.
        tier: Validation layer that rejected the input (``None`` when valid).
        message: Human-readable reason for rejection (``None`` when valid).
        position: Character index into the input at the failure point, when known.
        expected: Legal tokens at the failure point, when known.
        error: Underlying exception for the legacy tuple contract
            (``None`` when valid).
    """

    valid: bool
    tier: Tier | None = None
    message: str | None = None
    position: int | None = None
    expected: set[str] | None = None
    error: Exception | None = None

    def __bool__(self) -> bool:
        return self.valid

    def as_tuple(self) -> tuple[bool, Exception | None]:
        """Return the legacy ``(is_valid, exception)`` tuple."""
        return self.valid, self.error


def _tier_of(exc: Exception) -> Tier:
    """Classify a validation exception into its pipeline tier."""
    if isinstance(exc, LexError):
        return "lex"
    if isinstance(exc, RingSemanticsException):
        return "ring_semantics"
    return "grammar"


def _failure(exc: Exception) -> ValidationResult:
    """Build a ValidationResult for a rejected input."""
    if isinstance(exc, ParserException):
        return ValidationResult(
            valid=False,
            tier=_tier_of(exc),
            message=exc.message,
            position=exc.position,
            expected=exc.expected,
            error=exc,
        )
    return ValidationResult(
        valid=False,
        tier=_tier_of(exc),
        message=str(exc),
        position=getattr(exc, "error_index", None),
        error=exc,
    )


def _parse(mol: str) -> tuple[MolecularGraph | None, ValidationResult | None]:
    """Lex, parse and ring-check a SMILES string.

    Args:
        mol: SMILES string (surrounding whitespace is tolerated).

    Returns:
        ``(graph, None)`` on success, ``(None, failure)`` on rejection.

    Raises:
        TypeError: If ``mol`` is not a string.
        Exception: Internal pipeline errors propagate unchanged.
    """
    if not isinstance(mol, str):
        raise TypeError(f"SMILES input must be str, got {type(mol).__name__}")
    try:
        _parser.parse(_lexer.tokenize(mol.strip()))
        parser_manager.validate()  # ring semantics: no unclosed ring numbers
        graph = parser_manager.get_graph()
    except _VALIDATION_ERRORS as exc:
        parser_manager.clear()
        return None, _failure(exc)
    except Exception:
        parser_manager.clear()
        raise
    parser_manager.clear()
    return graph, None


def validate_smiles_detailed(mol: str, backend: str = "pysmiles") -> ValidationResult:
    """Validate a SMILES string, reporting the failing tier on rejection.

    Runs the full pipeline: lexing, grammar, ring-closure semantics and
    chemistry validation (valence + aromaticity via the backend engine).

    Args:
        mol: SMILES string to validate.
        backend: Chemistry policy — ``"pysmiles"`` (permissive, default) or
            ``"rdkit"`` (strict).

    Returns:
        A :class:`ValidationResult`; truthy iff the molecule is valid.

    Raises:
        TypeError: If ``mol`` is not a string.
        ImportError: If the chemistry backend is not installed.
    """
    graph, failure = _parse(mol)
    if failure is not None:
        return failure
    is_valid, chem_exc = ChemistryValidator(backend=backend).validate(graph)
    if not is_valid:
        return ValidationResult(
            valid=False,
            tier="chemistry",
            message=chem_exc.message if isinstance(chem_exc, ParserException) else str(chem_exc),
            error=chem_exc,
        )
    return ValidationResult(valid=True)


def validate_smiles(mol: str) -> tuple[bool, Exception | None]:
    """
    Validate SMILES string with full chemistry validation.

    Performs syntax validation, ring-semantics checks, graph generation and
    chemistry validation (valency and aromaticity).

    Args:
        mol: SMILES string to validate.

    Returns:
        Tuple of (is_valid, exception).
        - is_valid: True if molecule is valid, False otherwise.
        - exception: Exception if validation failed, None otherwise.

    Raises:
        TypeError: If ``mol`` is not a string.
        ImportError: If the chemistry backend is not installed.
    """
    return validate_smiles_detailed(mol).as_tuple()


def parse_smiles(mol: str) -> tuple[bool, MolecularGraph | None, Exception | None]:
    """
    Parse SMILES string and return the generated graph.

    Performs syntax validation, ring-semantics checks and graph generation,
    but does NOT perform chemistry validation.

    Args:
        mol: SMILES string to parse.

    Returns:
        Tuple of (is_valid, graph, exception).
        - is_valid: True if syntax is valid, False otherwise.
        - graph: MolecularGraph instance if valid, None otherwise.
        - exception: Exception if parsing failed, None otherwise.

    Raises:
        TypeError: If ``mol`` is not a string.
    """
    graph, failure = _parse(mol)
    if failure is not None:
        return False, None, failure.error
    return True, graph, None


__all__ = [
    "ValidationResult",
    "parse_smiles",
    "validate_smiles",
    "validate_smiles_detailed",
]
