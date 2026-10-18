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

import re
from dataclasses import dataclass, field
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
        warnings: Notes on accepted-but-irregular input (e.g. a cis/trans
            mark written on a double bond and dropped, or — with
            ``closed_shell=True`` — an atom left open-shell).
        repaired_smiles: For a chemistry failure, the input with only the
            broken atoms rewritten so that it validates (``None`` when no
            unambiguous repair exists). With ``closed_shell=True`` on an
            otherwise-valid molecule, holds the shell-closing repair instead
            (still never changes ``valid``).
        repairs: One line per edited atom, explaining the repair.
        failed_fragment: With ``remove_salts``, the component that failed
            (``None`` when valid or single-component).
        parent_smiles: With ``remove_salts``, the largest component of a
            multi-component input (of the repaired string when one exists).
        removed_fragments: With ``remove_salts``, the other components
            (counter-ions, solvents).
    """

    valid: bool
    tier: Tier | None = None
    message: str | None = None
    position: int | None = None
    expected: set[str] | None = None
    error: Exception | None = None
    warnings: list[str] = field(default_factory=list)
    repaired_smiles: str | None = None
    repairs: list[str] = field(default_factory=list)
    failed_fragment: str | None = None
    parent_smiles: str | None = None
    removed_fragments: list[str] = field(default_factory=list)

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


# A cis/trans mark directly in front of another bond symbol (``C=N\\#N``).
_DIRECTION_ON_BOND = re.compile(r"(?<!\\)[/\\](?=[-=#$:])")


def _direction_warnings(smiles: str) -> list[str]:
    """Warn about cis/trans marks the lexer dropped (they never sit in brackets)."""
    warnings = []
    unbracketed = re.sub(r"\[[^\]]*\]", lambda m: " " * len(m[0]), smiles)
    for match in re.finditer(r"\\\\", unbracketed):
        warnings.append(
            f"doubled backslash at position {match.start()} read as a single '\\' "
            "(a string-escaping artefact)"
        )
    for match in _DIRECTION_ON_BOND.finditer(unbracketed):
        bond = smiles[match.end()]
        warnings.append(
            f"cis/trans mark '{match[0]}' at position {match.start()} is written on a "
            f"'{bond}' bond and was ignored: only single bonds next to a double bond "
            "carry cis/trans direction"
        )
    return warnings


def _chemistry_message(exc: Exception | None) -> str:
    return exc.message if isinstance(exc, ParserException) else str(exc)


def _closed_shell_warning(atoms, n: int, deficit: int) -> str:
    element = atoms[n]["element"]
    plural = "s" if deficit != 1 else ""
    return (
        f"{element} (atom {n}) is open-shell: {deficit} unpaired electron{plural} "
        "(a radical; accepted under the permissive valence policy)"
    )


def validate_smiles_detailed(
    mol: str,
    backend: str = "pysmiles",
    *,
    physiological_ph: bool = True,
    remove_salts: bool = True,
    closed_shell: bool = False,
) -> ValidationResult:
    """Validate a SMILES string, reporting the failing tier on rejection.

    Runs the full pipeline: lexing, grammar, ring-closure semantics and
    chemistry validation (valence + aromaticity via the backend engine).  A
    chemistry failure also gets a repair suggestion when one exists.

    Args:
        mol: SMILES string to validate.
        backend: Chemistry policy — ``"pysmiles"`` (permissive, default) or
            ``"rdkit"`` (strict).
        physiological_ph: When a broken nitrogen can be fixed either by adding
            a charge or by dropping a hydrogen, pick its protonation state at
            body pH (7.4).  Off: keep the written charge and drop the hydrogen.
        remove_salts: Report the parent molecule and the removed counter-ions
            of a multi-component input, and name the component that failed.
            Every component is still validated.
        closed_shell: Opt-in closed-shell check (default off; off is
            byte-identical to this option not existing). When on, every atom
            left with an odd number of unpaired electrons — a radical — gets
            a warning in ``warnings``, plus, when the molecule is otherwise
            valid, a repair suggestion in ``repaired_smiles``/``repairs``
            (add the missing H, or add a negative charge, picked the same
            way ``physiological_ph`` picks for nitrogen). This **never**
            changes ``valid``: the permissive policy accepts radicals by
            design (see docs/validation.md, "Closed-shell check").

    Returns:
        A :class:`ValidationResult`; truthy iff the molecule is valid.

    Raises:
        TypeError: If ``mol`` is not a string.
        ImportError: If the chemistry backend is not installed.
    """
    from chem.adapter import open_shell_atoms, to_descriptor
    from chem.fragments import parent_index, split_fragments
    from chem.repair import suggest_closed_shell_repair, suggest_repair

    graph, failure = _parse(mol)
    if failure is not None:
        return failure
    smiles = mol.strip()
    validator = ChemistryValidator(backend=backend)
    atoms, bonds = to_descriptor(graph)

    result = ValidationResult(valid=True, warnings=_direction_warnings(smiles))
    is_valid, chem_exc = validator.validate_descriptor(atoms, bonds)

    fragments = split_fragments(smiles, atoms, bonds) if remove_salts else None
    if fragments is not None and len(fragments) < 2:
        fragments = None
    if fragments is not None and not is_valid:
        parent = parent_index(fragments)
        for k, fragment in enumerate(fragments):
            ok, exc = validator.validate_descriptor(fragment.atoms, fragment.bonds)
            if not ok:
                result.failed_fragment = fragment.smiles
                if k != parent:
                    chem_exc = ParserException(
                        rule="chemistry",
                        parameter=_chemistry_message(exc),
                        message=f"counter-ion {fragment.smiles!r}: {_chemistry_message(exc)}",
                    )
                break

    if not is_valid:
        result.valid = False
        result.tier = "chemistry"
        result.message = _chemistry_message(chem_exc)
        result.error = chem_exc
        repair = suggest_repair(smiles, atoms, bonds, physiological_ph=physiological_ph)
        if repair is not None:
            repaired_graph, _ = _parse(repair.smiles)
            if repaired_graph is not None:
                repaired = to_descriptor(repaired_graph)
                if validator.validate_descriptor(*repaired)[0]:
                    result.repaired_smiles = repair.smiles
                    result.repairs = repair.edits
                    if fragments is not None:
                        fragments = split_fragments(repair.smiles, *repaired)

    if fragments is not None:
        parent = parent_index(fragments)
        result.parent_smiles = fragments[parent].smiles
        result.removed_fragments = [f.smiles for k, f in enumerate(fragments) if k != parent]

    if closed_shell:
        try:
            flagged = open_shell_atoms(atoms, bonds)
        except Exception:  # noqa: BLE001 - aromaticity failure: nothing to add here
            flagged = []
        for n, deficit in flagged:
            result.warnings.append(_closed_shell_warning(atoms, n, deficit))
        if flagged and is_valid:
            # Only suggest a shell-closing repair when the molecule is
            # otherwise chemistry-valid: if it already failed for an
            # unrelated reason, repaired_smiles must keep meaning "this makes
            # the whole molecule valid", not conflate two unrelated fixes.
            closed_repair = suggest_closed_shell_repair(
                smiles, atoms, bonds, physiological_ph=physiological_ph
            )
            if closed_repair is not None:
                result.repaired_smiles = closed_repair.smiles
                result.repairs = closed_repair.edits

    return result


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
