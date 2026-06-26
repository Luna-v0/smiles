"""SMILES validator public API."""

from chem.graph_builder import GraphBuilder
from chem.structure import MolecularGraph
from chem.validator import ChemistryValidator
from syntax.lex import SmilesLex
from syntax.parser_manager import parser_manager
from syntax.yacc import SmilesParser

parser = SmilesParser()
lexer = SmilesLex()


def parse_smiles(mol: str) -> tuple[bool, MolecularGraph | None, Exception | None]:
    """
    Parse SMILES string and return the generated graph.

    Performs syntax validation and graph generation, but does NOT perform
    chemistry validation.

    Args:
        mol: SMILES string to parse.

    Returns:
        Tuple of (is_valid, graph, exception).
        - is_valid: True if syntax is valid, False otherwise.
        - graph: MolecularGraph instance if valid, None otherwise.
        - exception: Exception if parsing failed, None otherwise.
    """
    try:
        mol = mol.strip()  # tolerate surrounding whitespace (e.g. trailing newline)
        parser.parse(lexer.tokenize(mol))
        if parser_manager.has_open_cycles():
            from exceptions import ParserException
            exc = ParserException(
                rule="parse_smiles",
                parameter=mol,
                message="Unclosed ring number(s)",
            )
            parser_manager.clear()
            return False, None, exc
        graph = parser_manager.get_graph()
        parser_manager.clear()
        return True, graph, None
    except Exception as e:
        parser_manager.clear()
        return False, None, e


def validate_smiles(mol: str) -> tuple[bool, Exception | None]:
    """
    Validate SMILES string with full chemistry validation.

    Performs syntax validation, graph generation, and chemistry validation
    (valency and aromaticity).

    Args:
        mol: SMILES string to validate.

    Returns:
        Tuple of (is_valid, exception).
        - is_valid: True if molecule is valid, False otherwise.
        - exception: Exception if validation failed, None otherwise.
    """
    try:
        # Parse and build graph
        is_valid, graph, parse_exception = parse_smiles(mol)
        if not is_valid:
            return False, parse_exception

        if graph is None:
            return False, Exception("Failed to build graph")

        # Perform chemistry validation
        validator = ChemistryValidator()
        is_valid, chem_exception = validator.validate(graph)

        if not is_valid:
            return False, chem_exception

        return True, None
    except Exception as e:
        return False, e


__all__ = ['parse_smiles', 'validate_smiles', 'MolecularGraph', 'GraphBuilder', 'ChemistryValidator']
