"""Regression tests for the unified validation entry point (Phase 0).

The two historical entry points (``src.validate_smiles`` and
``syntax.yacc.validate_smiles``) must agree, internal errors must raise
rather than masquerade as "invalid molecule", and the structured
``ValidationResult`` must preserve the legacy tuple contract.
"""

import builtins
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src import ValidationResult, validate_smiles, validate_smiles_detailed
from syntax.yacc import validate_smiles as validate_smiles_yacc


class TestSingleEntryPoint:
    def test_valid_from_package_path(self):
        assert validate_smiles("CCO") == (True, None)

    def test_valid_from_yacc_path(self):
        assert validate_smiles_yacc("CCO") == (True, None)

    @pytest.mark.parametrize("smiles", ["CCO", "c1ccccc1", "C1CC", "FeCl", "N#N", "C(("])
    def test_entry_points_agree(self, smiles):
        assert validate_smiles(smiles)[0] == validate_smiles_yacc(smiles)[0]


class TestInternalErrorsPropagate:
    def test_non_string_input_raises(self):
        with pytest.raises(TypeError):
            validate_smiles(float("nan"))

    def test_missing_backend_raises_not_rejects(self, monkeypatch):
        """A missing pysmiles install must raise, not report 100% invalid."""
        real_import = builtins.__import__

        def _blocked(name, *args, **kwargs):
            if name.startswith("pysmiles"):
                raise ImportError("pysmiles not installed")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", _blocked)
        with pytest.raises(ImportError):
            validate_smiles("CCO")


class TestValidationResult:
    def test_valid_result(self):
        result = validate_smiles_detailed("CCO")
        assert isinstance(result, ValidationResult)
        assert bool(result) is True
        assert result.tier is None
        assert result.as_tuple() == (True, None)

    def test_lex_tier(self):
        result = validate_smiles_detailed("C?C")
        assert not result
        assert result.tier == "lex"
        assert result.position == 1

    def test_grammar_tier(self):
        result = validate_smiles_detailed("C((")
        assert not result
        assert result.tier == "grammar"

    def test_grammar_tier_position(self):
        result = validate_smiles_detailed("CC)C")
        assert not result
        assert result.tier == "grammar"
        assert result.position == 2

    def test_unclosed_ring_is_ring_semantics_tier(self):
        result = validate_smiles_detailed("C1CCC")
        assert not result
        assert result.tier == "ring_semantics"
        assert "1" in result.message

    def test_as_tuple_preserves_contract(self):
        result = validate_smiles_detailed("C((")
        valid, exc = result.as_tuple()
        assert valid is False
        assert isinstance(exc, Exception)
