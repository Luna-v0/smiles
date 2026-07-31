"""Regression tests for ring-closure semantics (Phase 3, OpenSMILES §3.4/§3.6).

Ring numbers are reused and matching them needs unbounded state, so these
rules are not expressible in the LALR grammar — the spec itself files them
under semantics (§2.1).  They are enforced by ``ParserManager`` during graph
construction and surface as ``tier == "ring_semantics"``.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))
from src import validate_smiles_detailed


class TestRingSemanticsRejections:
    def test_self_ring_bond(self):
        """§3.4: an atom may not close a ring onto itself."""
        result = validate_smiles_detailed("C11")
        assert not result
        assert result.tier == "ring_semantics"
        assert "1" in result.message and "same atom" in result.message

    def test_duplicate_bond_via_ring_closure(self):
        """§3.4: at most one bond between any atom pair, counting closures."""
        result = validate_smiles_detailed("C12CCCCC12")
        assert not result
        assert result.tier == "ring_semantics"
        assert "2" in result.message and "second bond" in result.message

    def test_duplicate_bond_with_chain_bond(self):
        """§3.4: a ring closure may not duplicate the chain bond either."""
        result = validate_smiles_detailed("C1C1")
        assert not result
        assert result.tier == "ring_semantics"

    def test_ring_bond_order_mismatch(self):
        """§3.4: bond symbols on both ends of a closure must agree.

        RDKit *accepts* C-1CCCCC=1 (it lets the closing symbol win).  The
        disagreement is intentional: OpenSMILES §3.4 says a ring-closure
        bond's symbol, when written on both ends, must be the same, and this
        validator is grounded in the spec, not in RDKit's behaviour.  This is
        the one case where spec-grounding beats the reference implementation.
        """
        result = validate_smiles_detailed("C-1CCCCC=1")
        assert not result
        assert result.tier == "ring_semantics"
        assert "disagree" in result.message and "1" in result.message

    def test_ring_bond_after_branch(self):
        """§3.6: branched_atom ::= atom ringbond* branch* — ring digits
        may not follow a branch on the same atom."""
        result = validate_smiles_detailed("C(O)1CCCCC1")
        assert not result
        assert result.tier == "ring_semantics"
        assert "branch" in result.message

    def test_unclosed_ring(self):
        result = validate_smiles_detailed("C1CCC")
        assert not result
        assert result.tier == "ring_semantics"
        assert "1" in result.message


class TestRingSemanticsAccepted:
    """Valid closure forms that the new checks must not reject."""

    @pytest.mark.parametrize("smiles", [
        "C1CCCCC1",         # plain
        "C-1CCCCC-1",       # explicit symbol on both ends, matching
        "C=1CCCCC=1",       # explicit double on both ends
        "C1CCCCC=1",        # symbol on one end only defines the order
        "C-1CCCCC1",        # symbol on the opening end only
        "C/1CCCCC\\1",      # directional bonds agree (both single)
        "c1ccccc1",         # aromatic closure
        "c:1:c:c:c:c:c:1",  # explicit aromatic bond closures
        "C%12CCCCC%12",     # two-digit ring number
        "C1CC1C1CC1",       # ring number reuse after closing
        "C1(O)CCCC1",       # ring bond *before* branch is legal (§3.6)
        "C1CC(CC1)F",       # closure inside a branch
        "c1ccc2ccccc2c1",   # fused rings
    ])
    def test_valid_closures_accepted(self, smiles):
        result = validate_smiles_detailed(smiles)
        assert result, f"{smiles} should be valid, got {result.tier}: {result.message}"
