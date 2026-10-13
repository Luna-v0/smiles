"""The OpenSMILES conformance gap table, frozen as an executable test.

This file is the conformance claim.  Each row is a divergence between the
validator's accepted language and the OpenSMILES specification
(http://opensmiles.org/opensmiles.html) that was measured on the pre-Phase-2
``main`` and then closed; the table is kept as a permanent regression test,
each row tagged with the spec section that decides it.

Two rows are *intentional disagreements with RDKit* and are marked as such:

- ``[C@@TH2](F)(Cl)(Br)I`` — valid here (spec §3.8 semantics; the plan
  mandates ``@@TH2`` parse), rejected by RDKit.
- ``C-1CCCCC=1`` — invalid here (§3.4: ring-closure bond symbols written on
  both ends must agree), accepted by RDKit (it lets one symbol win).

The claim being made is conformance to OpenSMILES, not agreement with RDKit.

``DOCUMENTED_DEVIATIONS`` is the other half: the divergences that are still
open and documented in ``docs/rdkit_comparison.md`` ("Known deviations from
the strict BNF").  E1-E3 are inputs the grammar accepts beyond the BNF; R1
and R2 are spec-valid inputs it rejects.  Those rows pin *current* behaviour,
not the spec verdict: if one fails, the deviation changed and the doc has to
change with it.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))
from src import validate_smiles_detailed

# (smiles, spec_verdict, spec_section, note)
GAP_TABLE = [
    # -- formerly TOO STRICT: spec-valid, now accepted ----------------------
    ("[CH4:1234]", True, "§3.10", "4-digit atom class"),
    ("C:C", True, "§3.2", "aromatic bond symbol"),
    ("c1:c:c:c:c:c:1", True, "§3.2", "explicit aromatic bonds in a ring"),
    ("[C@TH1](F)(Cl)(Br)I", True, "§3.8", "tetrahedral chirality class"),
    ("[C@@TH2](F)(Cl)(Br)I", True, "§3.8", "INTENTIONAL RDKit disagreement: RDKit rejects"),
    ("[S@AL1](F)(Cl)(Br)I", True, "§3.8", "allenal chirality class"),
    ("[Pt@SP1](F)(Cl)(Br)I", True, "§3.8", "square-planar chirality class"),
    ("[Co@TB15](F)(Cl)(Br)(I)(O)N", True, "§3.8", "trigonal-bipyramidal chirality class"),
    ("[Co@OH25](F)(Cl)(Br)(I)(O)N", True, "§3.8", "octahedral chirality class"),
    # -- formerly TOO LOOSE: spec-invalid, now rejected ---------------------
    ("C11", False, "§3.4", "self ring-bond"),
    ("C12CCCCC12", False, "§3.4", "duplicate bond between one atom pair"),
    ("C-1CCCCC=1", False, "§3.4", "INTENTIONAL RDKit disagreement: RDKit accepts"),
    ("C(O)1CCCCC1", False, "§3.6", "ring bond after branch on the same atom"),
]


@pytest.mark.parametrize(
    "smiles,spec_valid,section,note",
    GAP_TABLE,
    ids=[f"{row[2]}-{row[0]}" for row in GAP_TABLE],
)
def test_gap_table_row(smiles, spec_valid, section, note):
    result = validate_smiles_detailed(smiles)
    assert bool(result) == spec_valid, (
        f"OpenSMILES {section} says {smiles!r} is "
        f"{'valid' if spec_valid else 'invalid'} ({note}); validator said "
        f"{'valid' if result.valid else f'invalid at {result.tier}: {result.message}'}"
    )


def test_spec_invalid_rows_reject_at_ring_semantics_tier():
    """The TOO LOOSE rows are ring-semantics violations, not grammar errors."""
    for smiles, spec_valid, section, _ in GAP_TABLE:
        if spec_valid:
            continue
        result = validate_smiles_detailed(smiles)
        assert result.tier == "ring_semantics", (
            f"{smiles!r} ({section}) should reject at ring_semantics, "
            f"got {result.tier}"
        )


# (deviation, smiles, expected_tier, message_contains, note)
# expected_tier None means the input validates.
DOCUMENTED_DEVIATIONS = [
    # -- E1: bare atoms the BNF only allows inside brackets -------------------
    ("E1", "K", None, None, "bare element outside the organic subset"),
    ("E1", "W", None, None, "bare element outside the organic subset"),
    ("E1", "U", None, None, "bare element outside the organic subset"),
    ("E1", "se", None, None, "bare aromatic se"),
    ("E1", "as", None, None, "bare aromatic as"),
    ("E1", "c1ccsec1", None, None, "bare aromatic se inside a ring"),
    # -- E2: multi-digit hcount; chemistry still applies ----------------------
    ("E2", "[UH12]", None, None, "multi-digit hcount"),
    ("E2", "[CH12]", "chemistry", "valence", "parsed, then rejected on valence"),
    # -- E3: zero-padded charge digit runs below magnitude 15 -----------------
    ("E3", "[C+007]", None, None, "zero-padded charge"),
    ("E3", "[C+014]", None, None, "zero-padded charge just below the bound"),
    ("E3", "[C+0015]", "grammar", None, "the bound is on the value, not the digits"),
    # -- R1: ring closures may not span a dot ---------------------------------
    ("R1", "C1.C1", "grammar", None, "spec-valid ethane"),
    ("R1", "C%10.C%10", "grammar", None, "spec-valid ethane, two-digit ring number"),
    # -- R2: charge magnitude >= 15 is out of range ---------------------------
    ("R2", "[C+15]", "grammar", None, "the spec requires accepting +15"),
    ("R2", "[C-15]", "grammar", None, "the spec requires accepting -15"),
    ("R2", "[C+16]", "grammar", None, "BNF-valid two-digit charge"),
]


@pytest.mark.parametrize(
    "deviation,smiles,expected_tier,message_contains,note",
    DOCUMENTED_DEVIATIONS,
    ids=[f"{row[0]}-{row[1]}" for row in DOCUMENTED_DEVIATIONS],
)
def test_documented_deviation(deviation, smiles, expected_tier, message_contains, note):
    result = validate_smiles_detailed(smiles)
    expected = "valid" if expected_tier is None else f"invalid at {expected_tier}"
    got = "valid" if result.valid else f"invalid at {result.tier}: {result.message}"
    assert result.valid == (expected_tier is None) and result.tier == expected_tier, (
        f"documented deviation {deviation} ({note}) says {smiles!r} is "
        f"{expected}; validator said {got}. Update docs/rdkit_comparison.md "
        "if the change is intended."
    )
    if message_contains is not None:
        assert message_contains in result.message, (
            f"{smiles!r} ({deviation}) should fail on {message_contains}, "
            f"got: {result.message}"
        )
