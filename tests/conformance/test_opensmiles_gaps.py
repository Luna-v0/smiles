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
