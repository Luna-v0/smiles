"""Tracks how this system (our LALR grammar + pysmiles chemistry) is *better*
than pysmiles' own ``read_smiles``.

The value of keeping our parser is that it is a strict, OpenSMILES-compliant
syntax gate.  pysmiles' hand-written parser is permissive and waves through many
malformed strings (unmatched parentheses, dangling/leading bonds, empty input,
out-of-range charge/isotope).  Our pipeline rejects them, so they never reach —
and are never (wrongly) accepted by — the chemistry engine.

These tests assert *our* behaviour (always reject the malformed input) and also
document the gap versus pysmiles so a regression in either is caught.
"""

import logging
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))
logging.getLogger("pysmiles").setLevel(logging.CRITICAL)

from src import validate_smiles
import pysmiles

# Malformed SMILES that our grammar rejects.  The comment is pysmiles' verdict
# at the time of writing; most are *accepted* by pysmiles' parser.
STRICTER_THAN_PYSMILES = [
    "",          # empty string            (pysmiles: accepts)
    " ",         # whitespace              (pysmiles: accepts)
    "=CC",       # leading bond            (pysmiles: accepts)
    "#CC",       # leading bond            (pysmiles: accepts)
    "C=",        # dangling bond           (pysmiles: accepts)
    "C#",        # dangling bond           (pysmiles: accepts)
    "C(",        # unmatched '('           (pysmiles: accepts)
    "C(CC",      # unmatched '('           (pysmiles: accepts)
    "C((C))",    # unmatched '('           (pysmiles: accepts)
    "C()",       # empty branch            (pysmiles: accepts)
    "CC.",       # trailing dot            (pysmiles: accepts)
    "[C+16]",    # charge out of range >15 (pysmiles: rejects)
    "[99999C]",  # isotope out of range    (pysmiles: rejects)
]


def _pysmiles_accepts(smiles: str) -> bool:
    try:
        pysmiles.read_smiles(smiles)
        return True
    except Exception:
        return False


@pytest.mark.parametrize("smiles", STRICTER_THAN_PYSMILES)
def test_our_grammar_rejects_malformed(smiles):
    """Our pipeline rejects malformed SMILES (the syntax gate is ours)."""
    assert validate_smiles(smiles)[0] is False, f"should reject malformed {smiles!r}"


def test_we_catch_more_than_pysmiles_parser():
    """Our pipeline rejects all of these; pysmiles' parser accepts most."""
    ours_reject = sum(validate_smiles(s)[0] is False for s in STRICTER_THAN_PYSMILES)
    pysmiles_accept = sum(_pysmiles_accepts(s) for s in STRICTER_THAN_PYSMILES)
    assert ours_reject == len(STRICTER_THAN_PYSMILES)
    # pysmiles waves through the great majority of these malformed strings.
    assert pysmiles_accept >= 9
