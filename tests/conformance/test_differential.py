"""Two-way differential test: the validator vs the OpenSMILES grammar.

*Completeness*: every string generated from the spec BNF (ring-consistent,
see ``opensmiles_generator``) must be accepted by the lex, grammar and
ring-semantics tiers.  Chemistry is a policy layer outside the language, so
a ``tier == "chemistry"`` rejection is not a completeness failure.

*Soundness*: strings mutated toward the validity boundary that the parser
accepts (``parse_smiles`` — syntax + ring semantics, no chemistry) must be
derivable from the spec grammar (with the documented extensions, see
``spec_grammar``).

Known restrictions (validator rejects; strict BNF derives) — these can
never break soundness; the generator avoids producing them:

- **R1**: a dot may not appear while any ring number is open, so ring
  closures cannot span disconnected components (spec allows ``C1.C1``).
- **R2**: charge magnitude is capped below 15 (spec's BNF derives up to
  two digits, e.g. ``[C+16]``).
- **R3**: ring-closure semantics of §3.4/§3.6 (self bonds, duplicate
  bonds, mismatched closure symbols, ring digits after a branch) — these
  are rejections the spec itself mandates, enforced outside the CFG.

Both halves are seeded and deterministic, so they run in CI.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))
from src import parse_smiles, validate_smiles_detailed
from tests.conformance.opensmiles_generator import generate
from tests.conformance.spec_grammar import SpecSmilesAcceptor

SEED = 20260730
N_COMPLETENESS = 500
N_SOUNDNESS_BASES = 150

MUTATION_ALPHABET = "CNOScnos()[]1234567890=#$:/\\.@+-%H*"


def _mutants(rng, base: str):
    """A few single-edit mutations of ``base`` toward the validity boundary."""
    for _ in range(6):
        if not base:
            return
        roll = rng.random()
        position = rng.randrange(len(base))
        char = rng.choice(MUTATION_ALPHABET)
        if roll < 0.4:
            yield base[:position] + char + base[position:]
        elif roll < 0.8:
            yield base[:position] + char + base[position + 1:]
        else:
            yield base[:position] + base[position + 1:]


class TestCompleteness:
    """Spec-derivable strings must never be rejected by the syntax tiers."""

    @pytest.mark.parametrize("smiles", generate(SEED, N_COMPLETENESS))
    def test_generated_string_accepted(self, smiles):
        result = validate_smiles_detailed(smiles)
        assert result.valid or result.tier == "chemistry", (
            f"Completeness bug: spec-derivable {smiles!r} rejected at "
            f"tier={result.tier}: {result.message}"
        )


class TestSoundness:
    """Everything the parser accepts must be derivable from the spec grammar."""

    def test_accepted_mutants_are_spec_derivable(self):
        import random

        rng = random.Random(SEED)
        acceptor = SpecSmilesAcceptor(extensions=True)
        bases = generate(SEED + 1, N_SOUNDNESS_BASES)
        checked = accepted = 0
        violations = []
        for base in bases:
            for mutant in _mutants(rng, base):
                checked += 1
                ok, _, _ = parse_smiles(mutant)
                if not ok:
                    continue
                accepted += 1
                if not acceptor.accepts(mutant):
                    violations.append(mutant)
        assert checked > 500, "mutation harness produced too few cases"
        assert not violations, (
            f"Soundness bugs: parser accepts {len(violations)} strings not "
            f"derivable from the spec grammar (+documented extensions), e.g. "
            f"{violations[:10]!r}"
        )

    def test_generated_strings_are_spec_derivable(self):
        """Sanity: the generator's output itself is in the reference language."""
        acceptor = SpecSmilesAcceptor(extensions=True)
        bad = [s for s in generate(SEED + 2, 200) if not acceptor.accepts(s)]
        assert not bad, f"generator emitted non-derivable strings: {bad[:5]!r}"
