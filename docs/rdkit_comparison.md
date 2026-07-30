# OpenSMILES Conformance (and where that differs from RDKit)

The claim this validator makes is **conformance to the OpenSMILES
specification** (<http://opensmiles.org/opensmiles.html>): the accepted
language is the one the spec defines. That is a different — and narrower —
claim than "agreement with RDKit". RDKit is a *reference point* in the
comparisons below, not the ground truth; where its behaviour diverges from
the spec, this validator follows the spec.

The conformance claim is executable:

- `tests/conformance/test_opensmiles_gaps.py` — the measured gap table
  (every divergence found against the pre-conformance `main`), frozen as a
  spec-tagged regression test.
- `tests/conformance/test_differential.py` — a two-way differential test:
  strings generated from the spec BNF must be accepted (*completeness*),
  and mutated strings the parser accepts must be derivable from the spec
  grammar (*soundness*). Both halves are seeded and run in CI.

## Intentional disagreements with RDKit

These two are deliberate and are asserted by dedicated tests; do not "fix"
them toward RDKit.

| SMILES | Ours | RDKit | Why ours is correct |
|---|---|---|---|
| `C-1CCCCC=1` | reject | accept | §3.4: when a ring-closure bond symbol is written on **both** ends, they must match. RDKit silently lets one symbol win. |
| `[C@@TH2](F)(Cl)(Br)I` | accept | reject | §3.8 chirality classes; the conformance plan mandates `@@TH2` parse. RDKit's reader refuses the `@@` + class combination. |

## Unintentional / policy disagreements

These are not spec matters but validation-policy choices, documented in
[validation.md](validation.md):

- **Permissive valence (default backend)**: radicals and
  expanded-octet/hypervalent species (e.g. `N(C)(C)(C)C`, `[SH6]`) are
  accepted; RDKit's default sanitizer rejects them. Only physically
  impossible over-valence is rejected. Opt into RDKit-equivalent behaviour
  with `backend="rdkit"`.
- **Known deviations from the strict BNF** (listed with rationale in
  `tests/conformance/spec_grammar.py` / `test_differential.py`): bare
  single-letter elements outside the organic subset are accepted (E1),
  multi-digit hcount (E2), wide charge digit runs below magnitude 15 (E3);
  conversely, ring closures may not span a dot (`C1.C1` rejected, R1) and
  charge magnitude ≥ 15 is rejected (R2).

## Historical benchmark notes

On ZINC250k the validator accepts 100% (RDKit-canonical SMILES), so
"agreement with RDKit" is saturated there by construction — which is why it
is not the headline metric. The interesting, checkable property is that the
accepted language provably matches the OpenSMILES grammar, which the
differential tests above enforce.
