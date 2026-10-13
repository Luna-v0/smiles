# OpenSMILES Conformance (and where that differs from RDKit)

The claim this validator makes is **conformance to the OpenSMILES
specification** (<http://opensmiles.org/opensmiles.html>): the accepted
language is the one the spec defines. That is a different — and narrower —
claim than "agreement with RDKit". RDKit is a *reference point* in the
comparisons below, not the ground truth; where its behaviour diverges from
the spec, this validator follows the spec.

The conformance claim is backed by example-based regression tests, run with
`pytest` (the CI workflow only builds and deploys these docs):

- `tests/conformance/test_opensmiles_gaps.py` — the measured gap table
  (every divergence found against the pre-conformance `main`), frozen as a
  spec-tagged regression test. Its spec-invalid rows must also be rejected
  by the `ring_semantics` tier, not by the grammar. Its
  `DOCUMENTED_DEVIATIONS` table pins the deviations that are still open
  (see [Unintentional / policy disagreements](#unintentional-policy-disagreements)).
- `tests/syntax/test_ring_semantics.py` — the §3.4/§3.6 ring-closure rules
  (self bonds, duplicate bonds, mismatched closure symbols, ring digits
  after a branch, unclosed rings) and the closures that must stay accepted.
- `tests/chem/test_vs_pysmiles.py` — malformed strings the grammar must
  reject (leading or dangling bonds, unmatched parentheses, empty branches,
  charge magnitude ≥ 15).

These tests pin known cases. Nothing checks the accepted language against
the spec BNF as a whole any more: the seeded differential fuzzer that did
(a BNF string generator plus an independent reference acceptor) was removed
in October 2026, in commit `3b7f83f`.

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

- **Permissive valence (default backend)**: radicals and the
  expanded-octet/hypervalent species of period 3+ (e.g. `[SH6]`) are
  accepted; RDKit's default sanitizer rejects them. Second-row atoms keep
  their octet (`N(C)(C)(C)C` is rejected, as in RDKit). Opt into
  RDKit-equivalent behaviour with `backend="rdkit"`.
- **cis/trans mark on a non-single bond**: `C=N\#N` is accepted (the mark is
  dropped with a warning); RDKit rejects it as two bonds in a row.
- **Known deviations from the strict BNF.** The grammar accepts a little
  more than the spec: bare single-letter elements outside the organic
  subset (`K`, `W`, `U`) and bare aromatic `se`/`as` (E1), a multi-digit
  hcount such as `[UH12]` (E2; chemistry still applies, so `[CH12]` fails
  on valence), and zero-padded charge digit runs below magnitude 15 such as
  `[C+007]` (E3). It also rejects two spec-valid forms: ring closures may
  not span a dot (`C1.C1`, R1), and charge magnitude ≥ 15 is out of range
  (`[C+16]`, R2). Each of E1–E3, R1 and R2 is pinned by the
  `DOCUMENTED_DEVIATIONS` table in `tests/conformance/test_opensmiles_gaps.py`,
  which asserts the current behaviour above; R2 is also covered by
  `tests/chem/test_vs_pysmiles.py`.

## Historical benchmark notes

On ZINC250k the validator accepts 100% (RDKit-canonical SMILES), so
"agreement with RDKit" is saturated there by construction — which is why it
is not the headline metric. The property that matters is conformance to the
OpenSMILES grammar, which the tests listed at the top of this page check on
known cases.
