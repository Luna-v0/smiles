# Chemical Validation

This document describes the chemistry validation algorithms used by the SMILES parser to ensure molecular structures are chemically valid.

## Valence policy (decided: permissive, octet-bound for row 2)

The default chemistry backend (`ChemistryValidator(backend="pysmiles")`) is
**deliberately permissive**: it accepts radicals and the expanded-octet /
hypervalent species of period 3 and below — e.g. sulfate, `[SH6]`,
`F[P-](F)(F)(F)(F)F` — and rejects only physically impossible *over*-valence.
This is a considered choice, not an omission:

- OpenSMILES itself does not mandate a valence model; valence checking is a
  policy layered on top of the language.
- Generative-model output (the primary use case) legitimately contains
  radicals and hypervalent species that RDKit's default sanitizer rejects.

What counts as "impossible" (`chem/adapter.py`, `valence_violations`):

| Rule | Example rejected | Example accepted |
|---|---|---|
| Second-row atoms (B, C, N, O, F) cannot expand their octet: only the lowest valence of the isoelectronic element (N 3, N⁺ 4, O 2, O⁻ 1) | `N(C)(C)(C)C`, `[NH4]C`, `[OH3]C`, `C(N)=[NH2]` | `[NH4+]`, `[NH2+]=C(N)N` |
| …except the conventional neutral pentavalent-N shorthands (RDKit charge-separates them) | — | `CN(=O)=O`, `c1ccn(=O)cc1`, `C=N#N`, `N=N#N` |
| A closed-shell ion cannot bond | `[Cl-]C`, `[F-][B+3]…` | `[Na+].[Cl-]` |
| d/f-block metals get no valence limit; a neutral N/O bonded to one is a dative donor | `[NH4][Pt]…` | `[NH3][Pt]([NH3])(Cl)Cl`, `Cl[Pt-2](Cl)(Cl)Cl` |
| Main-group atoms of period 3+ use pysmiles' electron-count maximum | `O[AlH3](O)O` | `O=S(=O)(O)O` |

Until October 2026 neutral 4-valent N (`N(C)(C)(C)C`) and `[OH4]` were
accepted as "hypervalent".  Second-row atoms have no d orbitals to expand
into, and public datasets showed these strings are always a missing charge
or an extra hydrogen, so they are now rejected — with a repair suggestion.

## Repairs, body pH and salts

A chemistry failure comes with a repair when changing one atom fixes it
(`chem/repair.py`): drop hydrogens (`[AlH3]` → `[Al]`), add a charge
(`O=N([O-])C` → `O=[N+]([O-])C`), or move a bonded ion's charge onto its
partner (`[Cl-][Pt]` → `Cl[Pt-]`).  Only the edited atoms are rewritten, so
stereo and atom order are kept; `ValidationResult.repaired_smiles` holds the
result and `repairs` one explanation per atom.

When a nitrogen can be fixed either way (`=[NH2]` → `=[NH2+]` or `=[NH]`),
`physiological_ph=True` (default) picks its protonation state at pH 7.4 from a
small pKa table (`PKA_TABLE`: guanidines and amidines protonated, N-cyano
guanidines neutral, …).  With it off, the written charge is kept.

`remove_salts=True` (default) splits multi-component input into
`parent_smiles` (most heavy atoms) and `removed_fragments`.  Every component
is still validated: a broken counter-ion makes the record invalid and is named
in `failed_fragment`, so salt removal never hides a data error.

## cis/trans marks

`/` and `\` mark cis/trans direction.  OpenSMILES spells them as directional
*single* bonds, so `C=N\#N` is two bonds in a row and a syntax error there.
This validator treats the mark as an annotation instead: written in front of
another bond symbol it is accepted, ignored (a double or triple bond carries
no cis/trans direction of its own) and reported in
`ValidationResult.warnings`.

Callers who want RDKit-equivalent strictness can opt in with
`ChemistryValidator(backend="rdkit")` /
`validate_smiles_detailed(mol, backend="rdkit")`.

Note that `C11` — which the pysmiles backend used to wave through — is now
rejected *before* chemistry, by the ring-semantics tier (OpenSMILES §3.4:
an atom cannot close a ring onto itself). Ring-closure rules are language
semantics, not valence policy.

## Pipeline

`validate_smiles_detailed` runs four tiers in order (`validate_smiles` wraps
it and returns the `(is_valid, exception)` tuple). The first tier that
rejects is reported in `ValidationResult.tier`:

1. **`lex`**: the SLY lexer tokenises the string.
2. **`grammar`**: the LALR(1) parser checks OpenSMILES syntax.
3. **`ring_semantics`**: `ParserManager` enforces the ring-closure rules of
   OpenSMILES §3.4/§3.6 while it builds the molecular graph.
4. **`chemistry`**: `chem.adapter.to_descriptor` turns the graph into a
   backend-neutral atom and bond list, and `ChemistryValidator` hands it to
   the chosen backend. The backend never sees the SMILES string.
    - `pysmiles` (default): aromaticity perception and kekulization
      (`correct_aromatic_rings(strict=True)`, which fails on a region marked
      aromatic that cannot be kekulized), hydrogen filling, then the valence
      policy above (`valence_violations`).
    - `rdkit`: `Chem.SanitizeMol` on a molecule built atom by atom.

After the chemistry tier, a failure gets a repair suggestion and a
multi-component input is split into parent and counter-ions, as described
above. cis/trans warnings are collected for every input that parses.

## Legacy validator

Until October 2026 a hand-written validator applied the octet rule to
acyclic bracket atoms and Hückel's 4n+2 rule to aromatic cycles. That code
has been removed; its algorithms are kept for reference on
[Legacy Validation Algorithms](validation_algorithms.md).
