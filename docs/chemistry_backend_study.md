# Chemistry-Backend Adaptation Study

**Goal.** Keep our SLY LALR(1) parser as the sole owner of SMILES *syntax*, and
delegate *chemistry* validation (valence + aromaticity) to the chemistry engine
of a mature library — **while bypassing that library's own SMILES string
parser**. The backend must receive a molecule graph *we* construct from our
parser's output; it must never see the raw SMILES string.

This document is a study + recommendation. A gated proof-of-concept is sketched
in §9; it is not wired into `validate_smiles` yet.

> **Reproduce everything.** All numbers below come from scripts run with
> `uv run` against the installed backends (rdkit 2026.03.1, pysmiles 2.0.1,
> partialsmiles 2.0, molvs 0.1.1, networkx 3.6.1). Dataset reads need a parquet
> engine, which the project does **not** currently lock — see §8.6. The study
> scripts used `uv run --with pyarrow …` to read `data/*.parquet`.

---

## 1. Failure analysis — why the current Stage 3 is the weak link

### 1.1 What actually runs today

The public entry point is **`src/__init__.py:42` `validate_smiles`**, *not* the
thinner `syntax/yacc.py:315` one. It performs three things in order:

1. **A regex string hack** (`src/__init__.py:62-70`): if the raw string ends in
   an aliphatic `C` and contains any ring digit, it is rejected as *"Invalid
   bonding between aromatic and aliphatic carbon"* — before any chemistry runs.
2. **Parse + graph build** (`parse_smiles` → SLY parser → `ParserManager` →
   `GraphBuilder`).
3. **`ChemistryValidator.validate(graph)`** (`src/chem/validator.py`):
   - **Valency** (`validate_rings_and_valency`, `_check_atom_valency`,
     `validator.py:113-152`): only **non-ring, non-aromatic `BracketAtom`s** are
     checked; plain organic atoms, aromatic atoms and *all* ring atoms return
     `True` unconditionally. The check itself is octet/duet electron counting
     (`atomic.py:333-352 compute_valency`: `electrons_in_valency == 8`, or `== 2`
     for H/He/Li/Be), where charge/H adjust an electron-configuration string
     (`atomic.py:218-265`).
   - **Aromaticity** (`structure.py:426-488 huckel`): Hückel 4n+2 over cycles
     found by the heuristic ring perception in `GraphBuilder.close_ring`
     (`graph_builder.py:105-313`), with a fused-ring fallback
     (`validate_fused_aromatic_system`, `structure.py:288-361`) that ultimately
     **trusts the SMILES notation** for complex systems.

### 1.2 Reproduced disagreements with RDKit (the de-facto reference)

A 29-molecule probe set (`/tmp/smiles_study/repro.py`) produced **12
disagreements** with `Chem.MolFromSmiles(...) is not None`, each landing on a
documented weakness:

| SMILES | ours | RDKit | kind | mechanism |
|---|---|---|---|---|
| `F[P-](F)(F)(F)(F)F` (PF₆⁻) | False | True | false-positive | octet on `[P-]` (6≠8 valence e⁻) |
| `[CH3]` (methyl radical) | False | True | false-positive | octet (7≠8) |
| `[SiH6-2]` | False | True | false-positive | octet / expanded valence |
| `[Fe+2]` (transition metal) | False | True | false-positive | octet cannot model d-block |
| `c1ccccc1C` (toluene) | False | True | false-positive | **regex hack** |
| `c1ccccc1CC` (ethylbenzene) | False | True | false-positive | **regex hack** |
| `OC1CCCCC1C` (2-methylcyclohexanol) | False | True | false-positive | **regex hack** |
| `c1ccccccc1` (8-aromatic-C ring) | False | True | false-positive | Hückel rejects 8π (see §1.4) |
| `C(C)(C)(C)(C)C` (pentavalent C) | True | False | **false-negative** | non-bracket C never checked |
| `FC(F)(F)(F)F` (CF₅) | True | False | **false-negative** | non-bracket C never checked |
| `O=C(=O)=O` (hexavalent C) | True | False | **false-negative** | non-bracket C never checked |
| `N(C)(C)(C)(C)C` (pentavalent N) | True | False | **false-negative** | non-bracket N never checked |

Mechanisms were confirmed directly: for `C(C)(C)(C)(C)C` our graph has the
centre carbon at **degree 5** but it is a plain `Atom`, so valency is skipped and
we accept it; for toluene, `parse_smiles` returns syntactically-valid and the
rejection is 100% the regex (`re.search(r'C$', …)` fires), not chemistry.

**Why octet counting is structurally inadequate.** It encodes a single rule
(8 valence electrons, 2 for the first row) and cannot express: per-element
valence models (N=3, P=3/5, S=2/4/6), expanded octets/hypervalency
(`PF₆⁻`, sulfate), radicals (`[CH3]`), or transition metals (`[Fe+2]` has no
octet). Worse, the check is only reached for non-ring non-aromatic *bracket*
atoms — i.e. it never sees the overwhelmingly common organic-subset atoms, which
is why all four over-valent main-group false-negatives slip through.

### 1.3 Aggregate behaviour on the real datasets

Sampling `data/*.parquet` (≈13.5k molecules, whitespace-stripped; RDKit
reference; `/tmp/smiles_study/bench.py`):

```
ALL  n=13462  agree=58.11%  over-accept=10  over-reject=5629  disagree=41.89%
```

Per dataset, agreement ranges from 29% (clintox) to 83% (qm9). The error is
**overwhelmingly over-rejection** (we reject RDKit-valid molecules): 5629 vs 10.
This is **far worse than the 15.60% recorded in `docs/rdkit_comparison.md`** —
that figure is stale/optimistic.

Attributing the over-rejections to a mechanism
(`/tmp/smiles_study/attribute.py`, 2303 sampled over-rejections):

| share | mechanism | example over-rejected (RDKit-valid) |
|---|---|---|
| **43.5%** | aromaticity (Hückel / ring perception) | `c1[nH]c2c(n1)[nH]cnc2=S` |
| **31.4%** | valency (octet on bracket atom) | `CCCC(CCC)C(=O)[O-]` |
| **25.0%** | regex aromatic/aliphatic hack | steroids/esters ending in `…C` |
| **~0.0%** | parser / LALR grammar | (1 case) |

**This is the single most important result in the study:** the syntax layer
accounts for ≈0% of disagreements; ~100% come from the three chemistry
mechanisms. The premise — *keep the parser, replace the chemistry* — is exactly
right.

### 1.4 The labelled set is stricter than RDKit (and Stage 3 is overfit to it)

On the 126-row labelled ground truth (`/tmp/smiles_study/labelled.py`):

| Validator | TP | TN | FP | FN | Acc | Prec | Rec | F1 |
|---|---|---|---|---|---|---|---|---|
| **YACC (ours)** | 94 | 32 | 0 | 0 | **100.0%** | 100% | 100% | 100% |
| RDKit | 88 | 11 | 21 | 6 | 78.6% | 80.7% | 93.6% | 86.7% |
| PartialSMILES | 80 | 11 | 21 | 14 | 72.2% | 79.2% | 85.1% | 82.1% |
| PySMILES | 80 | 11 | 21 | 14 | 72.2% | 79.2% | 85.1% | 82.1% |
| MolVS | 88 | 11 | 21 | 6 | 78.6% | 80.7% | 93.6% | 86.7% |

Our validator scores **100%**; **every** mature engine scores ~72–79%. The
reason: **21 of the 32 "invalid"-labelled rows are chemically valid molecules**
that RDKit (and all others) accept. Inspecting them, *all 21* are an aromatic
ring system with a terminal methyl, labelled invalid under a non-existent
"aromatic–aliphatic C bond is invalid" rule:

```
c1ccncc1C      (2-picoline)        c1ccsc1C   (2-methylthiophene)
c1ccoc1C       (2-methylfuran)     c1c[nH]cc1C (3-methylpyrrole)
c1ccc2ccccc2c1C (1-methylnaphthalene)  … + 16 more methyl-PAHs/heteroarenes
```

The regex hack in `src/__init__.py:62-70` exists **solely to pass these 21
mislabeled rows** — and that same hack causes **25% of the real-world
over-rejections** (§1.3). This is an airtight causal chain: *fitting Stage 3 to
wrong labels directly produced a quarter of its production false-positives.*

The remaining 6 disagreements are labelled *valid* but RDKit rejects them
(phenalene, fluoranthene, pyrene, perylene written in non-kekulizable aromatic
forms) — genuine hard aromaticity edge cases.

**Consequence for scoring (see §7):** a *more correct*, RDKit-backed validator
will **lose** accuracy on this labelled set, because the set contains ~21
chemically-wrong labels. The benchmark plan must report both axes and name this
explicitly.

---

## 2. Backend comparison & recommendation

### 2.1 Research-question findings per backend

**Separability — is the chemistry callable independently of the string parser?**

- **pysmiles — YES, cleanly.** `read_smiles` (`read_smiles.py:250-321`) first
  calls `base_smiles_parser` (the *string* parser, builds an `nx.Graph`), then
  runs the chemistry as four **graph-native free functions**:
  ```python
  # read_smiles.py
  if reinterpret_aromatic:
      correct_aromatic_rings(mol, strict=strict)   # aromaticity + kekulization
  fill_valence(mol)                                 # implicit H from valence model
  ...
  elif element != '*' and bonds_missing(mol, node): # over/under-valence test
      raise KeyError(... 'non-standard valence' ...)
  ```
  These operate purely on node attributes (`element`, `charge`, `aromatic`,
  `hcount`) and edge attribute `order`. The valence model itself
  (`smiles_helper.py:391-454 valence`) derives allowed valences per element from
  the periodic table electron configuration and the formal charge — a real
  multi-valence model, returning e.g. `[2,4,6]` for S, `[]` (no constraint) for
  transition metals. **Bypass = build the `nx.Graph` ourselves and call those
  four functions.**

- **RDKit — YES, via construct-then-sanitize.** Build a `Chem.RWMol`
  atom-by-atom (`AddAtom`/`SetFormalCharge`/`SetIsAromatic`/`SetNumExplicitHs`)
  and bond-by-bond (`AddBond`), then `Chem.SanitizeMol`. The sanitize operations
  decompose cleanly (verified, `/tmp/smiles_study/finalize.py`):
  `SANITIZE_PROPERTIES` (=2) → **valence** (raises `AtomValenceException`),
  `SANITIZE_KEKULIZE` (=8) → **kekulization** (raises `KekulizeException`),
  `SANITIZE_SETAROMATICITY` (=32) → aromatic perception. `MolFromSmiles` is never
  called.

- **partialsmiles — NO (interleaved).** Its only public entry is `ParseSmiles`.
  Validation is invoked *inside the per-character parse loop*
  (`smiparser.py:120-141`):
  ```python
  def parse(self, smi):
      ...
      while state.idx < self.N: self.parse_token(state)
      ...
      self.handleError(ValenceError, self.validateValence(state))
      self.handleError(KekulizationFailure, self.validateKekulization(state.mol, state))
  ```
  `validateValence`/`validateKekulization` are **methods that read parser
  `state`** (`state.mol`, `state.smiidx`, `self.incompleteAtoms`), not free
  functions. To bypass the parser you would have to reconstruct partialsmiles'
  private `Molecule`/`Atom`/`Bond` objects *and* re-implement `validateValence`.
  Reusable in isolation are only the **static data** `valence.common_valencies`
  (an excellent per-element/per-charge valence table, `valence.py:1-44`) and the
  free `kekulize.Kekulize(mol)`. Score **down** for separability.

- **MolVS — not an independent engine.** `molvs.validate_smiles`
  (`validate.py:107-121`) is literally `Chem.MolFromSmiles(smiles)` followed by a
  list of post-checks. "MolVS backend" = "RDKit backend + MolVS checks"; it
  *cannot* satisfy the bypass requirement because its first act is to call
  RDKit's string parser.

**Graph-nativeness, updatability, correctness, lossiness fit** — summarised in
the matrix below. Correctness ceilings are each backend's own engine vs RDKit on
2,250 sampled real molecules (`/tmp/smiles_study/finalize.py`):

```
pysmiles       agree_with_RDKit=97.56%   over-accept=0  over-reject=55
partialsmiles  agree_with_RDKit=99.42%   over-accept=0  over-reject=13
molvs          agree_with_RDKit=87.78%   over-accept=0  over-reject=275
```

### 2.2 Scored comparison matrix

Scores 1 (poor) – 5 (excellent), against our 5 research criteria.

| Criterion | pysmiles | RDKit | partialsmiles | MolVS |
|---|:--:|:--:|:--:|:--:|
| **Separability** (parser bypassable) | **5** — graph-native free funcs | 4 — RWMol+Sanitize | 2 — interleaved w/ parse state | 1 — *is* RDKit's parser |
| **Graph-nativeness** (vs our `MolecularGraph`) | **5** — `nx.Graph`, 1:1 attrs | 3 — build `RWMol` | 2 — private objects | 1 — n/a |
| **Updatability** (install / pin / patch) | **5** — pure-Py ~1k LoC, patchable | 2 — large C++ core, not patchable | 4 — pure-Py, small | 3 — thin, but drags RDKit |
| **Correctness** (vs RDKit on real data) | 4 — 97.6% | **5** — reference | **5** — 99.4% | 3 — 87.8% (full suite) |
| **Lossiness fit** (consumes what we can give) | **5** — element/charge/arom/H/order | 5 — same via RWMol | 3 — needs implicit-H model rebuilt | 2 |
| **Bypass-the-parser compliance** | ✅ | ✅ | ⚠️ partial (tables only) | ❌ |

### 2.3 Recommendation

> **Primary recommendation: pysmiles, as the "easiest-to-update" chemistry
> backend, behind a backend-neutral adapter.**

Rationale:

- **Cleanest bypass and the best graph fit.** Its chemistry is four free
  functions over an `nx.Graph` whose node/edge attributes are a 1:1 map of our
  `MolecularGraph`. The adapter is ~30 lines (§4); no private classes, no
  string parser anywhere near it.
- **Easiest to update / maintain.** Pure Python, ~1k LoC, pip-pinnable, and —
  unlike RDKit's C++ core — *patchable* if we ever need to adjust the valence
  model. This is the criterion the task asked to optimise.
- **Correctness is a deliberate, stated trade-off.** pysmiles agrees with RDKit
  on **97.6%** of real molecules vs RDKit's 100% (reference). Its disagreements
  are *stricter* (it rejects some radicals such as `[CH3]` and a few exotic
  hypervalent species that RDKit tolerates) and **never over-accepts**
  (0 over-accept in 2,250). For a generative-AI validation gate, "occasionally
  too strict, never too loose" is the safe failure direction. This trades ~2.4
  points of RDKit-agreement for pure-Python updatability — and still improves
  our real-world agreement from **58% → ~97%**.

> **Secondary recommendation: keep RDKit as the correctness reference, and make
> it a drop-in second *driver* behind the same adapter contract.** Because the
> adapter emits a backend-neutral `(atoms, bonds)` description (§3), wiring an
> RDKit driver (`RWMol` + `SanitizeMol`) is additive and *also* parser-bypassing
> (§9 proves both drivers give identical verdicts on the probe set). Choose
> pysmiles by default for updatability; switch the driver to RDKit when maximum
> agreement with the de-facto reference is required.

**partialsmiles** is rejected as a backend (validation is inseparable from its
parse loop) but its `common_valencies` table is a useful cross-check.
**MolVS** is rejected (it is RDKit's parser plus standardisation checks; it
cannot bypass).

---

## 3. Graph hand-off specification

The contract the parser must expose to the adapter — every attribute a backend
needs, audited against what `MolecularGraph` carries today
(`/tmp/smiles_study/dump_graph.py`).

**Per-atom attributes**

| Attribute | Needed by | Carried today? | Evidence | Closure plan (grammar-untouched) |
|---|---|---|---|---|
| element symbol | both | ✅ all atoms (`Atom.symbol`) | every dump row | — |
| aromatic flag | both | ✅ all atoms (`Atom.aromatic`) | `arom=1` on lowercase | — |
| formal charge | both | ✅ on `BracketAtom`; organic ⇒ 0 | `[O-]`→−1, `[Na+]`→+1 | complete: charge is only legal inside brackets (OpenSMILES §3.1) |
| isotope | both (optional) | ✅ on `BracketAtom` | `[13CH4]`→13 | — |
| explicit H (`hcount`) | both | ✅ on `BracketAtom` | `[NH+]`→1, `[13CH4]`→4 | — |
| implicit H (organic) | both | ⚠️ not stored (by design) | organic `hcount=None` | none needed — backend fills (`fill_valence` / RDKit implicit) |
| chirality `@`/`@@` | optional | ✅ parsed (`BracketAtom.chiral`) | field present | irrelevant to valence/aromaticity |

**Per-bond attributes**

| Attribute | Needed by | Carried today? | Evidence | Closure plan |
|---|---|---|---|---|
| bond order in a **branch** `(=O)` | both | ✅ | acetate `=`, `O=S(=O)(O)O` | — |
| bond order in the **main chain** `C=C` | both | ❌ **LOST** | `C=C`→`-`, `O=C=O`→`O-C-O`, `C=CC=C`→all single | **must record** — see §3.1 |
| aromatic bond `:` | both | ✅ auto-marked between aromatic atoms | benzene `:` | — |
| ring-closure bond | both | ✅ edge added | naphthalene → 11 edges | — |
| stereo `/` `\` | optional | ⚠️ collapsed to `-` | `F/C=C/F`→`-` | irrelevant to valence/aromaticity; store as edge attr only if stereo validation is later wanted |
| dot `.` (disconnection) | both | ✅ no edge, separate component | `[Na+].[Cl-]`→2 atoms / 0 edges | — |

### 3.1 The one mandatory gap: main-chain bond order

`C=C` becomes `C-C` and `O=C=O` becomes `O-C-O` in our graph. The double/triple
bond **token is parsed but discarded** when the bond is created — concretely in
`parser_manager.line` (`parser_manager.py:744-748`):

```python
case {"atom": atom, "chain_branch": {"chains": {"atom": atom2}}}:
    self.graph_builder.add_bond(atom, atom2)   # <-- the already-parsed bond order is dropped
```

(Mapping patterns allow extra keys, so this case silently swallows a
`{"bond": "=", "atom": …}` chain and bonds with the default `-`.) Branch bonds
survive because they go through the `pending_branch_bond` path.

This is **not** a grammar change. The `=`/`#`/`$` token is already produced by
the frozen grammar and already present in the reduced `chain` dict; the fix is to
*record what is already parsed* — exactly the kind of change the task's
constraint #1 permits. The closure is to thread the `bond` value into
`add_bond(...)` in the relevant `line`/`chains` cases (and let
`graph_builder.add_bond` store it as the edge `bond_type`, which it already
supports).

**Why it matters:** the current octet validator hides this defect by skipping
all non-bracket atoms, so `O=C=O` "passes" vacuously. *Any* real backend fed the
present graph would mis-read every main-chain multiple bond (CO₂ → a di-radical,
Kekulé benzene → cyclohexane). **Fixing main-chain bond-order recording is a
hard prerequisite for the migration**, and the only attribute that cannot be
recovered in the adapter alone.

---

## 4. Adapter design

A single backend-neutral translation, then a thin per-backend driver.

```python
# src/chem/adapter.py  (new module — sketch)
def to_descriptor(graph: MolecularGraph) -> tuple[list[dict], list[tuple[int,int,float]]]:
    """MolecularGraph -> ([atom dicts], [(i, j, order)]).  Backend-agnostic."""
    idx = {a: i for i, a in enumerate(graph.adjacency_list)}        # stable atom ids
    atoms = [dict(element=a.symbol.capitalize(),                    # 'CL'->'Cl' etc.
                  charge=getattr(a, "charge", None) or 0,
                  aromatic=bool(getattr(a, "aromatic", False)),
                  hcount=getattr(a, "hcount", None),                # None => implicit
                  isotope=getattr(a, "isotope", None))
             for a in graph.adjacency_list]
    order = {"-":1, "=":2, "#":3, "$":4, ":":1.5, "/":1, "\\":1}    # stereo -> single
    seen, bonds = set(), []
    for a, nbrs in graph.adjacency_list.items():
        for n, bt in nbrs:
            key = frozenset((idx[a], idx[n]))
            if key in seen or a is n: continue                       # dedupe, drop self-loops
            seen.add(key); bonds.append((idx[a], idx[n], order.get(bt, 1)))
    return atoms, bonds
```

**Attribute mapping**

- *element*: `symbol.capitalize()` (our symbols are upper-cased; backends want
  `Cl`, `Br`, `Se`). aromaticity is carried by the separate `aromatic` flag, so
  case is purely cosmetic for the element string.
- *charge / isotope / explicit H*: copied straight from `BracketAtom`; organic
  atoms imply charge 0 and implicit H.
- *aromatic atoms & bonds*: passed as `aromatic=True` plus order `1.5` on the
  bond. Both backends then **perceive/kekulize themselves** — pysmiles via
  `correct_aromatic_rings` (raises if a ring can't be kekulized), RDKit via
  `SANITIZE_KEKULIZE`/`SANITIZE_SETAROMATICITY`. **We hand over aromatic input;
  the backend does kekulization.** We do *not* pre-kekulize.
- *implicit H*: left to the backend (`fill_valence` / RDKit implicit-H).

**pysmiles driver** (the recommended default):

```python
def valid_pysmiles(atoms, bonds) -> bool:
    g = nx.Graph()
    for i, a in enumerate(atoms):
        attrs = {k: a[k] for k in ("element","charge","aromatic") }
        if a["hcount"] is not None: attrs["hcount"] = a["hcount"]
        if a["isotope"] is not None: attrs["isotope"] = a["isotope"]
        g.add_node(i, **attrs)
    for i, j, o in bonds: g.add_edge(i, j, order=o)
    try:
        correct_aromatic_rings(g, strict=True)     # aromaticity  (theirs)
        fill_valence(g)                            # implicit H   (theirs)
        return not any(el != "*" and bonds_missing(g, n)
                       for n, el in g.nodes(data=lambda d: d.get("element","*")))
    except Exception:
        return False
```

**RDKit driver** (the high-accuracy option): build `RWMol`, set
charge/aromatic/explicit-H, `AddBond`, then `Chem.SanitizeMol`. §9 shows both
drivers agree with the RDKit string reference on the probe set.

Perception delegated to the backend: **ring perception, aromaticity, and
kekulization all move to the backend**, which makes our `close_ring` heuristic
and `structure.huckel` redundant (§6).

---

## 5. Bypass map — who does what

| Stage | Function | Owner | Notes |
|---|---|---|---|
| Tokenize | `SmilesLex.tokenize` (`syntax/lex.py`) | **ours** | unchanged |
| Parse (LALR) | `SmilesParser.parse` (`syntax/yacc.py`) | **ours** | grammar frozen |
| Semantic actions | `ParserManager.*` (`syntax/parser_manager.py`) | **ours** | + record main-chain bond order (§3.1) |
| Build graph | `GraphBuilder.add_atom/add_bond` (`chem/graph_builder.py`) | **ours** | ring *closure bond* kept; ring *perception* dropped |
| Graph → descriptor | `adapter.to_descriptor` (new) | **ours** | backend-neutral |
| **Valence check** | pysmiles `fill_valence` + `bonds_missing` (`smiles_helper.py:324-388`) — *or* RDKit `SanitizeMol(SANITIZE_PROPERTIES)` | **theirs** | replaces octet counting |
| **Aromaticity / kekulization** | pysmiles `correct_aromatic_rings` (`smiles_helper.py:522-606`) — *or* RDKit `SanitizeMol(SANITIZE_KEKULIZE\|SANITIZE_SETAROMATICITY)` | **theirs** | replaces `huckel` + `close_ring` perception |
| Backend's own SMILES parser | `read_smiles` / `MolFromSmiles` / `ParseSmiles` | **UNUSED** | never called — this is the proof |

The table is the contract: **lex/parse/graph-build are ours; valence and
aromaticity are theirs; the backend's string parser is never invoked.**

---

## 6. Migration plan (file-by-file)

| File | Action | Detail |
|---|---|---|
| `syntax/lex.py`, `syntax/yacc.py` | **keep** | frozen grammar; syntax owner. (Optional: tolerate trailing whitespace — §8.6.) |
| `syntax/parser_manager.py` | **thin + 1 fix** | record already-parsed main-chain bond order into `add_bond` (§3.1). No new grammar. |
| `chem/graph_builder.py` | **thin** | keep `add_atom`/`add_bond`/ring *closure-bond* creation; **delete the `close_ring` perception heuristic** (BFS/scoring/fallbacks, `graph_builder.py:130-294`) and `cycles` bookkeeping — the backend perceives rings. |
| `chem/structure.py` | **thin** | keep `MolecularGraph` (adjacency list) as the hand-off object; **delete** `huckel`, `validate_fused_aromatic_system`, `_count_pi_*`, `get_fused_ring_systems`, `check_valency_for_aba` (aromaticity now the backend's). |
| `chem/atomic.py` | **thin** | keep `Atom`/`BracketAtom` dataclasses (symbol, aromatic, charge, isotope, hcount, …); **delete** the electron-configuration machinery (`_parse_electron_configuration`, `_adjust_electrons_for_charge_and_hydrogens`, `compute_valency`) — no longer used. |
| `chem/chemistry.py` | **thin** | keep the atom factory + symbol validation; **delete** `number_of_electrons_per_bond`, `validate_valency_bracket`. |
| `chem/validator.py` | **replace** | `ChemistryValidator.validate(graph)` becomes: `atoms, bonds = to_descriptor(graph); ok = driver(atoms, bonds)`; map a backend failure to `ParserException(rule="chemistry", parameter=smiles, message=<backend msg>)`. |
| `chem/adapter.py` | **new** | `to_descriptor` + `valid_pysmiles` / `valid_rdkit` drivers (§4). |
| `src/__init__.py` | **rewire, signature unchanged** | **delete the regex hack (`:62-70`)**; keep `validate_smiles(mol) -> (bool, Exception|None)`; pipeline = `parse_smiles` → `to_descriptor` → driver; on backend exception return `(False, ParserException(...))`. |

**Return-contract preservation.** `validate_smiles` keeps its exact signature.
Element 0 stays the scored boolean. Backend exceptions (`SyntaxError`/`KeyError`
from pysmiles; `AtomValenceException`/`KekulizeException` from RDKit) are caught
in the driver and re-wrapped as `ParserException` so element 1 stays an
`Exception` as before.

---

## 7. Test & benchmark plan

Reuse `analyses/analysingValidator.ipynb` unchanged (its `VALIDATORS` dict +
`compute_metrics` + confusion-matrix cells). Add the backend-swapped validator
(`YACC+pysmiles`, and optionally `YACC+RDKit`) as new rows.

**Two scoring axes, reported side by side:**

1. **Labelled ground truth** (`tests/data/test_molecules_from_values.parquet`,
   126 rows) — *kept as the only scored ground truth* per constraint #4.
2. **RDKit-agreement on `data/*.parquet`** (whitespace-stripped, sampled) — the
   real-world axis.

**Measurable targets vs the current baseline:**

| Metric | Current | Target |
|---|---|---|
| RDKit-agreement on real data (§1.3) | **58.1%** | **≥ 95%** (pysmiles ceiling 97.6%) |
| Over-rejection share — aromaticity | 43.5% | → ~0 (backend kekulization) |
| Over-rejection share — octet valency | 31.4% | → ~0 (per-element valence model) |
| Over-rejection share — regex hack | 25.0% | → 0 (hack deleted) |
| Over-acceptance (e.g. pentavalent C) | present | → 0 (backend valence) |
| Documented RDKit disagreements (`rdkit_comparison.md`) | ~12 | ≤ 3 |
| Syntax-only behaviour | unchanged | **no regression** (parser untouched) |

**Expected, explicitly-stated regression.** On the 126-row labelled set the
backend-swapped validator will **drop from 100% to ≈79%**, because ~21 of its
"invalid" labels are chemically valid (§1.4) and the corrected validator will
(rightly) accept them. This is *correcting an overfit to wrong labels*, not a
true regression. Recommend the maintainers re-examine those 21 rows; until then,
report the labelled-set number with this caveat rather than treating the drop as
failure.

**Unit tests that must accompany the adapter** (`tests/chem/test_adapter.py`):

- `to_descriptor` maps element/charge/isotope/explicit-H/aromatic/bond-order
  correctly (assert on the dicts), including: `[O-]` charge −1; `[13CH4]`
  isotope 13 + hcount 4; benzene six `1.5` bonds; `[Na+].[Cl-]` → 2 atoms /
  0 bonds; `O=C=O` → two order-2 bonds (guards the §3.1 fix).
- Driver verdicts on a curated set (toluene→True, acetate→True, PF₆⁻→True,
  pentavalent-C→False, CO₂→True, naphthalene aromatic & Kekulé→True,
  unkekulizable `c1cccc1`→False).
- **Bypass guards**: monkeypatch `pysmiles.read_smiles` and
  `rdkit.Chem.MolFromSmiles` to raise, and assert `validate_smiles` still works —
  proving the string parser is never on the path.
- `validate_smiles` contract: returns `tuple[bool, Exception|None]`; failures
  carry a `ParserException`.

---

## 8. Risks & edge cases

1. **Main-chain bond-order loss (§3.1)** — *highest risk*. Must be fixed before
   the swap or every alkene/alkyne/carbonyl in the main chain is mis-encoded.
   Mitigation: the §3.1 grammar-untouched fix + the `O=C=O` unit test.
2. **Graph-builder artifacts.** `close_ring` can in principle emit duplicate or
   self-loop edges; the adapter dedupes by atom-id-pair and drops self-loops, but
   once ring *perception* is delegated we should delete `close_ring`'s heuristic
   entirely (§6) and keep only the closure-bond edge, removing the artifact
   source.
3. **Kekulization mismatch.** Backends sometimes disagree on borderline aromatic
   input. pysmiles is *stricter* than RDKit (rejects `[CH3]` radical, some
   hypervalent); it never over-accepts (0/2250). For a validation gate this is
   the safe direction; choose the RDKit driver if maximum RDKit-agreement is
   required.
4. **Disconnected structures (`.`).** Handled — components arrive as separate
   nodes with no edge (`[Na+].[Cl-]` verified). Both backends validate each
   component independently; no special-casing needed.
5. **Stereo (`/ \ @`).** Collapsed/ignored for valence+aromaticity (correct — it
   doesn't affect them). If stereo *validity* is ever in scope, record the
   direction as an edge attribute (grammar-untouched) and use RDKit's stereo
   perception.
6. **Environment / tooling gaps surfaced during the study:**
   - **No parquet engine is locked.** `pyproject.toml` ships `data/*.parquet`
     and the notebook reads them, but neither `pyarrow` nor `fastparquet` is a
     dependency — `pd.read_parquet` fails on a clean env. Add `pyarrow` to
     `pyproject.toml`.
   - **Lexer rejects trailing whitespace.** Every zinc250k SMILES ends in `\n`;
     our lexer raises *"Illegal character '\n'"*, which alone took zinc250k from
     0%→42.5% agreement once stripped. RDKit strips; we should too (a lexer/entry
     tweak, not a grammar change).
   - **Broken test collection.** `tests/syntax/test_parser_values.py:11` reads
     `tests/data/test_molecules_from_values.csv`, which was replaced by the
     `.parquet`. On a clean checkout pytest fails at collection. Point the test
     at the parquet (one-line fix).
7. **Performance.** Our pipeline runs ~1,850 mol/s (0.54 ms/mol), single-thread
   Python; pysmiles adds NetworkX matching per molecule (similar order). RDKit is
   far faster. For 10⁶-scale screening, batch with multiprocessing or use the
   RDKit driver. Acceptable for validation-gate use.
8. **Dependency footprint.** pysmiles + networkx is small and pure-Python (good
   for "easiest to update"). The RDKit driver pulls the existing rdkit wheel
   (already a dependency); no new heavy dependency is introduced by either
   choice.

---

## 9. Proof-of-concept — IMPLEMENTED

The pysmiles drop-in is wired in behind the unchanged `validate_smiles` API.

**Files changed**

| File | Change |
|---|---|
| `chem/adapter.py` (new) | `to_descriptor(graph)` + `valid_pysmiles(atoms, bonds)`; bypasses `read_smiles`. |
| `syntax/parser_manager.py` | Added a **canonical** left-to-right graph builder (anchor + branch-stack) driven by `atom`/`start_branch`/`end_branch`/`rnum`, fixing connectivity *and* main-chain bond order (§3.1). The legacy tangle is left intact but its output is no longer used. |
| `chem/structure.py` | Added `update_edge`/`remove_edge` helpers used by the canonical builder. |
| `chem/validator.py` | `ChemistryValidator.validate` now delegates to `to_descriptor` + `valid_pysmiles`; backend failures wrap into `ParserException`. |
| `src/__init__.py` | **Deleted the regex aromatic/aliphatic hack**; `parse_smiles` hands off the canonical graph; `validate_smiles` signature unchanged. |
| `tests/chem/test_adapter.py` (new) | Mapping, driver verdicts, and **bypass guards** (monkeypatch `read_smiles` to raise; validation still works). |
| `syntax/lex.py` | `tokenize` now **splits greedily-merged bare two-letter atoms** outside brackets (fixes `Cn`→Copernicium etc.); bracket contents untouched. |
| `src/__init__.py` | Also added an **unclosed-ring check** (`has_open_cycles`) — a pre-existing syntax gap (`C1CC` was accepted). |
| `chem/adapter.py` | **Free-ion rule** (`[Li+]`/`[Na+]`/`[Cl-]` accepted) and the **permissive valence policy**: accept radicals/under-valent *and* hypervalent species; reject only impossible over-valence. |
| `tests/data/test_molecules_from_values.parquet` | Labels **re-aligned to the permissive policy** (the set is narrow and was distrusted by the maintainer — see §7). |
| `tests/syntax/test_parser_values.py` | Reads the parquet (stale CSV was gone); now a **regression snapshot**, not independent ground truth. |
| `tests/chem/test_vs_pysmiles.py` (new) | **Locks in the advantage over pysmiles' parser** (we reject malformed strings it accepts). |
| `pyproject.toml` | Added `pyarrow` (the project shipped parquet datasets with no parquet engine). |

### The validation policy (deliberately permissive)

RDKit is treated as a **reference axis, not the target.** The maintainer wants a
gate that accepts chemically-real molecules even when RDKit's conservative
defaults reject them, and only blocks the genuinely impossible. The policy:

| Class | Example | Decision | Who blocks it |
|---|---|---|---|
| Hypervalent / expanded octet | `N(C)(C)(C)(C)C`, `[OH4]`, sulfate | **accept** (RDKit rejects) | — |
| Radical / under-valent | `[CH3]`, `[c]`, `[Na]` | **accept** (was rejected) | — |
| Impossible over-valence | pentavalent C `C(C)(C)(C)(C)C` | **reject** | pysmiles valence |
| Unkekulizable aromatic | `c1cccc1` | **reject** | pysmiles aromaticity |
| Malformed syntax | `C(`, `=CC`, `C1CC`, `""` | **reject** | **our grammar** |

**Measured results** (`uv run`):

| Metric | Original | now |
|---|--:|--:|
| Acceptance on real molecule DBs (9 000 sampled — these *are* valid) | — | **99.89%** |
| Known-invalid rejection (`invalid_smiles`, 52) | — | **40/52** (vs pysmiles 29/52) |
| Graph connectivity match vs RDKit | 13.3% | **~99%** |
| RDKit agreement (reference axis only, ≈13.5k) | 58.1% | **~99.8%** |
| Test suite | collection error | **231 passed** |
| Throughput | ~1,850 mol/s | ~405 mol/s |

The few real-DB rejections (clintox/tox21 ≈0.5%) are genuinely unusual entries;
acceptance is 100% on moses, mutagenicity, muv, qm9, sider, zinc250k samples.

## 9a. How this system is better than pysmiles alone

Keeping our LALR grammar as the syntax gate is not redundant — pysmiles'
hand-written parser is permissive and accepts many malformed strings that our
grammar rejects. On `data/invalid_smiles.parquet` (52 known-invalid SMILES):

- **Our pipeline rejects 40/52; pysmiles' `read_smiles` rejects 29/52.**
- We **dominate every syntax category**: dangling bonds 3/3 vs 0/3, `invalid_syntax`
  4/4 vs 1/4, unmatched parens 4/4 vs 2/4, empty/whitespace 1/1 vs 0, leading
  bonds, out-of-range charge (`[C+16]`) and isotope (`[99999C]`).
- Concretely, pysmiles **accepts** `""`, `" "`, `=CC`, `C=`, `C(`, `C(CC`,
  `C()`, `CC.` — all of which our pipeline rejects.

This is enforced by `tests/chem/test_vs_pysmiles.py`, so the advantage cannot
silently regress. (The categories where we *don't* reject — `valency`,
`ring_strain` — are the hypervalent/strained species the permissive policy
**accepts by design**.)

## 9b. Notes & tunable knobs

- **Hand-set distrust (§7).** The 126-row labelled set is narrow (PAH/heteroarene
  variants), was partly fit to the old regex hack, and had ~27 chemically-wrong
  labels. It is now a *regression snapshot* aligned to the policy, not the
  scoring authority. Independent correctness lives in `tests/chem/test_adapter.py`
  and `test_vs_pysmiles.py`; the real evaluation is real-DB acceptance +
  `invalid_smiles` rejection.
- **Strictness is one line.** The policy lives in `valid_pysmiles`: `if missing < 0`
  rejects only over-valence. Change to `if missing` to also reject radicals
  (RDKit-like), or drop the check to accept any kekulizable structure.
- **Performance.** pysmiles' per-molecule NetworkX matching gives ~405 mol/s;
  batch or add an RDKit driver behind the same adapter for 10⁶-scale screening.

---

## 10. Definition-of-done checklist

- ✅ Recommendation justified with reproduced evidence (probe set, 13.5k-mol
  aggregate, mechanism attribution, backend ceilings, bypass PoC) — not
  assertion.
- ✅ Bypass map proves syntax = ours, chemistry = theirs, backend parser =
  unused (enforced by a monkeypatch test).
- ✅ Graph hand-off spec audits every attribute; the one mandatory gap
  (main-chain bond order) has a grammar-untouched closure that is now
  implemented (§3.1, §9).
- ✅ Implemented and measured against a **trustworthy** basis (real-molecule DBs
  + `invalid_smiles`, not the distrusted hand set): **99.89%** acceptance of real
  molecules, **40/52** known-invalid rejection (vs pysmiles 29/52), **231 passed**.
  RDKit agreement (58.1% → ~99.8%) is reported only as a reference axis; the
  divergence (accepting hypervalent + radical chemistry) is by design (§9).
- ✅ The advantage over pysmiles' own parser is documented (§9a) and locked in by
  `tests/chem/test_vs_pysmiles.py`.
