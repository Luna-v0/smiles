"""Random generator of OpenSMILES-derivable strings (§2.2 BNF).

Used by the completeness half of the differential test: every generated
string is derivable from the spec grammar *and* ring-consistent, so the
validator must accept it at the lex, grammar and ring-semantics tiers
(chemistry is a policy layer outside the language and is exempt).

The generator therefore respects the semantic constraints the spec itself
imposes on ring bonds (§3.4/§3.6) — closures pair distinct, non-adjacent
atoms, bond symbols on both ends match, ring digits precede branches — and
avoids the validator's known restrictions listed in ``test_differential.py``
(ring bonds never span a dot; charge magnitude stays below 15).
"""

import itertools
import random

BOND_SYMBOLS = ["-", "=", "#", "$", ":", "/", "\\"]
ALIPHATIC_ORGANIC = ["B", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]
AROMATIC_ORGANIC = ["b", "c", "n", "o", "p", "s"]
BRACKET_ELEMENTS = [
    "C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H",
    "Fe", "Na", "K", "Se", "Si", "Sn", "Au", "Zn", "Mg", "Co", "Pt", "As",
]
BRACKET_AROMATIC = ["b", "c", "n", "o", "p", "s", "se", "as"]
CHIRALITY = (
    ["@", "@@", "@TH1", "@TH2", "@AL1", "@AL2", "@SP1", "@SP2", "@SP3"]
    + [f"@TB{i}" for i in (1, 7, 15, 20)]
    + [f"@OH{i}" for i in (1, 12, 25, 30)]
)


def random_bracket_atom(rng: random.Random) -> str:
    """One bracket atom exercising isotope/symbol/chiral/hcount/charge/class."""
    parts = ["["]
    if rng.random() < 0.3:
        parts.append(str(rng.choice([2, 13, 35, 123, 1000])))
    parts.append(rng.choice(BRACKET_ELEMENTS + BRACKET_AROMATIC + ["*"]))
    if rng.random() < 0.25:
        parts.append(rng.choice(CHIRALITY))
    if rng.random() < 0.35:
        parts.append("H")
        if rng.random() < 0.5:
            parts.append(str(rng.randint(0, 9)))
    if rng.random() < 0.3:
        sign = rng.choice(["+", "-"])
        roll = rng.random()
        if roll < 0.4:
            parts.append(sign)
        elif roll < 0.55:
            parts.append(sign * 2)
        else:
            parts.append(sign + str(rng.randint(0, 14)))
    if rng.random() < 0.25:
        parts.append(":" + str(rng.choice([0, 5, 42, 987, 1234])))
    parts.append("]")
    return "".join(parts)


def random_atom(rng: random.Random) -> str:
    roll = rng.random()
    if roll < 0.25:
        return random_bracket_atom(rng)
    if roll < 0.30:
        return "*"
    if roll < 0.60:
        return rng.choice(ALIPHATIC_ORGANIC)
    return rng.choice(AROMATIC_ORGANIC)


def _format_ring(number: int) -> str:
    return str(number) if number < 10 else f"%{number}"


def random_fragment(
    rng: random.Random,
    ring_numbers: "itertools.count",
    depth: int = 0,
    dots_allowed: bool = True,
) -> str:
    """One connected component: a backbone with branches and ring closures.

    Args:
        rng: Seeded random source.
        ring_numbers: Global allocator so ring numbers never collide across
            nesting levels (a reused number would close the wrong ring).
        depth: Current branch nesting depth (bounded).
        dots_allowed: False when an enclosing ring is open around this
            fragment — the validator rejects a dot while any ring is open
            (restriction R1), so no dot may appear anywhere inside.
    """
    n = rng.randint(1, 8)
    atoms = [random_atom(rng) for _ in range(n)]
    bond_before = [None] * n
    ring_marks = [[] for _ in range(n)]
    branches = [[] for _ in range(n)]
    used_pairs = set()

    for i in range(1, n):
        if rng.random() < 0.25:
            bond_before[i] = rng.choice(BOND_SYMBOLS)

    # Ring closures: distinct non-adjacent atom pairs (§3.4), digits placed
    # before branches (§3.6), bond symbols matching when on both ends.
    if n >= 3:
        for _ in range(rng.randint(0, 2)):
            i = rng.randrange(0, n - 2)
            j = rng.randrange(i + 2, n)
            if (i, j) in used_pairs:
                continue
            used_pairs.add((i, j))
            number = next(ring_numbers)
            if number > 99:
                break
            tag = _format_ring(number)
            roll = rng.random()
            if roll < 0.6:
                open_mark, close_mark = tag, tag
            elif roll < 0.8:
                bond = rng.choice(BOND_SYMBOLS)
                if rng.random() < 0.5:
                    open_mark, close_mark = bond + tag, tag
                else:
                    open_mark, close_mark = tag, bond + tag
            else:
                bond = rng.choice(BOND_SYMBOLS)
                open_mark, close_mark = bond + tag, bond + tag
            ring_marks[i].append(open_mark)
            ring_marks[j].append(close_mark)

    if depth < 2:
        for i in range(n):
            # A ring spanning atom i is still open while its branches parse,
            # so dots are forbidden in that whole subtree (restriction R1).
            spans_open_ring = any(a <= i < b for a, b in used_pairs)
            child_dots = dots_allowed and not spans_open_ring
            while rng.random() < 0.12:
                lead = ""
                roll = rng.random()
                if roll < 0.25:
                    lead = rng.choice(BOND_SYMBOLS)
                elif roll < 0.32 and child_dots:
                    lead = "."
                branches[i].append(
                    "(" + lead + random_fragment(rng, ring_numbers, depth + 1, child_dots) + ")"
                )

    pieces = []
    for i in range(n):
        if bond_before[i]:
            pieces.append(bond_before[i])
        pieces.append(atoms[i])
        pieces.extend(ring_marks[i])
        pieces.extend(branches[i])
    return "".join(pieces)


def random_smiles(rng: random.Random) -> str:
    """A full SMILES: one or more dot-separated fragments.

    Ring closures never span a dot (each fragment allocates and closes its
    own rings), matching the validator's documented restriction R1.
    """
    ring_numbers = itertools.count(1)
    fragments = [random_fragment(rng, ring_numbers)]
    while rng.random() < 0.15:
        fragments.append(random_fragment(rng, ring_numbers))
    return ".".join(fragments)


def generate(seed: int, count: int) -> list[str]:
    """Deterministic batch of spec-derivable SMILES strings."""
    rng = random.Random(seed)
    return [random_smiles(rng) for _ in range(count)]
