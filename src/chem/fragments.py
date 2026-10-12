"""Split a multi-component SMILES into its parent molecule and counter-ions.

Salts and mixtures are written as dot-separated components (``CC(=O)[O-].[Na+]``).
For modelling only the *parent* (the largest component) usually matters, so
the validator reports it separately — but every component is still
validated, so a malformed counter-ion is reported rather than silently
dropped.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

from chem.adapter import AtomDescriptor, BondDescriptor


@dataclass
class Fragment:
    """One connected component of a molecule.

    Attributes:
        smiles: The component's text, sliced from the original string.
        atoms: Atom descriptors of the component (re-indexed from 0).
        bonds: Bond descriptors of the component (re-indexed).
    """

    smiles: str
    atoms: List[AtomDescriptor]
    bonds: List[BondDescriptor]

    @property
    def heavy_atoms(self) -> int:
        return sum(1 for a in self.atoms if a["element"] not in ("H", "*"))

    @property
    def carbons(self) -> int:
        return sum(1 for a in self.atoms if a["element"] == "C")


def _top_level_dots(smiles: str) -> List[int]:
    """Positions of ``.`` outside brackets and branches."""
    dots, depth, in_bracket = [], 0, False
    for i, ch in enumerate(smiles):
        if ch == "[":
            in_bracket = True
        elif ch == "]":
            in_bracket = False
        elif not in_bracket and ch == "(":
            depth += 1
        elif not in_bracket and ch == ")":
            depth -= 1
        elif not in_bracket and depth == 0 and ch == ".":
            dots.append(i)
    return dots


def _components(n_atoms: int, bonds: List[BondDescriptor]) -> List[int]:
    """Component label of every atom (union-find)."""
    parent = list(range(n_atoms))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j, _ in bonds:
        parent[find(i)] = find(j)
    return [find(i) for i in range(n_atoms)]


def split_fragments(
    smiles: str, atoms: List[AtomDescriptor], bonds: List[BondDescriptor]
) -> Optional[List[Fragment]]:
    """Split a molecule into its dot-separated components.

    Returns:
        The fragments in source order, or ``None`` when the text pieces do not
        line up with the bonded components (a ring bond written across a
        dot, ``C1.C1``, or a dot inside a branch).
    """
    from chem.repair import atom_spans

    cuts = _top_level_dots(smiles)
    bounds: List[Tuple[int, int]] = []
    start = 0
    for cut in cuts + [len(smiles)]:
        bounds.append((start, cut))
        start = cut + 1

    spans = atom_spans(smiles)
    piece_of = []
    for atom_start, _ in spans:
        piece_of.append(next(k for k, (lo, hi) in enumerate(bounds) if lo <= atom_start < hi))

    labels = _components(len(atoms), bonds)
    pieces = {}
    for atom, (piece, label) in enumerate(zip(piece_of, labels)):
        if pieces.setdefault(piece, label) != label:
            return None  # one text piece holds two components
    if len(set(pieces.values())) != len(pieces) or len(pieces) != len(bounds):
        return None  # one component spans two pieces, or an empty piece

    fragments = []
    for piece, (lo, hi) in enumerate(bounds):
        members = [a for a in range(len(atoms)) if piece_of[a] == piece]
        remap = {old: new for new, old in enumerate(members)}
        fragments.append(
            Fragment(
                smiles=smiles[lo:hi],
                atoms=[atoms[a] for a in members],
                bonds=[(remap[i], remap[j], o) for i, j, o in bonds if i in remap],
            )
        )
    return fragments


def parent_index(fragments: List[Fragment]) -> int:
    """Index of the parent: most heavy atoms, then most carbons, then first."""
    return max(
        range(len(fragments)),
        key=lambda k: (fragments[k].heavy_atoms, fragments[k].carbons, -k),
    )
