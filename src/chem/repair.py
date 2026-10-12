"""Repair suggestions for SMILES that are chemically invalid as written.

Most invalid SMILES in real datasets describe a perfectly real molecule with
one atom written wrong: a missing charge (``O=N([O-])C`` for a nitro group),
one hydrogen too many (``[NH4][Pt]``), or an ionic charge on an atom that is
also drawn bonded (``[Cl-][Pt]``).  This module finds the smallest per-atom
edit that makes every atom's valence possible and returns the original string
with only those atoms rewritten, so stereo marks, atom order and formatting
are kept.

Some atoms can be fixed two ways: an over-valent ``=[NH2]`` becomes valid
either as ``=[NH2+]`` (add a charge) or ``=[NH]`` (drop a hydrogen).  The
string alone cannot tell which was meant.  With ``physiological_ph`` the
choice is the protonation state at body pH (7.4), decided from a small table
of conjugate-acid pKa values (:data:`PKA_TABLE`); without it the written
charge is kept and the hydrogen is dropped.
"""

from __future__ import annotations

import copy
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from chem.adapter import (
    AtomDescriptor,
    BondDescriptor,
    _NOBLE_GAS_ELECTRONS,
    valence_violations,
)

PHYSIOLOGICAL_PH = 7.4

#: Approximate pKa of the conjugate acid of each nitrogen class the repair can
#: meet.  A nitrogen whose class pKa is above :data:`PHYSIOLOGICAL_PH` is
#: protonated in the body.
PKA_TABLE: Dict[str, float] = {
    "guanidine": 13.6,
    "amidine": 12.4,
    "aliphatic amine": 10.6,
    "ammonia": 9.25,
    "imine": 7.0,
    "pyridine-type aromatic N": 5.2,
    "aniline-type N": 4.6,
    "N-cyano/acyl/sulfonyl amidine or guanidine": 0.0,
    "amide/sulfonamide/nitroso N": -1.0,
    "metal-coordinated N": -10.0,
}

# Elements that commonly form "-ate" complexes ([B-](F)(F)(F)F, [Al-]).
_ATE_FORMERS = {"B", "Al", "Ga", "In"}


@dataclass
class Repair:
    """A repaired SMILES and the edits that produced it.

    Attributes:
        smiles: The input string with only the edited atoms rewritten.
        edits: One human-readable line per edited atom.
    """

    smiles: str
    edits: List[str] = field(default_factory=list)


# --------------------------------------------------------------------------- #
# Source positions of atoms
# --------------------------------------------------------------------------- #

def atom_spans(smiles: str) -> List[Tuple[int, int]]:
    """Return the ``(start, end)`` character span of every atom, in source order.

    Atom order matches the descriptor order from ``to_descriptor`` (the parser
    registers atoms in source order).
    """
    from syntax.lex import SmilesLex

    spans: List[Tuple[int, int]] = []
    bracket_start: Optional[int] = None
    for token in SmilesLex().tokenize(smiles):
        if token.type == "[":
            bracket_start = token.index
        elif token.type == "]":
            spans.append((bracket_start, token.index + 1))
            bracket_start = None
        elif bracket_start is None and token.type in ("semi_symbol", "*"):
            spans.append((token.index, token.index + len(str(token.value))))
    return spans


_BRACKET = re.compile(
    r"^\[(?P<isotope>\d*)(?P<symbol>\*|[A-Za-z][a-z]?)(?P<chiral>@[^H+\-:\]]*)?"
    r"(?P<h>H\d*)?(?P<charge>[+-]+\d*)?(?P<cls>:\d+)?\]$"
)


def _charge_text(charge: int) -> str:
    if charge == 0:
        return ""
    sign = "+" if charge > 0 else "-"
    return sign if abs(charge) == 1 else f"{sign}{abs(charge)}"


def _h_text(hcount: int) -> str:
    if hcount == 0:
        return ""
    return "H" if hcount == 1 else f"H{hcount}"


# Organic-subset elements (written without brackets) and their normal
# valences, used to decide whether an edited atom can drop its brackets.
_ORGANIC_VALENCES = {
    "B": (3,), "C": (4,), "N": (3, 5), "O": (2,), "P": (3, 5), "S": (2, 4, 6),
    "F": (1,), "Cl": (1,), "Br": (1,), "I": (1,),
}


def _rewrite_atom(text: str, hcount: int, charge: int, bond_total: float) -> str:
    """Rewrite one atom token with a new hydrogen count and charge.

    A neutral, unlabelled, non-aromatic atom whose hydrogen count equals what
    the organic subset would imply is written bare (``[N]`` -> ``N``).
    """
    match = _BRACKET.match(text)
    if match is None:  # organic-subset atom
        parts = {"isotope": "", "symbol": text, "chiral": "", "cls": ""}
    else:
        parts = match.groupdict(default="")
    symbol = parts["symbol"]
    if (
        charge == 0
        and not (parts["isotope"] or parts["chiral"] or parts["cls"])
        and symbol in _ORGANIC_VALENCES
    ):
        implied = next((v for v in _ORGANIC_VALENCES[symbol] if v >= bond_total), None)
        if implied is not None and implied - bond_total == hcount:
            return symbol
    return (
        f"[{parts['isotope']}{symbol}{parts['chiral']}"
        f"{_h_text(hcount)}{_charge_text(charge)}{parts['cls']}]"
    )


# --------------------------------------------------------------------------- #
# Protonation state at physiological pH
# --------------------------------------------------------------------------- #

def _neighbours(bonds: List[BondDescriptor], n: int) -> List[Tuple[int, float]]:
    return [(j if i == n else i, o) for i, j, o in bonds if n in (i, j)]


def _is_withdrawing_neighbour(atoms, bonds, n: int, via: int) -> bool:
    """Whether atom ``via`` (bonded to N ``n``) withdraws N's lone pair.

    True for cyano carbon, acyl/sulfonyl carbon or sulfur (C=O, S=O), and
    nitroso/nitro nitrogen (N=O).
    """
    for other, order in _neighbours(bonds, via):
        if other == n:
            continue
        element = atoms[other]["element"]
        if order == 3 and element == "N":
            return True
        if order == 2 and element == "O":
            return True
    return False


def nitrogen_class(atoms, bonds, n: int) -> str:
    """Classify a nitrogen into a :data:`PKA_TABLE` row.

    Rules are deliberately few and ordered from most to least specific; they
    cover the over-valent nitrogens found in public datasets (guanidines,
    amidines, ammonia and amines bonded to metals).
    """
    from chem.adapter import _is_d_f_block

    neighbours = _neighbours(bonds, n)
    elements = [atoms[j]["element"] for j, _ in neighbours]
    if any(e != "*" and _is_d_f_block(e) for e in elements):
        return "metal-coordinated N"
    if atoms[n]["aromatic"]:
        return "pyridine-type aromatic N"
    if not neighbours:
        return "ammonia"

    double_c = [j for j, o in neighbours if o == 2 and atoms[j]["element"] == "C"]
    if double_c:
        centre = double_c[0]
        system_n = [j for j, _ in _neighbours(bonds, centre) if atoms[j]["element"] == "N"]
        if len(system_n) >= 2:
            for system_atom in system_n:
                for other, _ in _neighbours(bonds, system_atom):
                    if other != centre and _is_withdrawing_neighbour(atoms, bonds, system_atom, other):
                        return "N-cyano/acyl/sulfonyl amidine or guanidine"
            return "guanidine" if len(system_n) >= 3 else "amidine"
        return "imine"

    for j, _ in neighbours:
        if _is_withdrawing_neighbour(atoms, bonds, n, j):
            return "amide/sulfonamide/nitroso N"
        if atoms[j]["aromatic"]:
            return "aniline-type N"
        if any(o >= 2 for _, o in _neighbours(bonds, j)):
            return "aniline-type N"  # conjugated with a C=C / C=N
    return "aliphatic amine"


def protonated_at_ph(atoms, bonds, n: int, ph: float = PHYSIOLOGICAL_PH) -> Tuple[bool, str]:
    """Return ``(is_protonated, explanation)`` for nitrogen ``n`` at ``ph``."""
    kind = nitrogen_class(atoms, bonds, n)
    pka = PKA_TABLE[kind]
    state = "protonated" if pka > ph else "neutral"
    return pka > ph, f"{kind}, pKa ≈ {pka:g}: {state} at pH {ph:g}"


# --------------------------------------------------------------------------- #
# Repair search
# --------------------------------------------------------------------------- #

def _atom_ok(atoms, bonds, n: int) -> bool:
    try:
        return all(node != n for node, _ in valence_violations(atoms, bonds))
    except Exception:  # noqa: BLE001 - kekulization failure: not a fix
        return False


def _closed_shell_ion(atom: AtomDescriptor) -> bool:
    from pysmiles import PTE

    charge = atom.get("charge", 0)
    return bool(charge) and PTE[atom["element"]]["AtomicNumber"] - charge in _NOBLE_GAS_ELECTRONS


def _label(atoms, n: int) -> str:
    return f"{atoms[n]['element']} (atom {n})"


def suggest_repair(
    smiles: str,
    atoms: List[AtomDescriptor],
    bonds: List[BondDescriptor],
    physiological_ph: bool = True,
) -> Optional[Repair]:
    """Find the smallest per-atom edits that make an invalid molecule valid.

    Args:
        smiles: The original SMILES (edits are applied to this string).
        atoms: Atom descriptors from ``to_descriptor``.
        bonds: Bond descriptors from ``to_descriptor``.
        physiological_ph: Choose between "add a charge" and "drop a hydrogen"
            by the protonation state at pH 7.4.

    Returns:
        A :class:`Repair`, or ``None`` when the molecule is already valid, the
        failure is not a valence problem (e.g. an unkekulizable aromatic
        ring), or no single-atom edit fixes it.
    """
    try:
        violations = valence_violations(atoms, bonds)
    except Exception:  # noqa: BLE001 - aromaticity failure: no valence repair
        return None
    if not violations:
        return None

    work = copy.deepcopy(atoms)
    edits: List[str] = []
    changed = set()

    # 1. A closed-shell ion drawn bonded ([Cl-][Pt], [F-][B+3]): the bond is
    #    real, so the charge belongs on the atom it is bonded to.
    for n, _ in violations:
        neighbours = _neighbours(bonds, n)
        if _closed_shell_ion(work[n]) and len(neighbours) == 1:
            target = neighbours[0][0]
            charge = work[n]["charge"]
            work[target]["charge"] = work[target].get("charge", 0) + charge
            work[n]["charge"] = 0
            changed.update((n, target))
            edits.append(
                f"{_label(work, n)} is a closed-shell ion but is drawn bonded: moved "
                f"its {charge:+d} charge onto {_label(work, target)}"
            )

    # 2. Remaining over-valent atoms: drop hydrogens or change the charge.
    for n, reason in valence_violations(work, bonds):
        candidates: List[Tuple[str, AtomDescriptor]] = []
        hcount = work[n].get("hcount", 0)
        for k in range(1, hcount + 1):
            trial = copy.deepcopy(work)
            trial[n]["hcount"] = hcount - k
            if _atom_ok(trial, bonds, n):
                candidates.append((f"removed {k} H", trial[n]))
                break
        for delta in (+1, -1):
            if delta < 0 and work[n]["element"] not in _ATE_FORMERS:
                continue
            trial = copy.deepcopy(work)
            trial[n]["charge"] = work[n].get("charge", 0) + delta
            if _atom_ok(trial, bonds, n):
                candidates.append((f"set charge {trial[n]['charge']:+d}", trial[n]))

        if not candidates:
            return None
        choice = candidates[0]
        note = ""
        charge_options = [c for c in candidates if c[0].startswith("set charge")]
        h_options = [c for c in candidates if c[0].startswith("removed")]
        if charge_options and h_options:
            if physiological_ph and work[n]["element"] == "N":
                protonated, why = protonated_at_ph(work, bonds, n)
                choice = charge_options[0] if protonated else h_options[0]
                note = f" ({why})"
            else:
                choice = h_options[0]
                note = f" (alternative: {charge_options[0][0]})"
        work[n] = choice[1]
        changed.add(n)
        edits.append(f"{_label(work, n)}: {reason.split(': ', 1)[-1]}; {choice[0]}{note}")

    try:
        if valence_violations(work, bonds):
            return None
    except Exception:  # noqa: BLE001
        return None

    return Repair(smiles=_apply(smiles, atoms, work, bonds, changed), edits=edits)


def _apply(smiles, atoms, work, bonds, changed) -> str:
    """Rewrite the changed atoms in the original string."""
    from chem.adapter import _perceive

    spans = atom_spans(smiles)
    filled = _perceive(work, bonds)  # implicit H of edited organic atoms
    out = smiles
    for n in sorted(changed, key=lambda i: spans[i][0], reverse=True):
        start, end = spans[n]
        hcount = work[n].get("hcount", filled.nodes[n].get("hcount", 0))
        bond_total = sum(o for _, o in _neighbours(bonds, n))
        token = _rewrite_atom(smiles[start:end], hcount, work[n].get("charge", 0), bond_total)
        out = out[:start] + token + out[end:]
    return out
