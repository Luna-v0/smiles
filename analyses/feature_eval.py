"""Feature evaluation harness for the three validator feature groups:
salt removal, charge/hydrogen repairs (incl. ``physiological_ph``), and the
grammar deviations from strict OpenSMILES.

This is analysis code only -- it imports ``src``/``chem`` read-only (the
same private helpers ``chem/repair.py`` and ``chem/fragments.py`` already
import from each other) and never edits library behaviour.

Stages (each resumable, cached under ``analyses/_feature_eval_parts/``):

    main        full pass over the 9 datasets: validity, tier, salts/repair
                reporting fields, syntax-deviation flags, RDKit agreement.
    labelled    the 52 labelled-invalid records under all 4 configs
                (remove_salts x physiological_ph).
    ph_subset   re-run of the repaired-and-invalid subset with
                physiological_ph=False, to see what changes.
    salts_sample  a sample re-run with remove_salts=False, to confirm it
                does not change the valid/invalid verdict.
    salts_agree  RDKit SaltRemover / LargestFragmentChooser agreement on
                every multi-component record from `main`.
    corrupt     the known-answer corruption -> repair recovery test.

Usage (from the repo root, with the shared CPU budget respected):

    nice -n 19 .venv/bin/python analyses/feature_eval.py main
    nice -n 19 .venv/bin/python analyses/feature_eval.py labelled
    nice -n 19 .venv/bin/python analyses/feature_eval.py ph_subset
    nice -n 19 .venv/bin/python analyses/feature_eval.py salts_sample
    nice -n 19 .venv/bin/python analyses/feature_eval.py salts_agree
    nice -n 19 .venv/bin/python analyses/feature_eval.py corrupt
    nice -n 19 .venv/bin/python analyses/feature_eval.py all      # all of the above, in order
"""
from __future__ import annotations

import argparse
import glob
import logging
import os
import random
import re
import sys
import time
from multiprocessing import get_context

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO)
sys.path.insert(0, os.path.join(_REPO, "src"))

logging.getLogger("pysmiles").setLevel(logging.CRITICAL)

import numpy as np
import pandas as pd
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

PARTS = os.path.join(_REPO, "analyses", "_feature_eval_parts")
MAX_PROCS = 4
CHUNK = 50_000
MOSES_TRAIN_SAMPLE = 200_000
SEED = 42

REAL_DATASETS = [
    "tox21", "mutagenicity", "clintox", "sider",  # drug-like: salts & charges
    "zinc250k", "qm9", "muv", "moses_test", "moses_train",
]

# -- strict-OpenSMILES organic subset (bare, unbracketed) -------------------
_STRICT_BARE = {
    "B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I",  # organic subset
    "b", "c", "n", "o", "p", "s",                          # aromatic organic subset
    "*",
}
_BRACKET_SPAN = re.compile(r"\[[^\]]*\]")


def _niceify(_=None):
    try:
        os.nice(19)
    except Exception:
        pass


# --------------------------------------------------------------------------- #
# Data loading
# --------------------------------------------------------------------------- #

def _smiles_column(df):
    for c in ("smiles", "SMILES"):
        if c in df.columns:
            return c
    cand = [c for c in df.columns if "smi" in c.lower()]
    return cand[0] if cand else None


def load_dataset(name: str) -> list[str]:
    """Load one dataset's SMILES column. ``moses_train`` is a seeded 200k sample
    (per the brief: full moses_train would not add signal over moses_test and
    would roughly triple total runtime for no new finding)."""
    path = os.path.join(_REPO, "data", f"{name}.parquet")
    df = pd.read_parquet(path)
    col = _smiles_column(df)
    s = df[col].dropna().astype(str).str.strip()
    s = s[s != ""]
    if name == "moses_train" and len(s) > MOSES_TRAIN_SAMPLE:
        s = s.sample(MOSES_TRAIN_SAMPLE, random_state=SEED)
    return s.tolist()


def load_invalid() -> pd.DataFrame:
    path = os.path.join(_REPO, "data", "invalid_smiles.parquet")
    df = pd.read_parquet(path)
    return df


# --------------------------------------------------------------------------- #
# Syntax-deviation flags (post-hoc, no library change: strict-grammar
# emulation lives entirely here in the analysis harness).
# --------------------------------------------------------------------------- #

def e1_flag(smiles: str) -> bool:
    """Bare (unbracketed) atom outside the strict OpenSMILES organic subset."""
    from syntax.lex import SmilesLex
    try:
        depth = 0
        for tok in SmilesLex().tokenize(smiles):
            if tok.type == "[":
                depth += 1
            elif tok.type == "]":
                depth -= 1
            elif depth == 0 and tok.type == "semi_symbol":
                if tok.value not in _STRICT_BARE:
                    return True
        return False
    except Exception:
        return False


def e2_e3_flags(smiles: str) -> tuple[bool, bool]:
    """Multi-digit hcount (E2) / zero-padded charge magnitude (E3) in any bracket."""
    from chem.repair import _BRACKET
    e2 = e3 = False
    for m in _BRACKET_SPAN.finditer(smiles):
        bm = _BRACKET.match(m.group(0))
        if bm is None:
            continue
        h = bm.group("h") or ""
        if re.fullmatch(r"H\d{2,}", h):
            e2 = True
        charge = bm.group("charge") or ""
        if re.fullmatch(r"[+-]0\d+", charge):
            e3 = True
    return e2, e3


def syntax_flags(smiles: str) -> dict:
    e2, e3 = e2_e3_flags(smiles)
    return {"e1": e1_flag(smiles), "e2": e2, "e3": e3}


# --------------------------------------------------------------------------- #
# RDKit checks
# --------------------------------------------------------------------------- #

def rdkit_valid(smiles: str) -> bool:
    try:
        from rdkit import Chem
        return Chem.MolFromSmiles(smiles) is not None
    except Exception:
        return False


def rdkit_canonical(smiles: str) -> str | None:
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return None
        return Chem.MolToSmiles(m)
    except Exception:
        return None


# --------------------------------------------------------------------------- #
# Stage: main pass
# --------------------------------------------------------------------------- #

def process_one(args):
    dataset, smiles = args
    from src import validate_smiles_detailed

    t0 = time.perf_counter()
    try:
        r = validate_smiles_detailed(smiles)
        crashed, err = False, ""
    except Exception as exc:  # internal error: never let a worker die
        crashed, err = True, f"{type(exc).__name__}: {exc}"
        r = None
    elapsed_us = (time.perf_counter() - t0) * 1e6

    if crashed:
        return {
            "dataset": dataset, "smiles": smiles, "crashed": True, "error_text": err,
            "valid": False, "tier": None, "n_components": 1, "removed_fragments": "",
            "failed_fragment": "", "parent_heavy_atoms": None, "repaired": False,
            "n_repairs": 0, "repair_kinds": "", "repaired_smiles": "",
            "cistrans_warn": False, "backslash_warn": False,
            "e1": False, "e2": False, "e3": False,
            "rdkit_valid": False, "rdkit_valid_repaired": None, "elapsed_us": elapsed_us,
        }

    flags = syntax_flags(smiles)
    warnings = r.warnings or []
    cistrans = any("cis/trans mark" in w for w in warnings)
    backslash = any("doubled backslash" in w for w in warnings)
    n_components = 1 + len(r.removed_fragments) if r.parent_smiles is not None else 1
    kinds = []
    for edit in r.repairs:
        if "removed" in edit and "H" in edit:
            kinds.append("drop_H")
        elif "set charge" in edit:
            kinds.append("add_charge")
        elif "moved" in edit:
            kinds.append("move_ion_charge")
        else:
            kinds.append("other")
    parent_heavy = None
    if r.parent_smiles is not None:
        parent_heavy = sum(
            1 for ch in re.findall(r"[A-Za-z]", r.parent_smiles) if ch.isupper() or ch in "cnops"
        )  # rough heavy-atom proxy, refined later from fragments if needed
    rdkv_rep = rdkit_valid(r.repaired_smiles) if r.repaired_smiles else None

    return {
        "dataset": dataset, "smiles": smiles, "crashed": False, "error_text": "",
        "valid": bool(r.valid), "tier": r.tier, "n_components": n_components,
        "removed_fragments": "|".join(r.removed_fragments), "failed_fragment": r.failed_fragment or "",
        "parent_heavy_atoms": parent_heavy, "repaired": r.repaired_smiles is not None,
        "n_repairs": len(r.repairs), "repair_kinds": ",".join(kinds),
        "repaired_smiles": r.repaired_smiles or "",
        "cistrans_warn": cistrans, "backslash_warn": backslash,
        "e1": flags["e1"], "e2": flags["e2"], "e3": flags["e3"],
        "rdkit_valid": rdkit_valid(smiles), "rdkit_valid_repaired": rdkv_rep,
        "elapsed_us": elapsed_us,
    }


def run_main(procs: int = MAX_PROCS):
    os.makedirs(PARTS, exist_ok=True)
    ctx = get_context("fork")
    for name in REAL_DATASETS:
        smiles_list = load_dataset(name)
        nchunks = (len(smiles_list) + CHUNK - 1) // CHUNK
        for ci in range(nchunks):
            part = os.path.join(PARTS, f"main__{name}__{ci:03d}.parquet")
            if os.path.exists(part):
                continue
            chunk = smiles_list[ci * CHUNK:(ci + 1) * CHUNK]
            t0 = time.time()
            with ctx.Pool(procs, initializer=_niceify, maxtasksperchild=2000) as pool:
                rows = pool.map(process_one, [(name, s) for s in chunk], chunksize=200)
            # Write to a temp name in the same dir, then atomically rename into place,
            # so a kill mid-write can never leave a partial file that os.path.exists()
            # treats as a finished chunk on resume.
            tmp_part = part + ".tmp"
            pd.DataFrame(rows).to_parquet(tmp_part, index=False)
            os.replace(tmp_part, part)
            dt = time.time() - t0
            print(f"[main] {name} chunk {ci+1}/{nchunks} ({len(chunk):,}) in {dt:.0f}s "
                  f"[{len(chunk)/dt:,.0f} mol/s]", flush=True)

    parts = sorted(glob.glob(os.path.join(PARTS, "main__*.parquet")))
    full = pd.concat((pd.read_parquet(p) for p in parts), ignore_index=True)
    out_path = os.path.join(_REPO, "analyses", "feature_eval_main.parquet")
    full.to_parquet(out_path, index=False)
    print(f"wrote {out_path} ({len(full):,} rows)", flush=True)


# --------------------------------------------------------------------------- #
# Stage: labelled invalids under all 4 configs
# --------------------------------------------------------------------------- #

def _run_config(smiles: str, remove_salts: bool, physiological_ph: bool):
    from src.validation import validate_smiles_detailed
    try:
        r = validate_smiles_detailed(smiles, remove_salts=remove_salts, physiological_ph=physiological_ph)
        return {
            "valid": bool(r.valid), "tier": r.tier, "repaired": r.repaired_smiles is not None,
            "repaired_smiles": r.repaired_smiles or "", "failed_fragment": r.failed_fragment or "",
            "n_components": 1 + len(r.removed_fragments) if r.parent_smiles is not None else 1,
            "crashed": False,
        }
    except Exception as exc:
        return {"valid": False, "tier": None, "repaired": False, "repaired_smiles": "",
                "failed_fragment": "", "n_components": 1, "crashed": True, "error": str(exc)}


def run_labelled():
    df = load_invalid()
    rows = []
    for _, rec in df.iterrows():
        smi = rec.get("smiles")
        cat = rec.get("error_category")
        smi_str = "" if smi is None else str(smi)
        for salts in (True, False):
            for ph in (True, False):
                res = _run_config(smi_str, salts, ph)
                rows.append({"smiles": smi_str, "error_category": cat,
                             "remove_salts": salts, "physiological_ph": ph, **res})
    out = pd.DataFrame(rows)
    out_path = os.path.join(_REPO, "analyses", "feature_eval_labelled.parquet")
    out.to_parquet(out_path, index=False)
    print(f"wrote {out_path} ({len(out)} rows)", flush=True)


# --------------------------------------------------------------------------- #
# Stage: physiological_ph=False re-run of the repaired-and-was-invalid subset
# --------------------------------------------------------------------------- #

def _ph_off_one(args):
    dataset, smiles, ph_on_repaired = args
    from src.validation import validate_smiles_detailed
    try:
        r = validate_smiles_detailed(smiles, physiological_ph=False)
        rep_off = r.repaired_smiles or ""
    except Exception:
        rep_off = ""
    return {
        "dataset": dataset, "smiles": smiles,
        "repaired_smiles_ph_on": ph_on_repaired, "repaired_smiles_ph_off": rep_off,
        "same_either_way": (ph_on_repaired == rep_off),
    }


def run_ph_subset(procs: int = MAX_PROCS):
    main_path = os.path.join(_REPO, "analyses", "feature_eval_main.parquet")
    main = pd.read_parquet(main_path)
    sub = main[(~main["valid"]) & (main["repaired"])][["dataset", "smiles", "repaired_smiles"]]
    ctx = get_context("fork")
    args = list(sub.itertuples(index=False, name=None))
    with ctx.Pool(procs, initializer=_niceify) as pool:
        rows = pool.map(_ph_off_one, args, chunksize=100)
    out = pd.DataFrame(rows)
    out_path = os.path.join(_REPO, "analyses", "feature_eval_ph_subset.parquet")
    out.to_parquet(out_path, index=False)
    print(f"wrote {out_path} ({len(out)} rows); differ={len(out) - out['same_either_way'].sum()}", flush=True)


# --------------------------------------------------------------------------- #
# Stage: remove_salts=False sample, to confirm it doesn't change the verdict
# --------------------------------------------------------------------------- #

def _salts_off_one(args):
    dataset, smiles, valid_on = args
    from src.validation import validate_smiles_detailed
    try:
        r = validate_smiles_detailed(smiles, remove_salts=False)
        valid_off = bool(r.valid)
    except Exception:
        valid_off = None
    return {"dataset": dataset, "smiles": smiles, "valid_salts_on": valid_on,
            "valid_salts_off": valid_off, "changed": valid_off is not None and valid_off != valid_on}


def run_salts_sample(n: int = 50_000, procs: int = MAX_PROCS):
    main_path = os.path.join(_REPO, "analyses", "feature_eval_main.parquet")
    main = pd.read_parquet(main_path)
    sample = main.sample(min(n, len(main)), random_state=SEED)[["dataset", "smiles", "valid"]]
    ctx = get_context("fork")
    args = list(sample.itertuples(index=False, name=None))
    with ctx.Pool(procs, initializer=_niceify) as pool:
        rows = pool.map(_salts_off_one, args, chunksize=200)
    out = pd.DataFrame(rows)
    out_path = os.path.join(_REPO, "analyses", "feature_eval_salts_sample.parquet")
    out.to_parquet(out_path, index=False)
    print(f"wrote {out_path} ({len(out)} rows); changed={out['changed'].sum()}", flush=True)


# --------------------------------------------------------------------------- #
# Stage: RDKit SaltRemover / LargestFragmentChooser agreement
# --------------------------------------------------------------------------- #

def _salts_agree_one(args):
    dataset, smiles, parent_smiles, removed = args
    from rdkit import Chem
    from rdkit.Chem import SaltRemover
    from rdkit.Chem.MolStandardize import rdMolStandardize

    out = {"dataset": dataset, "smiles": smiles, "our_parent": parent_smiles,
           "rdkit_parse_ok": False, "lfc_agrees": None, "lfc_parent_canon": None,
           "our_parent_canon": None, "salt_remover_removed_n": None}
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return out
    out["rdkit_parse_ok"] = True
    try:
        lfc = rdMolStandardize.LargestFragmentChooser()
        chosen = lfc.choose(mol)
        lfc_canon = Chem.MolToSmiles(chosen)
        out["lfc_parent_canon"] = lfc_canon
    except Exception:
        lfc_canon = None
    our_parent_mol = Chem.MolFromSmiles(parent_smiles) if parent_smiles else None
    our_canon = Chem.MolToSmiles(our_parent_mol) if our_parent_mol is not None else None
    out["our_parent_canon"] = our_canon
    if lfc_canon is not None and our_canon is not None:
        out["lfc_agrees"] = (lfc_canon == our_canon)
    try:
        remover = SaltRemover.SaltRemover()
        stripped = remover.StripMol(mol, dontRemoveEverything=True)
        out["salt_remover_removed_n"] = mol.GetNumAtoms() - stripped.GetNumAtoms()
    except Exception:
        pass
    return out


def _fragment_shape_one(args):
    dataset, smiles = args
    from src.validation import validate_smiles_detailed, _parse
    from chem.adapter import to_descriptor
    from chem.fragments import split_fragments, parent_index

    out = {"dataset": dataset, "smiles": smiles, "parent_smiles": "",
           "parent_heavy": None, "parent_carbons": None, "n_fragments": None,
           "max_other_heavy": None, "tie_on_heavy_atoms": False, "failed_fragment": ""}
    try:
        r = validate_smiles_detailed(smiles)
        out["parent_smiles"] = r.parent_smiles or ""
        out["failed_fragment"] = r.failed_fragment or ""
        graph, failure = _parse(smiles)
        if failure is not None:
            return out
        atoms, bonds = to_descriptor(graph)
        frags = split_fragments(smiles, atoms, bonds)
        if frags is None or len(frags) < 2:
            return out
        p = parent_index(frags)
        out["n_fragments"] = len(frags)
        out["parent_heavy"] = frags[p].heavy_atoms
        out["parent_carbons"] = frags[p].carbons
        others = [f.heavy_atoms for k, f in enumerate(frags) if k != p]
        out["max_other_heavy"] = max(others) if others else 0
        out["tie_on_heavy_atoms"] = out["max_other_heavy"] >= out["parent_heavy"]
    except Exception:
        pass
    return out


def run_salts_agree(procs: int = MAX_PROCS):
    main_path = os.path.join(_REPO, "analyses", "feature_eval_main.parquet")
    main = pd.read_parquet(main_path)
    multi = main[main["n_components"] >= 2].copy()
    ctx = get_context("fork")

    # Fragment shape (parent choice, ties, failed-fragment attribution).
    args_shape = list(multi[["dataset", "smiles"]].itertuples(index=False, name=None))
    with ctx.Pool(procs, initializer=_niceify) as pool:
        shape_rows = pool.map(_fragment_shape_one, args_shape, chunksize=100)
    shape = pd.DataFrame(shape_rows)
    shape_path = os.path.join(_REPO, "analyses", "feature_eval_fragment_shape.parquet")
    shape.to_parquet(shape_path, index=False)
    print(f"wrote {shape_path} ({len(shape)} rows); "
          f"ties={shape['tie_on_heavy_atoms'].sum()}", flush=True)

    # RDKit agreement on parent choice.
    args = list(shape[["dataset", "smiles", "parent_smiles", "smiles"]]
                .itertuples(index=False, name=None))
    args = [(d, s, p, "") for d, s, p, _ in args]
    with ctx.Pool(procs, initializer=_niceify) as pool:
        rows = pool.map(_salts_agree_one, args, chunksize=100)
    out = pd.DataFrame(rows)
    out_path = os.path.join(_REPO, "analyses", "feature_eval_salts_agree.parquet")
    out.to_parquet(out_path, index=False)
    print(f"wrote {out_path} ({len(out)} rows)", flush=True)


# --------------------------------------------------------------------------- #
# Stage: known-answer corruption -> repair recovery test
# --------------------------------------------------------------------------- #

def current_hcount(atoms, bonds, n):
    if "hcount" in atoms[n] and atoms[n]["hcount"] is not None:
        return atoms[n]["hcount"]
    from chem.adapter import _perceive
    mol = _perceive(atoms, bonds)
    return mol.nodes[n].get("hcount", 0)


def rewrite_single_atom(smiles, bonds, n, new_hcount, new_charge):
    from chem.repair import atom_spans, _rewrite_atom, _neighbours
    spans = atom_spans(smiles)
    start, end = spans[n]
    bond_total = sum(o for _, o in _neighbours(bonds, n))
    token = _rewrite_atom(smiles[start:end], new_hcount, new_charge, bond_total)
    return smiles[:start] + token + smiles[end:]


_HFLEX_ELEMENTS = {"C", "N", "O", "S", "P"}


def pick_corruption_targets(smiles, atoms, bonds, rng):
    """Return {kind: atom_index} for the corruption kinds applicable to this molecule."""
    out = {}
    charged = [n for n, a in enumerate(atoms) if a.get("charge", 0) != 0]
    if charged:
        out["strip_charge"] = rng.choice(charged)
    hflex_add = [n for n, a in enumerate(atoms)
                 if a.get("element") in _HFLEX_ELEMENTS and a.get("charge", 0) == 0]
    if hflex_add:
        out["add_H"] = rng.choice(hflex_add)
    hflex_rm = [n for n in hflex_add if current_hcount(atoms, bonds, n) >= 1]
    if hflex_rm:
        out["remove_H"] = rng.choice(hflex_rm)
    return out


def build_corruption(smiles, atoms, bonds, kind, n):
    if kind == "strip_charge":
        h = atoms[n].get("hcount") if "hcount" in atoms[n] else current_hcount(atoms, bonds, n)
        return rewrite_single_atom(smiles, bonds, n, h, 0)
    if kind == "add_H":
        h = current_hcount(atoms, bonds, n)
        return rewrite_single_atom(smiles, bonds, n, h + 1, atoms[n].get("charge", 0))
    if kind == "remove_H":
        h = current_hcount(atoms, bonds, n)
        return rewrite_single_atom(smiles, bonds, n, max(0, h - 1), atoms[n].get("charge", 0))
    raise ValueError(kind)


def _corrupt_one(args):
    pool_name, dataset, smiles, kind, seed = args
    from src.validation import validate_smiles_detailed
    from chem.adapter import to_descriptor
    from chem.structure import MolecularGraph  # noqa: F401  (import check only)

    rng = random.Random(seed)
    out = {"pool": pool_name, "dataset": dataset, "original_smiles": smiles, "kind": kind,
           "corrupted_smiles": None, "detected": None, "repaired": None,
           "exact_recovery_literal": None, "ph_alt_matches_original": None,
           "orig_charge_sign": None, "error": ""}
    try:
        from src.validation import _parse
        graph, failure = _parse(smiles)
        if failure is not None:
            out["error"] = "original did not parse"
            return out
        atoms, bonds = to_descriptor(graph)
        targets = pick_corruption_targets(smiles, atoms, bonds, rng)
        if kind not in targets:
            out["error"] = "no eligible atom"
            return out
        n = targets[kind]
        if kind == "strip_charge":
            c = atoms[n].get("charge", 0)
            out["orig_charge_sign"] = "positive" if c > 0 else ("negative" if c < 0 else "zero")
        corrupted = build_corruption(smiles, atoms, bonds, kind, n)
        out["corrupted_smiles"] = corrupted
    except Exception as exc:
        out["error"] = f"corruption build failed: {exc}"
        return out

    try:
        r_on = validate_smiles_detailed(corrupted, physiological_ph=True)
    except Exception as exc:
        out["error"] = f"validate crashed: {exc}"
        return out
    out["detected"] = not r_on.valid
    if not out["detected"]:
        return out  # corruption did not actually break the molecule
    out["repaired"] = r_on.repaired_smiles is not None
    if not out["repaired"]:
        return out

    orig_canon = rdkit_canonical(smiles)
    rep_on_canon = rdkit_canonical(r_on.repaired_smiles)
    out["exact_recovery_literal"] = (
        orig_canon is not None and rep_on_canon is not None and orig_canon == rep_on_canon
    )
    try:
        r_off = validate_smiles_detailed(corrupted, physiological_ph=False)
        rep_off_canon = rdkit_canonical(r_off.repaired_smiles) if r_off.repaired_smiles else None
    except Exception:
        rep_off_canon = None
    out["ph_alt_matches_original"] = (
        orig_canon is not None and rep_off_canon is not None and orig_canon == rep_off_canon
    )
    return out


def build_corruption_pools(target_n=3500, scan_cap=40):
    """Scan the drug-like + clean datasets for `target_n` charged / neutral,
    single-component, valid molecules each, seeded and deterministic."""
    from src.validation import validate_smiles_detailed, _parse
    from chem.adapter import to_descriptor

    charged_src = ["tox21", "clintox", "sider", "mutagenicity", "zinc250k"]
    neutral_src = ["moses_test", "qm9", "zinc250k"]
    rng = random.Random(SEED)

    def scan(sources, want_charged):
        pool = []
        scanned = 0
        for name in sources:
            smis = load_dataset(name)
            rng.shuffle(smis)
            for smi in smis:
                if len(pool) >= target_n or scanned >= target_n * scan_cap:
                    break
                scanned += 1
                if "." in smi:
                    continue
                graph, failure = _parse(smi)
                if failure is not None:
                    continue
                atoms, bonds = to_descriptor(graph)
                has_charge = any(a.get("charge", 0) != 0 for a in atoms)
                if has_charge != want_charged:
                    continue
                r = validate_smiles_detailed(smi)
                if not r.valid:
                    continue
                pool.append((name, smi))
            if len(pool) >= target_n:
                break
        return pool

    charged_pool = scan(charged_src, True)
    neutral_pool = scan(neutral_src, False)
    return charged_pool, neutral_pool


def run_corrupt(target_n=3500, procs: int = MAX_PROCS):
    charged_pool, neutral_pool = build_corruption_pools(target_n)
    print(f"[corrupt] charged pool={len(charged_pool)}  neutral pool={len(neutral_pool)}", flush=True)

    jobs = []
    for i, (name, smi) in enumerate(charged_pool):
        for kind in ("strip_charge", "add_H", "remove_H"):
            jobs.append(("charged", name, smi, kind, SEED * 1_000_003 + i))
    for i, (name, smi) in enumerate(neutral_pool):
        for kind in ("add_H", "remove_H"):
            jobs.append(("neutral", name, smi, kind, SEED * 7_000_033 + i))

    ctx = get_context("fork")
    with ctx.Pool(procs, initializer=_niceify, maxtasksperchild=2000) as pool:
        rows = pool.map(_corrupt_one, jobs, chunksize=50)
    out = pd.DataFrame(rows)
    out_path = os.path.join(_REPO, "analyses", "feature_eval_corruption.parquet")
    out.to_parquet(out_path, index=False)
    print(f"wrote {out_path} ({len(out)} rows)", flush=True)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=[
        "main", "labelled", "ph_subset", "salts_sample", "salts_agree", "corrupt", "all",
    ])
    ap.add_argument("--procs", type=int, default=MAX_PROCS)
    ap.add_argument("--corrupt-n", type=int, default=3500)
    ap.add_argument("--salts-sample-n", type=int, default=50_000)
    args = ap.parse_args()
    _niceify()

    stages = {
        "main": lambda: run_main(args.procs),
        "labelled": run_labelled,
        "ph_subset": lambda: run_ph_subset(args.procs),
        "salts_sample": lambda: run_salts_sample(args.salts_sample_n, args.procs),
        "salts_agree": lambda: run_salts_agree(args.procs),
        "corrupt": lambda: run_corrupt(args.corrupt_n, args.procs),
    }
    if args.stage == "all":
        for name in ["main", "labelled", "ph_subset", "salts_sample", "salts_agree", "corrupt"]:
            print(f"=== stage: {name} ===", flush=True)
            stages[name]()
    else:
        stages[args.stage]()


if __name__ == "__main__":
    main()
