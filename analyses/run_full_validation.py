"""Run every validator over every molecule in data/*.parquet (in parallel) and
cache the per-molecule results, so the analysis notebook can load them instantly.

Memory-safe and RESUMABLE: results are computed in chunks and each chunk is
written to analyses/_validation_parts/ immediately, so a kill (OOM, timeout)
loses at most one chunk -- just re-run to resume. Workers are recycled
(maxtasksperchild) to bound memory growth from the chemistry backends.

Usage (from the repo root):
    uv run python analyses/run_full_validation.py               # full ~2.25M
    uv run python analyses/run_full_validation.py --sample 2000 # quick test
    uv run python analyses/run_full_validation.py --procs 14

Outputs (under analyses/):
    full_validation_results.parquet     # one row per molecule: dataset + bool per validator
    invalid_validation_results.parquet  # the labelled-invalid set (with category)
"""
import argparse
import glob
import os
import re
import sys
import time
from multiprocessing import get_context

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO)
sys.path.insert(0, os.path.join(_REPO, "src"))

import logging
logging.getLogger("pysmiles").setLevel(logging.CRITICAL)

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
import partialsmiles as ps
import pysmiles
from molvs import validate_smiles as molvs_validate

from src import validate_smiles, parse_smiles, parser, lexer, parser_manager
from chem.validator import ChemistryValidator
from chem.adapter import to_descriptor, valid_rdkit

_oc = ChemistryValidator()

VALIDATORS = [
    "YACC (original)", "py-yacc-smiles", "rd-yacc-kit",
    "RDKit", "PartialSMILES", "PySMILES", "MolVS",
]

PARTS = os.path.join(_REPO, "analyses", "_validation_parts")
CHUNK = 150_000


def _yacc_original(m):
    try:
        if re.search(r"C$", m) and not re.search(r"C\d$", m) and re.search(r"\d", m):
            return False
        parser.parse(lexer.tokenize(m))
        g = parser_manager.graph_builder.get_graph()
        parser_manager.clear()
        ok, _ = _oc.validate_rings_and_valency(g)
        if not ok:
            return False
        ok, _ = _oc.validate_aromaticity(g)
        return bool(ok)
    except Exception:
        parser_manager.clear()
        return False


def _py_yacc_smiles(m):
    try:
        return bool(validate_smiles(m)[0])
    except Exception:
        return False


def _rd_yacc_kit(m):
    try:
        ok, g, _ = parse_smiles(m)
        if not ok:
            return False
        a, b = to_descriptor(g)
        return bool(valid_rdkit(a, b)[0])
    except Exception:
        return False


def _rdkit(m):
    try:
        return Chem.MolFromSmiles(m) is not None
    except Exception:
        return False


def _partialsmiles(m):
    try:
        ps.ParseSmiles(m, partial=False)
        return True
    except Exception:
        return False


def _pysmiles(m):
    try:
        pysmiles.read_smiles(m)
        return True
    except Exception:
        return False


def _molvs(m):
    try:
        return len(molvs_validate(m)) == 0
    except Exception:
        return False


def validate_all(smiles):
    m = smiles.strip()
    return (_yacc_original(m), _py_yacc_smiles(m), _rd_yacc_kit(m),
            _rdkit(m), _partialsmiles(m), _pysmiles(m), _molvs(m))


def _smiles_column(df):
    for c in ("smiles", "SMILES"):
        if c in df.columns:
            return c
    cand = [c for c in df.columns if "smi" in c.lower()]
    return cand[0] if cand else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=None)
    ap.add_argument("--procs", type=int, default=min(14, max(1, os.cpu_count() - 2)))
    ap.add_argument("--chunk", type=int, default=CHUNK)
    args = ap.parse_args()
    os.makedirs(PARTS, exist_ok=True)
    ctx = get_context("fork")

    files = sorted(glob.glob(os.path.join(_REPO, "data", "*.parquet")),
                   key=lambda p: os.path.getsize(p))  # smallest first
    inv_records = []
    for path in files:
        name = os.path.basename(path).replace(".parquet", "")
        df = pd.read_parquet(path)
        col = _smiles_column(df)
        if col is None:
            continue
        s = df[col].dropna().astype(str).str.strip()
        if name == "invalid_smiles":
            cats = df.get("error_category")
            inv_records = [(smi, str(cats.loc[i]) if cats is not None else "?")
                           for i, smi in zip(s.index, s)]
            continue
        if args.sample and len(s) > args.sample:
            s = s.sample(args.sample, random_state=42)
        s = s.tolist()
        nchunks = (len(s) + args.chunk - 1) // args.chunk
        for ci in range(nchunks):
            part = os.path.join(PARTS, f"{name}__{ci:03d}.parquet")
            if os.path.exists(part):
                continue
            chunk = s[ci * args.chunk:(ci + 1) * args.chunk]
            t0 = time.time()
            with ctx.Pool(args.procs, maxtasksperchild=2000) as pool:
                res = pool.map(validate_all, chunk, chunksize=500)
            out = pd.DataFrame(res, columns=VALIDATORS)
            out.insert(0, "dataset", name)
            out.to_parquet(part, index=False)
            dt = time.time() - t0
            print(f"{name} chunk {ci+1}/{nchunks} ({len(chunk):,}) in {dt:.0f}s "
                  f"[{len(chunk)/dt:,.0f} mol/s]", flush=True)

    # ---- concatenate all chunk parts ----
    parts = sorted(glob.glob(os.path.join(PARTS, "*.parquet")))
    full = pd.concat((pd.read_parquet(p) for p in parts), ignore_index=True)
    full["dataset"] = full["dataset"].astype("category")
    out_path = os.path.join(_REPO, "analyses", "full_validation_results.parquet")
    full.to_parquet(out_path, index=False)
    print(f"\nwrote {out_path}  ({len(full):,} rows from {len(parts)} chunks)", flush=True)

    # ---- invalid set ----
    if inv_records:
        inv_smiles = [r[0] for r in inv_records]
        with ctx.Pool(args.procs) as pool:
            inv_res = pool.map(validate_all, inv_smiles, chunksize=8)
        inv = pd.DataFrame(inv_res, columns=VALIDATORS)
        inv.insert(0, "smiles", inv_smiles)
        inv.insert(1, "error_category", [r[1] for r in inv_records])
        inv_path = os.path.join(_REPO, "analyses", "invalid_validation_results.parquet")
        inv.to_parquet(inv_path, index=False)
        print(f"wrote {inv_path}  ({len(inv)} rows)", flush=True)


if __name__ == "__main__":
    main()
