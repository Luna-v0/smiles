"""Build an Excel workbook focused on RDKit and rd-yacc-kit: one sheet per
validator listing every molecule it marks INVALID together with the *reason*,
plus a final sheet of the molecules where the two disagree.

rd-yacc-kit uses RDKit's chemistry on the graph OUR parser builds, so
RDKit-vs-rd-yacc-kit disagreements isolate parser/graph effects from chemistry.

Universe: every molecule in data/*.parquet (~2.25M, from the cached results)
plus the labelled-invalid set. SMILES are re-attached to the cache by per-dataset
alignment (verified by a spot-check).

Output: analyses/validator_invalids.xlsx
"""
import os
import sys

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO)
sys.path.insert(0, os.path.join(_REPO, "src"))

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from src import parse_smiles
from chem.adapter import to_descriptor, valid_rdkit

FULL = os.path.join(_REPO, "analyses", "full_validation_results.parquet")
INV = os.path.join(_REPO, "analyses", "invalid_validation_results.parquet")
OUT = os.path.join(_REPO, "analyses", "validator_invalids.xlsx")


def rdkit_reason(smi):
    """Why RDKit rejects `smi` ('' if it accepts)."""
    if Chem.MolFromSmiles(smi) is not None:
        return ""
    relaxed = Chem.MolFromSmiles(smi, sanitize=False)
    if relaxed is None:
        return "SMILES parse error (syntax)"
    try:
        Chem.SanitizeMol(relaxed)
        return "rejected by full parse"
    except Exception as exc:  # noqa: BLE001
        return str(exc)


def rd_yacc_kit_reason(smi):
    """Why rd-yacc-kit rejects `smi` ('' if it accepts)."""
    ok, g, exc = parse_smiles(smi)
    if not ok:
        return f"syntax (our grammar): {exc}"
    atoms, bonds = to_descriptor(g)
    valid, reason = valid_rdkit(atoms, bonds)
    if valid:
        return ""
    return f"chemistry (RDKit sanitize): {reason}".replace("rdkit sanitize failed: ", "")


# ---- assemble the molecule universe with SMILES + the two verdicts ----
cache = pd.read_parquet(FULL).reset_index(drop=True)
cache["dataset"] = cache["dataset"].astype(str)
order = cache["dataset"].drop_duplicates().tolist()
smiles = pd.Series(index=cache.index, dtype=object)
for name in order:
    idxs = cache.index[cache["dataset"] == name]
    df = pd.read_parquet(os.path.join(_REPO, "data", f"{name}.parquet"))
    col = "smiles" if "smiles" in df.columns else ("SMILES" if "SMILES" in df.columns else
           [c for c in df.columns if "smi" in c.lower()][0])
    s = df[col].dropna().astype(str).str.strip().tolist()
    assert len(s) == len(idxs), f"alignment mismatch {name}"
    smiles.loc[idxs] = s
cache["smiles"] = smiles.values

# spot-check
chk = cache.sample(min(400, len(cache)), random_state=0)
mism = sum((Chem.MolFromSmiles(r.smiles) is not None) != bool(r.RDKit) for r in chk.itertuples())
print(f"alignment spot-check: {mism}/{len(chk)} mismatches")
assert mism == 0

real = cache[["dataset", "smiles", "RDKit", "rd-yacc-kit"]]
inv = pd.read_parquet(INV)
inv_min = inv[["smiles", "RDKit", "rd-yacc-kit"]].copy()
inv_min.insert(0, "dataset", "invalid_smiles")
allm = pd.concat([real, inv_min], ignore_index=True)

# compute reasons only for molecules rejected by at least one of the two
need = allm[(~allm["RDKit"]) | (~allm["rd-yacc-kit"])].copy()
print(f"computing reasons for {len(need)} rejected molecules...")
need["RDKit_reason"] = [rdkit_reason(s) if not v else "" for s, v in zip(need.smiles, need.RDKit)]
need["rd-yacc-kit_reason"] = [rd_yacc_kit_reason(s) if not v else "" for s, v in zip(need.smiles, need["rd-yacc-kit"])]

rdkit_sheet = (need[~need["RDKit"]][["dataset", "smiles", "RDKit_reason", "rd-yacc-kit"]]
               .rename(columns={"rd-yacc-kit": "rd-yacc-kit_verdict"}).reset_index(drop=True))
rdyk_sheet = (need[~need["rd-yacc-kit"]][["dataset", "smiles", "rd-yacc-kit_reason", "RDKit"]]
              .rename(columns={"RDKit": "RDKit_verdict"}).reset_index(drop=True))
disagree = (need[need["RDKit"] != need["rd-yacc-kit"]]
            [["dataset", "smiles", "RDKit", "RDKit_reason", "rd-yacc-kit", "rd-yacc-kit_reason"]]
            .reset_index(drop=True))

readme = pd.DataFrame({
    "sheet": ["RDKit", "rd-yacc-kit", "disagreements"],
    "what": [
        f"Every molecule RDKit marks INVALID, with RDKit's reason. ({len(rdkit_sheet)} rows)",
        f"Every molecule rd-yacc-kit marks INVALID, with the reason. ({len(rdyk_sheet)} rows)",
        f"Molecules where RDKit and rd-yacc-kit disagree (= parser/graph effects). ({len(disagree)} rows)",
    ],
    "note": [
        "reason = RDKit SanitizeMol message, or 'SMILES parse error' for syntax.",
        "reason is tagged 'syntax (our grammar)' or 'chemistry (RDKit sanitize)'.",
        "rd-yacc-kit uses RDKit chemistry, so differences come from our parser/graph.",
    ],
})

with pd.ExcelWriter(OUT, engine="openpyxl") as xl:
    readme.to_excel(xl, sheet_name="README", index=False)
    rdkit_sheet.to_excel(xl, sheet_name="RDKit", index=False)
    rdyk_sheet.to_excel(xl, sheet_name="rd-yacc-kit", index=False)
    disagree.to_excel(xl, sheet_name="disagreements", index=False)

print(f"\nwrote {OUT}")
print(f"  RDKit invalid:       {len(rdkit_sheet)}")
print(f"  rd-yacc-kit invalid: {len(rdyk_sheet)}")
print(f"  disagreements:       {len(disagree)}")
print("\n--- disagreements (RDKit vs rd-yacc-kit) ---")
with pd.option_context("display.max_colwidth", 40, "display.width", 200):
    print(disagree.to_string())
print("\n--- reason breakdown ---")
print("RDKit reasons:", dict(rdkit_sheet["RDKit_reason"].str[:40].value_counts().head(6)))
print("rd-yacc-kit reasons:", dict(rdyk_sheet["rd-yacc-kit_reason"].str[:40].value_counts().head(6)))
