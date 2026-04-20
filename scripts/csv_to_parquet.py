import sys
from pathlib import Path
import pandas as pd

CSV_FILES = [
    "data/clintox.csv",
    "data/invalid_smiles.csv",
    "data/moses_test.csv",
    "data/moses_train.csv",
    "data/mutagenicity.csv",
    "data/muv.csv",
    "data/qm9.csv",
    "data/sider.csv",
    "data/tox21.csv",
    "data/zinc250k.csv",
    "analyses/validation_results_full.csv",
    "tests/data/test_molecules_from_values.csv",
]

root = Path(__file__).resolve().parent.parent
for rel in CSV_FILES:
    src = root / rel
    dst = src.with_suffix(".parquet")
    if not src.exists():
        print(f"skip (missing): {rel}")
        continue
    df = pd.read_csv(src, low_memory=False)
    df.to_parquet(dst, engine="pyarrow", compression="snappy", index=False)
    before = src.stat().st_size
    after = dst.stat().st_size
    print(f"{rel}: {before/1e6:.2f} MB -> {after/1e6:.2f} MB ({after/before*100:.1f}%)")
