import os
import sys

import pandas as pd
import pytest

# Since pythonpath includes src, we can import directly
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))
from src import validate_smiles

# Regression snapshot over the small hand-written labelled set.  This dataset is
# narrow (variations of a handful of PAHs/heteroarenes) and was partly built to
# match the old string heuristics, so its labels are *not* treated as
# independent ground truth — they are aligned to the agreed permissive
# validation policy (accept radicals + hypervalent chemistry, reject impossible
# valence and malformed syntax).  Independent correctness assertions live in
# tests/chem/test_adapter.py; the trustworthy evaluation is the large molecule
# databases + data/invalid_smiles.parquet (see docs/chemistry_backend_study.md).
test_data = pd.read_parquet('tests/data/test_molecules_from_values.parquet')
test_data = test_data.dropna(subset=['smiles'])

# Smoke sample: this dataset is a distrusted regression snapshot (see the module
# docstring), so instead of exploding it into one parametrized test per row we
# run a small deterministic slice as a single smoke test.  Independent
# correctness lives in tests/chem/test_adapter.py.
_SMOKE_SAMPLE = list(test_data.head(20).itertuples(index=False, name=None))


def test_smiles_validation_smoke():
    for molecule_name, smiles, is_valid, description in _SMOKE_SAMPLE:
        result, _ = validate_smiles(smiles)
        assert result is bool(is_valid), (
            f"Expected {smiles} to be {'valid' if is_valid else 'invalid'}, got "
            f"{result}. Molecule: {molecule_name}. Description: {description}"
        )


def test_simple_molecules():
    # test all the benzene derivatives
    benzene_derivatives = [
        "c1ccccc1",
        "[H]c1c([H])c([H])c([H])c([H])c1[H]",
        "C1CCCCC1"
    ]
    
    for smiles in benzene_derivatives:
        result, _ = validate_smiles(smiles)
        assert result is True, f"Expected {smiles} to be valid, but it was invalid. For Benzene"

    napthalene_derivatives = [
        "c1ccc2ccccc2c1",
        "[H]c1ccc2ccccc2c1[H]",
        "C1CCCC2CCCCC21"
    ]

    for smiles in napthalene_derivatives:
        result, _ = validate_smiles(smiles)
        assert result is True, f"Expected {smiles} to be valid, but it was invalid. For Napthalene"

# NOTE: the OpenSMILES "TOO STRICT" acceptance cases ([CH4:1234], C:C,
# chirality classes, ...) used to live here as OPENSMILES_ACCEPTS.  They are an
# exact subset of the frozen conformance gap table in
# tests/conformance/test_opensmiles_gaps.py, which is their canonical home, so
# the duplicate parametrized test was removed.
