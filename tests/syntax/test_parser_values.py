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


@pytest.mark.parametrize(
    "molecule_name,smiles,is_valid,description",
    list(test_data.itertuples(index=False, name=None)),
)
def test_smiles_validation(molecule_name, smiles, is_valid, description):
    result, _ = validate_smiles(smiles)
    if is_valid:
        assert result is True, f"Expected {smiles} to be valid, but it was invalid. Molecule: {molecule_name}"
    else:
        assert result is False, f"Expected {smiles} to be invalid, but it was valid. Molecule: {molecule_name}. Description: {description}"
