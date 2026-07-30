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


# OpenSMILES features the grammar used to reject (the "TOO STRICT" rows of the
# conformance gap table).  Each is tagged with the spec section that makes it
# valid.  [C@@TH2] is valid per §3.8 even though RDKit rejects it — the
# intentional-disagreement rationale lives in tests/conformance/.
OPENSMILES_ACCEPTS = [
    ("[CH4:1234]", "§3.10 atom class, 4-digit"),
    ("C:C", "§3.2 aromatic bond symbol"),
    ("c1:c:c:c:c:c:1", "§3.2 explicit aromatic bonds in ring"),
    ("[C@TH1](F)(Cl)(Br)I", "§3.8 tetrahedral chirality class"),
    ("[C@@TH2](F)(Cl)(Br)I", "§3.8 chirality (RDKit rejects; spec accepts)"),
    ("[S@AL1](F)(Cl)(Br)I", "§3.8 allenal chirality class"),
    ("[Pt@SP1](F)(Cl)(Br)I", "§3.8 square-planar chirality class"),
    ("[Co@TB15](F)(Cl)(Br)(I)(O)N", "§3.8 trigonal-bipyramidal chirality class"),
    ("[Co@OH25](F)(Cl)(Br)(I)(O)N", "§3.8 octahedral chirality class"),
]


@pytest.mark.parametrize("smiles,spec_ref", OPENSMILES_ACCEPTS)
def test_opensmiles_features_accepted(smiles, spec_ref):
    result, error = validate_smiles(smiles)
    assert result is True, f"{smiles} is valid per OpenSMILES {spec_ref}, got: {error}"
