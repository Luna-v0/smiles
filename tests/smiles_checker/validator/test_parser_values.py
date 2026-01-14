import pytest
import pandas as pd
# Since pythonpath includes src, we can import directly
import sys
import os
# Add src to path for this import
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))
from src import validate_smiles

# Load test data from CSV
test_data = pd.read_csv('tests/smiles_checker/validator/test_molecules_from_values.csv')

@pytest.mark.parametrize("molecule_name,smiles,is_valid,description", test_data.values)
def test_smiles_validation(molecule_name, smiles, is_valid, description):
    result, _ = validate_smiles(smiles)
    if is_valid:
        assert result is True, f"Expected {smiles} to be valid, but it was invalid. Molecule: {molecule_name}"
    else:
        assert result is False, f"Expected {smiles} to be invalid, but it was valid. Molecule: {molecule_name}. Description: {description}"
