from validator.yacc import validate_smiles


def test_benzene():
    assert validate_smiles("c1ccccc1") == (
        True,
        None,
    ), "Failed to validate aromatic ring with lowercase syntax"
    assert validate_smiles("C1=CC=CC=C1") == (
        True,
        None,
    ), "Failed to validate aromatic ring with uppercase syntax"
    assert validate_smiles("C1=CC=C=C1C") == (
        True,
        None,
    ), "Failed to validate aromatic ring with uppercase syntax with mapping at the end"
    assert validate_smiles("[H]c1c([H])c([H])c([H])c([H])c1[H]") == (
        True,
        None,
    ), "Failed to validate aromatic ring with bracket hydrogens"
    assert (
        validate_smiles("c1ccccc1C")[0] == False
    ), "Cannot join aromatic carbon ring with a Non-aromatic carbon"
    assert validate_smiles("c1[c]cccc1") == (
        True,
        None,
    ), "Failed to validate aromatic ring with lowercase syntax and bracketed carbon"


def test_naftalina():
    assert validate_smiles("c1ccc2ccccc2c1") == (
        True,
        None,
    ), "Failed to validate naphthalene structure"
    assert validate_smiles("C1=CC=CC=C1C2=CC=CC=C2") == (
        True,
        None,
    ), "Failed to validate naphthalene structure with uppercase syntax"
    assert validate_smiles("c1ccccc1C2=CC=CC=C2") == (
        True,
        None,
    ), "Failed to validate naphthalene structure with aromatic carbon and non-aromatic carbon"
