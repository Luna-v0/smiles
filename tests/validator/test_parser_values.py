from validator.yacc import validate_smiles

def test_benzene():
    assert validate_smiles("c1ccccc1")[0] == True
    assert validate_smiles("C1=CC=CC=C1")[0] == True
    assert validate_smiles("C1=CC=C=C1C")[0] == True
    assert validate_smiles("[H]c1c([H])c([H])c([H])c([H])c1[H]")[0] == True
    assert validate_smiles("c1ccccc1C")[0] == False
    assert validate_smiles("c1[c]cccc1")[0] == True
