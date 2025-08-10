from fastapi import FastAPI
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

app = FastAPI()

@app.get("/smiles-to-mol")
def smiles_to_mol(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": "Invalid SMILES"}
    formula = CalcMolFormula(mol)
    return {"formula": formula}

