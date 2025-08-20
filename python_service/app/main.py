from fastapi import FastAPI
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import descriptors.qed as qed
import descriptors.all_descriptors as all_ds

app = FastAPI()

@app.get("/smiles_to_mol")
def smiles_to_mol(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": "Invalid SMILES"}
    formula = CalcMolFormula(mol)
    return {"formula": formula}


@app.get("/v1/calc_qed")
def calc_qed(smiles: str):
    return qed.calc_qed(smiles)

@app.get("/v1/descriptors")
def calc_descriptors(smiles: str):
    return all_ds.calc_all_descriptors(smiles)



