from fastapi import FastAPI, Query
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import descriptors.qed as qed
import descriptors.all_descriptors as all_ds


#Types
from typing import List

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
def calc_descriptors(smiles: List[str] = Query(...)) -> List:
    print(smiles)
    return all_ds.calc_all_descriptors(smiles)



