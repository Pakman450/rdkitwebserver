import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'app')))

from main import app
from fastapi.testclient import TestClient

client = TestClient(app)

def test_smiles_to_mol_valid():
    response = client.get("/smiles-to-mol", params={"smiles": "CCO"})
    assert response.status_code == 200
    assert response.json() == {"formula": "C2H6O"}

def test_smiles_to_mol_invalid():
    response = client.get("/smiles-to-mol", params={"smiles": "invalid"})
    assert response.status_code == 200
    assert "error" in response.json()
