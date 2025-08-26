
import pytest
from main import app
from fastapi.testclient import TestClient

import json

from lib.celery_worker import celery_app, calculate_descriptors_chunk
celery_app.conf.update(
    task_always_eager=True,
    task_eager_propagates=True  # propagate exceptions to the test
)


client = TestClient(app)


@pytest.fixture
def sample_file():
    with open("./examples/smiles.txt", "r") as f:
        data = f.readlines()  # read entire file as a single string
    
    lines = [line.rstrip("\n") for line in data] 
    return lines  # encode to bytes



def test_calculate_descriptors_chunk(sample_file):
    
    l_smiles = []
    for mol in sample_file:
        l_smiles.append(mol)

    # run the task synchronously
    result = calculate_descriptors_chunk(l_smiles)

    # reading in the made json file
    proccessed_lines = []
    with open(f"./{result}","r") as f:
        proccessed_lines = json.load(f)
    
    assert "json" in result and "None" in result
    assert len(proccessed_lines) == 10

    for i,line in enumerate(proccessed_lines):
        assert len(proccessed_lines[i].keys()) == 217


