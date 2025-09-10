
import pytest
from main import app
from fastapi.testclient import TestClient

import json
from typing import List

from lib.celery_worker import celery_app, \
                              calculate_descriptors_chunk_smi, \
                              calculate_descriptors_chunk_sdf
celery_app.conf.update(
    task_always_eager=True,
    task_eager_propagates=True  # propagate exceptions to the test
)


client = TestClient(app)


@pytest.fixture
def sample_file_smi():
    with open("./examples/smiles.txt", "r") as f:
        data = f.readlines()  # read entire file as a single string
    
    lines = [line.rstrip("\n") for line in data] 
    return lines  # encode to bytes

@pytest.fixture
def sample_file_sdf():
    with open("./examples/short_structures.sdf", "r") as f:
        data = f.read()  # read entire file as one string

    # Split into molecule blocks at "$$$$"
    blocks = data.split("$$$$")

    # Strip empty blocks and whitespace
    blocks = [block.strip() for block in blocks if block.strip()]

    return blocks

def test_calculate_descriptors_chunk_smi(sample_file_smi: List[str]):
    

    # run the task synchronously
    result = calculate_descriptors_chunk_smi(sample_file_smi)

    # reading in the made json file
    proccessed_lines = []
    with open(f"./{result}","r") as f:
        proccessed_lines = json.load(f)
    
    assert "json" in result and "None" in result
    assert len(proccessed_lines) == 10

    for i,line in enumerate(proccessed_lines):
        assert len(proccessed_lines[i].keys()) == 219


def test_calculate_descriptors_chunk_sdf(sample_file_sdf: List[str]):

    # check if there are three blocks
    assert len(sample_file_sdf) == 7

    # run the task synchronously
    result = calculate_descriptors_chunk_sdf(sample_file_sdf)

    # reading in the made json file
    proccessed_json= []
    with open(f"./{result}","r") as f:
        proccessed_json = json.load(f)
    
    assert "json" in result and "None" in result
    assert len(proccessed_json) == 7

    for i,line in enumerate(proccessed_json):
        assert len(proccessed_json[i].keys()) == 219