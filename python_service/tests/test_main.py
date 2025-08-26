
import pytest
import io
import uuid
from pathlib import Path
import json

from main import app
from fastapi.testclient import TestClient

# client = TestClient(app)

@pytest.fixture
def client():
    # This fixture provides a TestClient for all tests
    with TestClient(app) as c:
        yield c

@pytest.fixture
def sample_file():
    with open("./examples/smiles.txt", "rb") as f: 
        contents = f.read()
        return io.BytesIO(contents)

@pytest.fixture
def size_dataset():
    with open("./examples/smiles.txt", "r") as f:
        data = f.readlines()  # read entire file as a single string
    
    lines = [line.rstrip("\n") for line in data] 
    return lines 

@pytest.fixture
def response_test(client, sample_file):
    response = client.post(
        "/v1/descriptors",
        files={"file": ("smiles.txt", sample_file, "text/plain")}
    )
    return response



def test_main(response_test, size_dataset):

    import main
    results = response_test.json()

    
    assert response_test.status_code == 200

    job_id_str = results["job_id"]
    job_id = uuid.UUID(job_id_str) 
    assert isinstance(job_id, uuid.UUID)

    # num of mols / chunk size = 10 / 5
    assert results['num_chunks'] == len(size_dataset) / main.CHUNK_SIZE 



def test_get_job_status(client):

    import main

    # create a fake merged file
    job_id = "12345678-1234-1234-1234-1234567890ab"
    merged_file = main.RESULTS_DIR / f"{job_id}_merged.json" 
    merged_file.write_text("{}")

    response = client.get(f"/v1/descriptors/{job_id}")
    
    assert response.status_code == 200
    
    data = response.json()

    assert data["status"] == "done"
    assert "merged_file" in data
    assert data["merged_file"] == f"{main.RESULTS_DIR}/{job_id}_merged.json"
