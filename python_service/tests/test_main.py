import pytest
import uuid
from pathlib import Path
from io import BytesIO
from fastapi.testclient import TestClient

from main import app, RESULTS_DIR, CHUNK_SIZE

EXAMPLES_FILE = Path("./examples/smiles.txt")

@pytest.fixture(scope="session")
def client():
    with TestClient(app) as c:
        yield c

@pytest.fixture
def file_request_data(client):
    """
    Reads the smiles.txt file and posts it to the /v1/descriptors endpoint.
    Returns both the response and the text lines for assertions.
    """
    lines = EXAMPLES_FILE.read_text().splitlines()
    binary_file = BytesIO(EXAMPLES_FILE.read_bytes())

    response = client.post(
        "/v1/descriptors",
        files={"file": ("smiles.txt", binary_file, "text/plain")}
    )

    return {"response": response, "lines": lines}

def test_main(file_request_data):
    response = file_request_data["response"]
    lines = file_request_data["lines"]

    assert response.status_code == 200

    results = response.json()
    job_id = uuid.UUID(results["job_id"])
    assert isinstance(job_id, uuid.UUID)

    assert results["num_chunks"] == len(lines) / CHUNK_SIZE

def test_get_job_status_done(client):
    job_id = "12345678-1234-1234-1234-1234567890ab"
    merged_file = RESULTS_DIR / f"{job_id}_merged.json"
    merged_file.write_text("{}")

    response = client.get(f"/v1/descriptors/{job_id}")
    assert response.status_code == 200

    data = response.json()
    assert data["status"] == "done"
    assert data["merged_file"] == str(merged_file)

def test_get_job_status_pending(client):
    job_id = "PENDING"

    response = client.get(f"/v1/descriptors/{job_id}")
    merged_file = RESULTS_DIR / f"{job_id}_merged.json"
    assert response.status_code == 200

    data = response.json()
    assert data["status"] == "pending"
    assert "merged_file" not in data

