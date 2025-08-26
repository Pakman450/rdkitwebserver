
import pytest
from main import app
from fastapi.testclient import TestClient

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



def test_task_directly(sample_file):
    
    l_smiles = []
    for mol in sample_file:
        l_smiles.append(mol)

    # run the task synchronously
    result = calculate_descriptors_chunk(l_smiles)
    
    print(result)

    assert "descriptors" in result


# def test_descriptors_main(sample_file):
#     response = client.post(
#         "/v1/descriptors",
#         files={"file": ("test.csv", sample_file, "text/csv")}
#     )
#     data = response.json()
#     job_id = data["job_id"]

#     # Get the response and the keys of the json response
#     assert list(data.keys()) == ['job_id', "num_chunks"]
#     assert response.status_code == 200
    

#     # Get task result (runs synchronously in eager mode)
#     task = AsyncResult(job_id, app=celery_app)
#     assert task.ready()
#     result = task.result
#     assert "descriptors" in result
