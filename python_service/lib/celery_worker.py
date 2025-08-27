# app/celery_worker.py
import uuid
from pathlib import Path
import json
from celery import Celery
import lib.descriptors as all_ds

BROKER_URL = "redis://localhost:6379/0"
BACKEND_URL = "redis://localhost:6379/1"

celery_app = Celery(
    "tasks",
    broker=BROKER_URL,
    backend=BACKEND_URL
)

RESULTS_DIR = Path("./results")
RESULTS_DIR.mkdir(exist_ok=True)

@celery_app.task
def calculate_descriptors_chunk_smi(smiles_chunk, job_id=None):
    """
    Calculate descriptors for a chunk of SMILES and save JSON.
    """
    result = all_ds.calc_all_descriptors_smi(smiles_chunk)
    chunk_file = RESULTS_DIR / f"{job_id}_{uuid.uuid4()}.json"
    with open(chunk_file, "w") as f:
        json.dump(result, f)
    return str(chunk_file)

@celery_app.task
def calculate_descriptors_chunk_sdf(sdf_chunk, job_id=None):
    """
    Calculate descriptors for a chunk of SMILES and save JSON.
    """
    result = all_ds.calc_all_descriptors_sdf(sdf_chunk)
    chunk_file = RESULTS_DIR / f"{job_id}_{uuid.uuid4()}.json"
    with open(chunk_file, "w") as f:
        json.dump(result, f,indent=2)
    return str(chunk_file)



@celery_app.task
def merge_chunks(chunk_files, merged_file):
    """
    Merge multiple JSON chunk files into a single JSON.
    """
    combined = []
    for fpath in chunk_files:
        with open(fpath) as f:
            combined.extend(json.load(f))
    with open(merged_file, "w") as f:
        json.dump(combined, f,indent=2)
    return str(merged_file)
