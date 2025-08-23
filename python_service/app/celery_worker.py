import uuid
from pathlib import Path
import json
from celery import Celery
import descriptors as all_ds

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
def calculate_descriptors_task(smiles_list):
    # Calculate descriptors in memory
    result = all_ds.calc_all_descriptors(smiles_list)
    
    # Write result to disk
    out_file = RESULTS_DIR / f"{uuid.uuid4()}.json"
    with open(out_file, "w") as f:
        json.dump(result, f)
    
    # Return file path as result
    return str(out_file)
