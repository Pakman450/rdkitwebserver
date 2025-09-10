# app/celery_worker.py
import uuid
from pathlib import Path
import json
from celery import Celery
import lib.descriptors as all_ds
import os
import csv
from typing import List


#BROKER_URL = "redis://localhost:6379/0"

BROKER_URL ="pyamqp://guest:guest@localhost//",  # RabbitMQ broker
BACKEND_URL = "redis://localhost:6379/1" # redis backend

celery_app = Celery(
    "tasks",
    broker=BROKER_URL,
    backend=BACKEND_URL
)

RESULTS_DIR = Path("./results")
RESULTS_DIR.mkdir(exist_ok=True)

@celery_app.task
def calculate_descriptors_chunk_smi(smiles_chunk: List[str], job_id: str =None) -> str:
    """
    Calculate descriptors for a chunk of SMILES and save JSON.
    """
    result = all_ds.calc_all_descriptors_smi(smiles_chunk)
    chunk_file = RESULTS_DIR / f"{job_id}_{uuid.uuid4()}.json"
    with open(chunk_file, "w") as f:
        json.dump(result, f, indent=2)
    return str(chunk_file)

@celery_app.task
def calculate_descriptors_chunk_sdf(sdf_chunk: List[str], job_id: str =None) -> str:
    """
    Calculate descriptors for a chunk of SMILES and save JSON.
    """
    result = all_ds.calc_all_descriptors_sdf(sdf_chunk)
    chunk_file = RESULTS_DIR / f"{job_id}_{uuid.uuid4()}.json"
    with open(chunk_file, "w") as f:
        json.dump(result, f,indent=2)
    return str(chunk_file)


@celery_app.task
def merge_chunks(chunk_files: List[str], merged_file: str) -> str:
    """
    Merge multiple JSON chunk files into a single JSON with indentation,
    without loading all chunks into RAM. Also, the output also writes
    a csv file. s
    """

    basename, extension = os.path.splitext(merged_file)
    csv_f = open(f"{basename}.csv", "w+")



    with open(merged_file, "w") as out_f:
        out_f.write("[\n")
        first = True

        for fpath in chunk_files:
            with open(fpath) as f:
                data = json.load(f)

                if first:
                    csv_writer = csv.DictWriter(csv_f, data[0].keys())
                    csv_writer.writeheader()
                

                for item in data:
                    if not first:
                        out_f.write(",\n")
                    
                    csv_writer.writerow(item)
                    
                    # Pretty-print individual item with 2-space indent
                    item_str = json.dumps(item, indent=2)
                    # Add extra spaces to align with outer array
                    indented_item = "  " + item_str.replace("\n", "\n  ")
                    out_f.write(indented_item)
                    first = False

        out_f.write("\n]\n")

    return str(merged_file)