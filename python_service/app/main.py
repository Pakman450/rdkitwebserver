# app/main.py
from fastapi import FastAPI, UploadFile, File
from celery_worker import calculate_descriptors_chunk, merge_chunks, celery_app
from celery import chord
from pathlib import Path
import uuid

app = FastAPI()

CHUNK_SIZE = 5
RESULTS_DIR = Path("./results")
RESULTS_DIR.mkdir(exist_ok=True)

@app.post("/v1/descriptors")
async def submit_descriptors_job(file: UploadFile = File(...)):

    #read in smiles file
    contents = await file.read()
    smiles_list = contents.decode("utf-8").splitlines()

    # Create job id and chunkify smiles file
    job_id = str(uuid.uuid4())
    chunks = [smiles_list[i:i + CHUNK_SIZE] for i in range(0, len(smiles_list), CHUNK_SIZE)]

    # Create tasks for each chunk
    chunk_tasks = [calculate_descriptors_chunk.s(chunk, job_id=job_id) for chunk in chunks]

    # Define path for final merged file
    merged_file = RESULTS_DIR / f"{job_id}_merged.json"

    # Use chord: merge_chunks runs after all chunk tasks finish
    # first argument of merge_chunks is injected automatically by chord(arg)
    job_result = chord(chunk_tasks)(merge_chunks.s(str(merged_file)))

    return {"job_id": job_id, "num_chunks": len(chunks)}


@app.get("/v1/descriptors/{job_id}")
def get_job_status(job_id: str):
    merged_file = RESULTS_DIR / f"{job_id}_merged.json"
    
    # Check if merge task result exists in Celery backend
    # Find the chord result if needed
    if merged_file.exists():
        return {"status": "done", "merged_file": str(merged_file)}
    
    # Otherwise, check pending chunk tasks
    # Celery stores task states in backend
    # This is optional; for simplicity, we just check if merged file exists
    return {"status": "pending"}
