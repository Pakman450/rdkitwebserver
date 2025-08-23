from fastapi import FastAPI, UploadFile, File
from celery_worker import calculate_descriptors_task, celery_app
from celery.result import AsyncResult
import json
from pathlib import Path

app = FastAPI()

@app.post("/v1/descriptors")
async def submit_descriptors_job(file: UploadFile = File(...)):
    contents = await file.read()
    smiles_list = contents.decode("utf-8").splitlines()
    
    # Enqueue the task
    job = calculate_descriptors_task.delay(smiles_list)
    
    return {"job_id": job.id}

@app.get("/v1/descriptors/{job_id}")
def get_job_status(job_id: str):
    job = AsyncResult(job_id, app=celery_app)
    
    if job.state == "PENDING":
        return {"status": "pending"}
    elif job.state == "FAILURE":
        return {"status": "failed", "error": str(job.result)}
    elif job.state == "SUCCESS":
        # job.result is the file path
        file_path = Path(job.result)
        if file_path.exists():
            # Optional: return content or just path
            with open(file_path) as f:
                result = json.load(f)
            return {"status": "done", "result": result}
        else:
            return {"status": "done", "result": None, "warning": "file missing"}
    
    return {"status": job.state}
