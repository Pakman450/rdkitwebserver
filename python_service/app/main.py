from fastapi import FastAPI, UploadFile, File
from celery_worker import calculate_descriptors_task, celery_app
from celery.result import AsyncResult

app = FastAPI()

# Submit a file for descriptor calculation
@app.post("/v1/descriptors")
async def submit_descriptors_job(file: UploadFile = File(...)):
    contents = await file.read()
    smiles_list = contents.decode("utf-8").splitlines()
    
    # Enqueue the job
    job = calculate_descriptors_task.delay(smiles_list)
    
    return {"job_id": job.id}

# Check the status of a job
@app.get("/v1/descriptors/{job_id}")
def get_job_status(job_id: str):
    job = AsyncResult(job_id, app=celery_app)
    if job.state == "PENDING":
        return {"status": "pending"}
    elif job.state == "SUCCESS":
        return {"status": "done", "result": job.result}
    elif job.state == "FAILURE":
        return {"status": "failed", "error": str(job.result)}
    return {"status": job.state}
