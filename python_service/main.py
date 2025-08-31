# main.py



from lib.routers import submit,job_id

from fastapi import FastAPI




app = FastAPI(
    title="Descriptor API",
    description="Compute molecular descriptors in parallel using Celery.",
    version="1.0.0"
)

# Include each router separately
app.include_router(submit.router, prefix="/v1")
app.include_router(job_id.router, prefix="/v1/job")

@app.get("/ping")
def ping():
    """
    Check if the service is alive.

    Returns:
        dict: `{"message": "pong"}`
    """
    return {"message": "pong"}