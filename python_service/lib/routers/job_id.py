
from fastapi import APIRouter, UploadFile, File
from pathlib import Path

from lib.schemas import JobIdResponse


router = APIRouter()

RESULTS_DIR = Path("./results")
RESULTS_DIR.mkdir(exist_ok=True)

@router.get("/{job_id}",
            summary="get status for a job id",
            description="""
                Submit a job id that you have received from job submission, 
                and get a status for that submitted job.
            """,
            response_model=JobIdResponse,
            responses={
                200: {
                    "description": "Job successfully submitted",
                    "content": {
                        "application/json": {
                            "example": {
                                "status": "done",
                                "merged_file": "123e4567-e89b-12d3-a456-426614174000_merged.json"
                            }
                        }
                    },
                }
            },
        )
def get_job_status(job_id: str)->JobIdResponse:
    merged_file = RESULTS_DIR / f"{job_id}_merged.json"
    
    # Check if merge task result exists in Celery backend
    # Find the chord result if needed
    if merged_file.exists():
        return JobIdResponse(status="done", merged_file=str(merged_file))
    
    # Otherwise, check pending chunk tasks
    # Celery stores task states in backend
    # This is optional; for simplicity, we just check if merged file exists
    return JobIdResponse(status="pending", merged_file=None)