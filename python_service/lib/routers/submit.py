from fastapi import APIRouter, UploadFile, File
from pathlib import Path
from celery import chord
import uuid

from lib.schemas import JobResponse
from lib.stream_smi import stream_smi_lines
from lib.stream_sdf import stream_sdf_lines
from lib.celery_worker import calculate_descriptors_chunk_smi, calculate_descriptors_chunk_sdf, merge_chunks, celery_app
router = APIRouter()

CHUNK_SIZE = 5
RESULTS_DIR = Path("./results")
RESULTS_DIR.mkdir(exist_ok=True)

@router.post("/smi/chunk/descriptors", 
            summary="Submit descriptors computation job via chunks",
            description="""
            Upload .smi/.txt file containing SMILES strings.
            The service splits molecules into chunks, runs descriptor calculations in parallel,  
            and returns a `job_id` you can use to check status.
            """,
            response_model=JobResponse,
            responses={
                200: {
                    "description": "Job successfully submitted",
                    "content": {
                        "application/json": {
                            "example": {
                                "job_id": "123e4567-e89b-12d3-a456-426614174000",
                                "num_chunks": 3
                            }
                        }
                    },
                }
            },
            )
async def submit_smi_chunk_descriptors(file: UploadFile = File(...)) -> JobResponse:
    job_id = str(uuid.uuid4())
    merged_file = RESULTS_DIR / f"{job_id}_merged.json"

    chunk_tasks = []
    buffer = []
    num_chunks = 0

    # Read uploaded file line by line
    async for line in stream_smi_lines(file):
    
        if not line:
            continue

        buffer.append(line)
        
        if len(buffer) >= CHUNK_SIZE:

            chunk_tasks.append(calculate_descriptors_chunk_smi.s(buffer, job_id=job_id))
            num_chunks += 1
            buffer = []

    if buffer:

        chunk_tasks.append(calculate_descriptors_chunk_smi.s(buffer, job_id=job_id))
        num_chunks += 1

    if chunk_tasks:
        chord(chunk_tasks)(merge_chunks.s(str(merged_file)))

    return JobResponse(job_id=job_id, num_chunks=num_chunks)


@router.post("/sdf/descriptors",
            summary="Submit descriptors computation job",
            description="""
            Upload .sdf file containing sdf mol blocks.
            The service splits molecules into chunks, runs descriptor calculations in parallel,  
            and returns a `job_id` you can use to check status.
            """,
            response_model=JobResponse,
            responses={
                200: {
                    "description": "Job successfully submitted",
                    "content": {
                        "application/json": {
                            "example": {
                                "job_id": "123e4567-e89b-12d3-a456-426614174000",
                                "num_chunks": 3
                            }
                        }
                    },
                }
            },
            )
async def submit_sdf_descriptors(file: UploadFile = File(...)) -> JobResponse:


    job_id = str(uuid.uuid4())
    merged_file = RESULTS_DIR / f"{job_id}_merged.json"

    chunk_tasks = []
    buffer = []
    num_chunks = 0

    # Read uploaded file line by line
    async for sdf_mol in stream_sdf_lines(file):
    
        if not sdf_mol:
            continue

        buffer.append(sdf_mol)
        
        if len(buffer) >= CHUNK_SIZE:

            chunk_tasks.append(calculate_descriptors_chunk_sdf.s(buffer, job_id=job_id))
            num_chunks += 1
            buffer = []

    if buffer:

        chunk_tasks.append(calculate_descriptors_chunk_sdf.s(buffer, job_id=job_id))
        num_chunks += 1

    if chunk_tasks:
        chord(chunk_tasks)(merge_chunks.s(str(merged_file)))

    return JobResponse(job_id=job_id, num_chunks=num_chunks)