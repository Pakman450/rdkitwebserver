# main.py
from fastapi import FastAPI, UploadFile, File
from lib.celery_worker import calculate_descriptors_chunk_smi, calculate_descriptors_chunk_sdf, merge_chunks, celery_app
from lib.format_checker import get_file_format

from celery import chord
from pathlib import Path
import uuid
import aiofiles

from pydantic import BaseModel

class JobResponse(BaseModel):
    job_id: str
    num_chunks: int


class JobIdResponse(BaseModel):
    status: str
    merged_file: str | None = None


async def stream_lines(file: UploadFile):
    while True:
        line_bytes = await file.read()
        if not line_bytes:
            break
        line = line_bytes.decode().strip()
        if not line:
            continue
        yield line   # yields one line at a time
        # once the caller consumes it, it's gone from memory


app = FastAPI(
    title="Descriptor API",
    description="Compute molecular descriptors in parallel using Celery.",
    version="1.0.0"
)

# NOTE: this is set to 12500 because
# current set up as 8 concurrency threads for workers 
# Also, set up for 100K SDF mols
CHUNK_SIZE = 444958
RESULTS_DIR = Path("./results")
RESULTS_DIR.mkdir(exist_ok=True)

@app.post("/v1/smi/descriptors", 
            summary="Submit descriptors computation job",
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
async def submit_smi_descriptors(file: UploadFile = File(...)) -> JobResponse:
    """
    Submit a SMILES file for descriptor calculation.

    - **file**: Upload a `.smi` or `.txt` file with one SMILES per line.
    - Returns a unique `job_id` and the number of chunks created.
    """

    # check file format
    print(get_file_format(file.filename))

    # read in smiles file
    contents = await file.read()
    smiles_list = contents.decode("utf-8").splitlines()

    # Create job id and chunkify smiles file
    job_id = str(uuid.uuid4())
    chunks = [smiles_list[i:i + CHUNK_SIZE] for i in range(0, len(smiles_list), CHUNK_SIZE)]

    # Create tasks for each chunk
    chunk_tasks = [calculate_descriptors_chunk_smi.s(chunk, job_id=job_id) for chunk in chunks]

    # Define path for final merged file
    merged_file = RESULTS_DIR / f"{job_id}_merged.json"

    # Use chord: merge_chunks runs after all chunk tasks finish
    # first argument of merge_chunks is injected automatically by chord(arg)
    job_result = chord(chunk_tasks)(merge_chunks.s(str(merged_file)))

    return JobResponse(job_id=job_id,num_chunks=len(chunks))

@app.post("/v1/smi/chunk/descriptors", 
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
    async for line in stream_lines(file):
    
        if not line:
            continue

        buffer.append(line)

        print(len(buffer))
        
        if len(buffer) >= CHUNK_SIZE:
            print("hehe")
            chunk_tasks.append(calculate_descriptors_chunk_smi.s(buffer, job_id=job_id))
            num_chunks += 1
            buffer = []

    if buffer:

        chunk_tasks.append(calculate_descriptors_chunk_smi.s(buffer, job_id=job_id))
        num_chunks += 1

    if chunk_tasks:
        chord(chunk_tasks)(merge_chunks.s(str(merged_file)))

    return JobResponse(job_id=job_id, num_chunks=num_chunks)


#TODO make this batchable for sdf file format
@app.post("/v1/sdf/descriptors",
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

    contents = await file.read()
    sdf_text = contents.decode("utf-8")

    # Split into molecule blocks
    sdf_list = [block.strip() for block in sdf_text.split("$$$$") if block.strip()]

    # Create job id and chunkify smiles file
    job_id = str(uuid.uuid4())
    chunks = [sdf_list[i:i + CHUNK_SIZE] for i in range(0, len(sdf_list), CHUNK_SIZE)]

    # Create tasks for each chunk
    chunk_tasks = [calculate_descriptors_chunk_sdf.s(chunk, job_id=job_id) for chunk in chunks]

    # Define path for final merged file
    merged_file = RESULTS_DIR / f"{job_id}_merged.json"

    # Use chord: merge_chunks runs after all chunk tasks finish
    # first argument of merge_chunks is injected automatically by chord(arg)
    job_result = chord(chunk_tasks)(merge_chunks.s(str(merged_file)))

    return JobResponse(job_id=job_id,num_chunks=len(chunks))

@app.get("/v1/job/{job_id}",
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



@app.get("/ping")
def ping():
    """
    Check if the service is alive.

    Returns:
        dict: `{"message": "pong"}`
    """
    return {"message": "pong"}