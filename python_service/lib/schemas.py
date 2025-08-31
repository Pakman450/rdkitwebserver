from pydantic import BaseModel

class JobResponse(BaseModel):
    job_id: str
    num_chunks: int


class JobIdResponse(BaseModel):
    status: str
    merged_file: str | None = None