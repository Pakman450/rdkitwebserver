from typing import AsyncGenerator
from fastapi import  UploadFile

async def stream_sdf_lines(file: UploadFile, chunk_size: int = 64 * 1024) -> AsyncGenerator[str, None]:
    """
    Stream an uploaded text file line-by-line without loading it all into RAM.
    Reads fixed-size byte chunks, splits on newlines, and yields decoded lines.
    """
    pending = b""
    while True:
        chunk = await file.read(chunk_size)  # returns b"" at EOF
        if not chunk:
            break
        pending += chunk
        # split into full lines; keep the trailing partial line in `pending`
        *full_lines, pending = pending.split(b"$$$$")
        for raw in full_lines:
            mol = raw.strip()  # removes \r too
            if mol:
                yield mol.decode("utf-8", errors="ignore")

    # flush any leftover bytes as the last line
    if pending:
        mol = pending.strip()
        if mol:
            yield mol.decode("utf-8", errors="ignore")
