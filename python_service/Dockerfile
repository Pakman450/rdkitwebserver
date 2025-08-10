# Start from official Python image
FROM python:3.10-slim

# Install system dependencies needed for RDKit
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip

# Set working directory
WORKDIR /app

# Copy requirements and install
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy service code
COPY . .

# Expose FastAPI port
EXPOSE 5000

# Set working directory
WORKDIR /app/app

# Start FastAPI server
CMD ["uvicorn", "main:app", "--host" , "0.0.0.0", "--port", "5000"]
