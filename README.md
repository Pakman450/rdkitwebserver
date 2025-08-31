# rdkitwebserver


## To start development

We need multiple environments to set up first

### Make sure you have python packages

Make sure you are not in a conda environment.


These commands will read in `requirements.txt` to download python packages 
for python functionality
```
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### docker environment

This command with go into the express backend and the
python backend to start up the docker containers
```
docker compose up --build
```

### to start development mode
for python service
``` 
cd python_service && uvicorn main:app --host 0.0.0.0 --port 5000 
```

for celery
```
cd python_service && celery -A lib.celery_worker worker --loglevel=info
```

for express
``` 
cd backend && nodemon index.js
```



## python backend
### to develop
First, you need to run redis-server
```
sudo systemctl start redis-server
```

You can see the status of the redis-server
```
sudo systemctl status redis-server
```

Then, start the fastapi app by using uvicorn
```
cd python_server && uvicorn main:app --reload
```

Also, you have to keep celery on too
```
cd python_service && celery -A lib.celery_worker worker --loglevel=info
```

if you run:

```
cd python_server/app/examples && curl -X POST "http://localhost:8000/v1/smi/descriptors" -F "file=@smiles.txt"
```

you will get a job id back. The job will run in the background, where
the celery worker will do the work. 

You can see the status of the job by doing:

```
curl "http://localhost:8000/v1/job/JOBID"
```

where if completed, will return results

### Ensure these criteria before making changes
-[] add test cases

-[] submit benchmarking



# Benchmarking
## no aiofiles
### 444957 smiles divided up by 8 threads
time: 1358 sec (55620 chunk size)
### 444957 smiles all computed under 1 thread
time: 3162 sec (444958 chunk size) 

### 444957 smiles divided up by 8 threads line by line
time: 1652 sec (55620 chunk size)
### 444957 smiles all computed under 1 thread line by line
time: 3247 sec (444958 chunk size) 
