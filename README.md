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
cd python_service/app && uvicorn main:app --host 0.0.0.0 --port 5000 
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
cd python_server/app && uvicorn main:app --reload
```

Also, you have to keep celery on too
```
cd python_server/app && celery -A celery_worker.celery_app worker --loglevel=info
```

if you run:

```
cd python_server/app/examples && curl -X POST "http://localhost:8000/v1/descriptors" -F "file=@smiles.txt"
```

you will get a job id back. The job will run in the background, where
the celery worker will do the work. 

You can see the status of the job by doing:

```
curl "http://localhost:8000/v1/descriptors/JOBID"
```

where if completed, will return results

