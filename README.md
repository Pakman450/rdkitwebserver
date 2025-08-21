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
