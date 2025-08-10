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
```

### docker environment

This command with go into the express backend and the
python backend to start up the docker containers
```
docker compose up --build
```