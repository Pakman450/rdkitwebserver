#!/bin/bash

celery -A lib.celery_worker  worker --loglevel=info