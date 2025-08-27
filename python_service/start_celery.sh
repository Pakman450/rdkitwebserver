#!/bin/bash

# To purge any tasks left in the queue
celery -A lib.celery_worker purge -f

celery -A lib.celery_worker  worker --loglevel=info
