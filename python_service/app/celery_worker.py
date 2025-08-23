from celery import Celery
import descriptors.all_descriptors as all_ds  # your existing descriptor calculation code

celery_app = Celery(
    "tasks",
    broker="redis://localhost:6379/0",  # or your Redis server
    backend="redis://localhost:6379/1"
)

@celery_app.task
def calculate_descriptors_task(smiles_list):
    return all_ds.calc_all_descriptors(smiles_list)
