from celery import Celery
from app.config import settings

# Create Celery app
app = Celery('tasks', 
             broker=settings.CELERY_BROKER_URL, 
             backend=settings.CELERY_RESULT_BACKEND)

# Optional: Configure Celery
app.conf.update(
    task_track_started=True,
    task_time_limit=30 * 60,  # 30 minutes
)

# Import tasks to register them
from app import tasks

# This is important for Celery to discover tasks
if __name__ == '__main__':
    app.start()
