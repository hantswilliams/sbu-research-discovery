# sbu-research-discovery

## getting it started from scratch...
- first run the app for the first time locally to create the database and tables
    - `python app.py`
- then run the following command to populate the database with the data from the csv file
    - `python populate_db.py`
- then can run the cron job to update the database with the latest data
    - `python cron_job.py`