# 1

- need to update db from sqlite to mysql / fallsprevention db so we can make sure that the cronb job secript that can be run properly 
    - currently in the github/actions folder the cronjob is setup to connect to the sqlite db but we need to change it to the mysql db
    - need to update the db connection string in the cronjob script to connect to the mysql db
    - in addition need to refactor the code to use the mysql db instead of the sqlite db in the main codebase