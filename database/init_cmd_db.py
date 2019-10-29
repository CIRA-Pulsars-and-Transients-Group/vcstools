import os
import sqlite3 as lite
try:
    DB_FILE = os.environ['CMD_VCS_DB_FILE']
except:
    print("environmental variable JOB_DB_FILE must be defined")
    
con = lite.connect(DB_FILE)
with con:
    cur = con.cursor()
    # ProcessVCS Table: stores process_vcs.py commands
    # Arguments - all arguments to the task
    # Obsid - Observation ID (gpsseconds)
    # Username - user who ran the job
    # Started - start time
    cur.execute("CREATE TABLE ProcessVCS(Rownum integer primary key autoincrement, Arguments TEXT, Obsid INT, UserId TEXT, Started date)")
with con:
    cur = con.cursor()
    # Commands Table: stores all commands of jobs that go to the queue
    # JobID - ID of job which ran the command
    # Command  - e.g. obsdownload, cotter, 
    # Arguments - all arguments to the task
    # Started - start time
    # Ended - end time NULL implies that the job is running or was recently terminated 
    # Exit - Command exit code (NULL and (Ended != NULL) means the job was terminated.
    cur.execute("CREATE TABLE Commands(Rownum integer primary key autoincrement, trackvcs INT, JobId INT, Command TEXT, Channels TEXT, Arguments TEXT, Started date, Ended date, Exit INT, FOREIGN KEY(trackvcs) REFERENCES ProcessVCS(Rownum))")

