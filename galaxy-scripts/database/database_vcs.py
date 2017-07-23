#!/usr/bin/env python
import os, datetime, logging
import sqlite3 as lite
from optparse import OptionParser #NB zeus does not have argparse!
from dateutil.parser import parse
from dateutil.relativedelta import relativedelta

DB_FILE = os.environ['CMD_VCS_DB_FILE']

def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d
    
def database_command(args, obsid):
        DB_FILE = os.environ['CMD_VCS_DB_FILE']
        args_string = ""
        for a in args:
                if not a == args[0]:
                        args_string = args_string + str(a) + " "
                        
        con = lite.connect(DB_FILE)
        con.isolation_level = 'EXCLUSIVE'
        con.execute('BEGIN EXCLUSIVE')
        with con:
                cur = con.cursor()
                
                cur.execute("INSERT INTO ProcessVCS(Arguments, Obsid, UserId, Started) VALUES(?, ?, ?, ?)", (args_string, obsid, os.environ['USER'], datetime.datetime.now()))
                vcs_command_id = cur.lastrowid
        return vcs_command_id
    
def add_database_function():
    batch_line ='function run\n' +\
                '{\n' +\
                '    # run command and add relevant data to the job database\n' +\
                '    # 1st parameter is command to be run including aprun (e.g. wsclean)\n' +\
                '    # 2nd parameter is parameters to that command (e.g. "-j $ncpus")\n' +\
                '    # 3rd parameter is vcs row number\n' +\
                '    command="${1##*/}"\n' +\
                '    command="${command##* }"\n' +\
                '    rownum=`database_vcs.py -m "s" -v "$3" -c $command -a "$2"\n' +\
                '    $1 $2\n' +\
                '    errcode=$?\n' +\
                '    database_vcs.py -m "e" -r $rownum --errorcode $errcode\n' +\
                '    echo "database_vcs.py -m "e" -r $rownum --errorcode $errcode"\n' +\
                '    if [ "$errcode" != "0" ]; then\n' +\
                '        exit $errcode\n' +\
                '    fi\n' +\
                '}\n'
    return batch_line
    
    
def database_script_start(vcs_id, command, arguments):
    import datetime
    con = lite.connect(DB_FILE)
    with con:
        cur = con.cursor()
        #Commands(Rownum integer primary key autoincrement, trackvcs INT, JobId INT, Command TEXT, Channels TEXT, Arguments TEXT, Started date, Ended date, Exit INT, FOREIGN KEY(trackvcs) REFERENCES ProcessVCS(Rownum))
        cur.execute("INSERT INTO Commands (trackvcs, Command, Arguments, Started) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (vcs_id, command, arguments, datetime.datetime.now()))
        row_id = cur.lastrowid
    return row_id

def database_script_stop(rownum, errorcode):    
    end_time = datetime.datetime.now()

    con = lite.connect(DB_FILE)
    with con:
        cur = con.cursor()
        cur.execute("SELECT Ended FROM Commands WHERE Rownum=?", (rownum,))
        ended = cur.fetchone()[0]
        if ended is not None:
            logging.warn("Overwriting existing completion time: %s" % ended)
        cur.execute("UPDATE Commands SET Ended=?, Exit=? WHERE Rownum=?", (end_time, errorcode, rownum))



"""
batch_line = "run prepsubband '-ncpus 8 -lodm " + str(dm_start) +\
                " -dmstep " + str(dm_line[2]) + " -numdms 500 -numout " + str(numout) +\
                " -o " + str(obsid) + " /group/mwaops/vcs/" + str(obsid) + \
                "/pointings/" + str(pointing) + "/" + str(obsid) + "*.fits' " +\
                work_dir + ' blindsearch ' + obsid
"""

if __name__ == '__main__':

    parser = OptionParser(usage = "usage: %prog <options>" +
    """
    Script used to manage the VCS database by recording the scripts process_vcs.py uses and prints the databse
    """)
    parser.add_option("-m", "--mode", dest="mode", metavar="mode", default='v', type=str, help='This script has three modes: "v" used to view the database, "s" used to start a record of a script on the database and "e" used to record the end time and error code of a script on the database. Default mode is v')
    
    view_options = OptionGroup(parser, 'View Options')
    view_options.add_option("--recent", dest="recent", metavar="HOURS", default=None, type=float, help="print only jobs started in the last N hours")
    view_options.add_option("-n", "--number", dest="n", metavar="N", default=20, type=int, help="number of jobs to print [default=%default]")
    view_options.add_option("--all", dest="all", action="store_true", help="print all lines of the database")
    view_options.add_option("-s", "--startrow", dest="startrow", default=0, type=int, help="ignore any row earlier than this one")
    view_options.add_option("-e", "--endrow", dest="endrow", default=None, type=int, help="ignore any row later than this one")
    view_options.add_option("-u", "--user", dest="user", default=None, type=str, help="Only prints one user's jobs.")
    view_options.add_option("-o", "--obsid", dest="obsid", default=None, type=str, help="Only prints one obsid's jobs.")
    
    start_options = OptionGroup(parser, 'Script Start Options')
    start_options.add_option("-v", "--vcs_id", dest="vcs_id", default=None, type=str, help="The row number of the process vcs command of the databse")
    start_options.add_option("-c", "--command", dest="command", default=None, type=str, help="The script name being run. eg volt_download.py.")
    start_options.add_option("-a", "--argument", dest="argument", default=None, type=str, help="The arguments that script used.")
    
    end_options = OptionGroup(parser, 'Script End Options')
    end_options.add_option("--errorcode", dest="errorcode", default=None, type=str, help="Error code of scripts.")
    end_options.add_option("-r", "--rownum", dest="rownum", default=None, type=str, help="The row number of the script.")
    opts, args = parser.parse_args()
    
    
    if args.mode == "s":
        if len(args) != 4:
            parser.error("incorrect number of arguments")
        database_script_start(args.vcs_id, args.command, args.arguments)
    elif args.mode == "e":
        if len(args) != 3:
            parser.error("incorrect number of arguments")
        database_script_stop(args.rownum, args.errorcode)
    else:
        con = lite.connect(DB_FILE)
        #con = lite.connect(DB_FILE, detect_types=lite.PARSE_DECLTYPES|lite.PARSE_COLNAMES) # return datetime as datetime objects
        con.row_factory = dict_factory

        if len(args) != 0:
            parser.error("Incorrect number of arguments")

        query = "SELECT * FROM ProcessVCS"
            
        print opts.user
        if opts.user:
            query += " WHERE UserId='" + str(opts.user) + "'"

        if opts.obsid:
            query += " WHERE Arguments LIKE '%" + str(opts.obsid) + "%'"

        if opts.recent is not None:
            query += ''' WHERE Started > "%s"''' % str(datetime.datetime.now() - relativedelta(hours=opts.recent))
            logging.debug(query)

        with con:
            cur = con.cursor()
            print query
            cur.execute(query)
            rows = cur.fetchall()

        if opts.startrow or opts.endrow:
            rows = rows[opts.startrow:]
            if opts.endrow is not None:
                rows = rows[:opts.endrow+1]
        elif not (opts.all or opts.recent):
            rows = rows[-opts.n:]

        print 'RowNum','Obsid','Started','UserID','Arguments'

        for row in rows:
            print str(row['Rownum']).rjust(4),
            print row['Obsid'],
            print row['Started'][:19],
            print row['UserId'].ljust(10),
            print row['Arguments']
            print "\n"

