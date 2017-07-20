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
                '    command="${1##*/}"
                '    rownum=`cmd_start.py $command -a "$2"\n' +\
                '    $1 $2\n' +\
                '    errcode=$?\n' +\
                '    cmd_stop.py $rownum -e $errcode\n' +\
                '    echo "cmd_stop.py $rownum -e $errcode"\n' +\
                '    if [ "$errcode" != "0" ]; then\n' +\
                '        exit $errcode\n' +\
                '    fi\n' +\
                '}\n'
    return batch_line
    
    
def database_script_start():
    con = lite.connect(DB_FILE)
    with con:
        cur = con.cursor()
        cur.execute("INSERT INTO Commands (Datadir, Project, Obsid, JobId, TaskId, Command, Channels, Arguments, UserId, Started) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (opts.datadir, opts.project, opts.obsid, JOB_ID, TASK_ID, args[0], opts.chans, opts.arguments, os.environ['USER'], datetime.datetime.now()))
        row_id = cur.lastrowid
    return row_id

def database_script_stop():
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    if opts.time and opts.errfile:
        logging.warn("--time and --errfile both set. --time will be used")

    if opts.time:
        end_time = parse(opts.time)
    elif opts.errfile:
        opts.errfile = os.path.expanduser(os.path.expandvars(opts.errfile))
        end_time = datetime.datetime.fromtimestamp(os.path.getmtime(opts.errfile))
    else:
        end_time = datetime.datetime.now()

    con = lite.connect(DB_FILE)
    with con:
        cur = con.cursor()
        if not opts.force:
            cur.execute("SELECT Ended FROM Commands WHERE Rownum=?", (args[0],))
            ended = cur.fetchone()[0]
            if ended is not None:
                if opts.no_overwrite:
                    raise RuntimeError, "Ended is already set"
                else:
                    logging.warn("Overwriting existing completion time: %s" % ended)
        cur.execute("UPDATE Commands SET Ended=?, Exit=? WHERE Rownum=?", (end_time, opts.exit, args[0]))



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
    view_options = OptionGroup(parser, 'View Options')
    view_options.add_option("-r", "--recent", dest="recent", metavar="HOURS", default=None, type=float, help="print only jobs started in the last N hours")
    view_options.add_option("-n", "--number", dest="n", metavar="N", default=20, type=int, help="number of jobs to print [default=%default]")
    view_options.add_option("-a", "--all", dest="all", action="store_true", help="print all lines of the database")
    view_options.add_option("-s", "--startrow", dest="startrow", default=0, type=int, help="ignore any row earlier than this one")
    view_options.add_option("-e", "--endrow", dest="endrow", default=None, type=int, help="ignore any row later than this one")
    view_options.add_option("-u", "--user", dest="user", default=None, type=str, help="Only prints one user's jobs.")
    view_options.add_option("-o", "--obsid", dest="obsid", default=None, type=str, help="Only prints one obsid's jobs.")
    start_options = OptionGroup(parser, 'Script Start Options')
    end_options = OptionGroup(parser, 'Script End Options')
    opts, args = parser.parse_args()

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

