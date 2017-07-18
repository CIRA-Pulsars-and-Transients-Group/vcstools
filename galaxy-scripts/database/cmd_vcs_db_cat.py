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
    
def add_database_function():
    batch_line ='function run\n' +\
                '{\n' +\
                '    # run command and add relevant data to the job database\n' +\
                '    # 1st parameter is command to be run (e.g. wsclean)\n' +\
                '    # 2nd parameter is parameters to that command (e.g. "-j $ncpus")\n' +\
                '    # 3rd parameter is datadir\n' +\
                '    # 4th parameter is project\n' +\
                '    # 5th parameter is obsid\n' +\
                '    # 6th parameter is chans [optional]\n' +\
                '    if [ -z "$6" ]; then\n' +\
                '        rownum=`cmd_start.py $1 -a "$2" -d $3 -p $4 -o $5`\n' +\
                '    else\n' +\
                '        rownum=`cmd_start.py $1 -a "$2" -d $3 -p $4 -o $5 -c $6`\n' +\
                '    fi\n' +\
                '    $1 $2\n' +\
                '    errcode=$?\n' +\
                '    cmd_stop.py $rownum -e $errcode\n' +\
                '    echo "cmd_stop.py $rownum -e $errcode"\n' +\
                '    if [ "$errcode" != "0" ]; then\n' +\
                '        exit $errcode\n' +\
                '    fi\n' +\
                '}\n'
    return batch_line

"""
batch_line = "run prepsubband '-ncpus 8 -lodm " + str(dm_start) +\
                " -dmstep " + str(dm_line[2]) + " -numdms 500 -numout " + str(numout) +\
                " -o " + str(obsid) + " /group/mwaops/vcs/" + str(obsid) + \
                "/pointings/" + str(pointing) + "/" + str(obsid) + "*.fits' " +\
                work_dir + ' blindsearch ' + obsid
"""

parser = OptionParser(usage = "usage: %prog <options>" +
"""
Print rows from the job database
""")
parser.add_option("-r", "--recent", dest="recent", metavar="HOURS", default=None, type=float, help="print only jobs started in the last N hours")
parser.add_option("-n", "--number", dest="n", metavar="N", default=20, type=int, help="number of jobs to print [default=%default]")
parser.add_option("-a", "--all", dest="all", action="store_true", help="print all lines of the database")
parser.add_option("-s", "--startrow", dest="startrow", default=0, type=int, help="ignore any row earlier than this one")
parser.add_option("-e", "--endrow", dest="endrow", default=None, type=int, help="ignore any row later than this one")
parser.add_option("-u", "--user", dest="user", default=None, type=str, help="Only prints one user's jobs.")
parser.add_option("-o", "--obsid", dest="obsid", default=None, type=str, help="Only prints one obsid's jobs.")
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

"""
print 'RowNum','JobId','TaskId','Datadir','Project','Obsid','Channels','Start','End','Exit','User'

for row in rows:
    print str(row['Rownum']).rjust(4),
    print str(row['JobId']).rjust(7),
    print str(row['TaskID']).rjust(4),
    print row['Datadir'] + os.sep,
    print row['Project'] + os.sep,
    print row['Obsid'],
    if row['Channels'] is not None:
        print row['Channels'].ljust(7),
    else:
        print '-------',
    print row['Started'][:19],
    if row['Ended'] is not None:
        print row['Ended'][11:19].ljust(8),
    else:
        print '--------',
    print row['Command'].ljust(15),
    if row['Exit'] is not None:
        print str(row['Exit']).rjust(3),
    else:
        print "---",
    print row['UserId'].ljust(10),
    if opts.args:
        print row['Arguments']
    print "\n"
"""
