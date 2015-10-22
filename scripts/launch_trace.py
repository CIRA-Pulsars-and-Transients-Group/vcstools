#!/mnt/gleam/chen/pyws/bin/python
#
#    ICRAR - International Centre for Radio Astronomy Research
#    (c) UWA - The University of Western Australia, 2014
#    Copyright by UWA (in the framework of the ICRAR)
#    All rights reserved
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#    MA 02111-1307  USA
#

#******************************************************************************
#
#
# Who       When        What
# --------  ----------  -------------------------------------------------------
# cwu      06/May/2014  Created
#
import subprocess, os, sys, signal
import calendar, time
from threading import Timer

SET_TIMER = False # if running on the debug queue, please set it to true
EXPIRY_LENGTH = 3500 # one min to an hour (on Magnus debug queue)
wait_for_cpu_trace = True

def usage():
    print 'python launch_trace.py app'
    print 'e.g. python launch_trace.py ls -l'

def time_expired(sp):
    """
    If time is up, we need to save the CPU log
    """
    global wait_for_cpu_trace
    if (wait_for_cpu_trace):
        wait_for_cpu_trace = False # data race, who cares!?
        sp.send_signal(signal.SIGTERM)
        print "(Expired) Tracing return code:", sp.wait()

def trace(log_path=None, proc_path=None):
    #log_path = "/group/partner1024/cwu/chiles/perf_logs"
    if (log_path is None):
        log_path = os.path.dirname(os.path.realpath(__file__))
    #proc_path = "/home/cwu/chiles/processing"
    cmdlist = sys.argv[1:]
    start_time = calendar.timegm(time.gmtime())
    #print "CPU recording start_time: ", start_time
    cpu_logfile = '%s_cpu.log' % str(start_time)
    app_logfile = '%s_app.log' % str(start_time)
    h_app_logfile = open(log_path + '/' + app_logfile, 'w')

    sp = subprocess.Popen(cmdlist, stdout = h_app_logfile)

    if (proc_path is None):
        trace_io = False
    else:
        strace_path = proc_path + '/ioprofiler-trace.sh'
        if (os.path.exists(strace_path)):
            trace_io = True
            io_logfile = '%s_io.log' % str(start_time)
            cmd0 = "{0} -p {1} -o {2}/{3}".format(strace_path, sp.pid, log_path, io_logfile)
            print cmd0
            cmdlist0 = cmd0.split()
            sp0 = subprocess.Popen(cmdlist0)
        else:
            trace_io = False

    cmd1 = 'python /home/fkirsten/software/galaxy-scripts/scripts/trace_cpu_mem.py -o %s/%s -p %d' % (log_path, cpu_logfile, sp.pid)
    print cmd1
    cmdlist1 = cmd1.split()
    sp1 = subprocess.Popen(cmdlist1)

    if (SET_TIMER):
        t = Timer(EXPIRY_LENGTH, time_expired, [sp1])
        print "Timer set to {0} seconds".format(EXPIRY_LENGTH)
        t.start()

    print "Waiting...."
    print "Application return code:", sp.wait()
    global wait_for_cpu_trace
    if (wait_for_cpu_trace):
        wait_for_cpu_trace = False
        sp1.send_signal(signal.SIGTERM)
        print "Tracing return code:", sp1.wait()

    if (trace_io):
        sp0.send_signal(signal.SIGTERM) # terminate strace, which does not end
                                        # automatically

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        usage()
        sys.exit(1)
    trace()
