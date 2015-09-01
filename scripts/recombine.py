#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import subprocess
import os
import sys
import glob
import getopt
from mpi4py import MPI

def usage(opts={}):
    print "recombine.py -s <start time> -o <obsid> -w <working directory> -e <recombine executable>\n"

testsize = 327680000


if __name__ == '__main__':


    the_options = {'recombine': "recombine", 'start': int(0), 'root' : "./", 'obsid' : int(0), 'testit' : 0, 'skip' : " ", 'read_pfb' : 1}
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hce:o:s:tw:")
    except getopt.GetoptError:
        usage(the_options)
        sys.exit(-1)
    finally:
        if len(sys.argv) < 1:
            usage(the_options)
            sys.exit(-1)
        
    for opt,arg in opts:    
    
        if (opt == "-h"):
            usage(the_options)
        elif (opt == "-e"):
            the_options['recombine'] = arg
        elif (opt == "-s"):
            the_options['start'] = int(arg)
        elif (opt == "-t"):
            the_options['testit'] = 1
        elif (opt == "-o"):
            the_options['obsid'] = int(arg)
        elif (opt == "-w"):
            the_options['root'] = arg
        elif (opt == "-c"):
            the_options['skip'] = "-c"
            the_options['read_pfb'] = 0

    if (the_options['start'] == 0):
        usage(the_options)
        sys.exit(-1)
    if (the_options['obsid'] == 0):
        usage(the_options)
        sys.exit(-1)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    time_to_combine = the_options['start']+rank
    files_glob = "%s/combined/%d_%d_ch*.dat" % (the_options['root'],the_options['obsid'],time_to_combine)
    broken = 24;
    for to_check in sorted(glob.glob(files_glob)):
        file_statinfo = os.stat(to_check)
        if (file_statinfo.st_size == testsize):
            broken = broken-1
        else:    
            os.remove(to_check)

    if (broken > 0):
        f=[]
        for vcs in range(1,17):
            for stream in [0,1]:
                file_to_combine = "%d_%d_vcs%02d_%1d.dat " % (the_options['obsid'],time_to_combine,vcs,stream)
                f.append(file_to_combine)


        recombine_line = "%s %s -o %s -t %d -m %s/%d.metafits -i %s/combined -f " % (the_options['recombine'],the_options['skip'],the_options['obsid'],time_to_combine,the_options['root'],the_options['obsid'],the_options['root'])



        for f_to_r in f:
            recombine_line = "{0} {1}/raw/{2}".format(recombine_line, the_options['root'],f_to_r)

        recombine_line = "%s\n" % recombine_line
        log_name="{0}/recombine_{1}.log".format(working_dir,time_to_combine)
        with open(log_name, 'w') as log:
            subprocess.call(recombine_line,shell=True,stdout=log,stderr=log)
        
    comm.Barrier()
       
   



