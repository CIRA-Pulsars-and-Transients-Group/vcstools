#!/usr/bin/env python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import subprocess
import os
import sys
import glob
import getopt
from mpi4py import MPI
import logging
import argparse

logger = logging.getLogger(__name__)



def usage(opts={}):
    logger.info("recombine.py -s <start time> -o <obsid> -w <working directory> -e <recombine executable>\n")

testsize = 327680000


if __name__ == '__main__':

    from sys import argv

    the_options = {'recombine': "recombine", 'start': int(0), 'root' : "./", 'obsid' : int(0), 'testit' : 0, 'skip' : " ", 'read_pfb' : 1}
  
    loglevels = dict(DEBUG=logging.DEBUG,
                    INFO=logging.INFO,
                    WARNING=logging.WARNING,
                    ERROR = logging.ERROR)    
 
    
    #Arguments
    parser = argparse.ArgumentParser(description="""Combines raw VCS data""")
    parser.add_argument("-o", "--obsid", type=int, help="Input Observation ID")
    parser.add_argument("-s", "--start", type=int, help="GPS time to test")
    parser.add_argument("-w", "--data_dir", type=str, help="Directory containing the raw data")
    parser.add_argument("-e", "--recombine", type=str, help="filename of the recombine function") 

    args = parser.parse_args()
    
    #logger.setLevel(loglevels[args.loglvl])
    #logger.setLevel(DEBUG)
    #ch = logging.StreamHandler()
    #ch.setLevel(loglevels[args.loglvl])
    #ch.setLevel(DEBUG)
    #formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    #ch.setFormatter(formatter)
    #logger.addHandler(ch)


    """
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hce:o:s:tw:")
        logger.debug("options: {0}".format(opts))
        logger.debug("arguments: {0}".format(args))
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
    """
    #Setting the_options to the arg inputs
    the_options["obsid"] = args.obsid
    the_options["start"] = args.start
    the_options["root"] = args.data_dir
    the_options["recombine"] = args.recombine
    
    logger.info("Input Options:")
    usage(the_options)

    if (the_options['start'] == 0):
        usage(the_options)
        sys.exit(-1)
    if (the_options['obsid'] == 0):
        usage(the_options)
        sys.exit(-1)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    time_to_combine = the_options['start']+rank
    #files_glob = "{0}/combined/{1}_{2}_ch*.dat" % (the_options['root'],the_options['obsid'],time_to_combine)
    #Need to change the 
    files_glob = "{0}/recombine_test_delete_this/{1}_{2}_ch*.dat".format(the_options['root'], the_options['obsid'], time_to_combine)
 
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
                file_to_combine = "{0}_{1}_vcs{2}02_{3}1.dat ".format(the_options['obsid'],time_to_combine,vcs,stream)
                f.append(file_to_combine)


        recombine_line = "{0} {1} -o {2} -t {3} -m {4}/{5}_metafits_ppds.fits -i {6}/recombine_test_delete_this -f ".format(the_options['recombine'],the_options['skip'],the_options['obsid'],time_to_combine,the_options['root'],the_options['obsid'],the_options['root'])



        for f_to_r in f:
            recombine_line = "{0} {1}/raw/{2}".format(recombine_line, the_options['root'],f_to_r)

        recombine_line = "{0}\n".format(recombine_line)
        #print "this is what is called: " + recombine_line
        #subprocess.call("which recombine",shell=True)
        subprocess.call(recombine_line,shell=True)

        #log_name="{0}/recombine_{1}.log".format(the_options['root'],time_to_combine)

        #with open(log_name, 'w') as log:
        #    subprocess.call(recombine_line,shell=True,stdout=log,stderr=log)

    comm.Barrier()
       
   



