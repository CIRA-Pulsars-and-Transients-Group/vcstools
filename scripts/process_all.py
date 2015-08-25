#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import subprocess
import os
import sys
import glob
import time
import getopt

import urllib
import urllib2
import json

import datetime

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

# Append the service name to this base URL, eg 'con', 'obs', etc.
BASEURL = 'http://ngas01.ivec.org/metadata/'

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

def sfreq(freqs):

    if len(freqs) != 24:
        print "There are not 24 coarse chans defined for this obs. Got: %s" % freqs
        return

    freqs.sort()   # It should already be sorted, but just in case...
    lowchans = [f for f in freqs if f <= 128]
    highchans = [f for f in freqs if f > 128]
    highchans.reverse()
    freqs = lowchans + highchans
    return freqs

def getmeta(service='obs', params=None):
    """
    Function to call a JSON web service and return a dictionary:
    Given a JSON web service ('obs', find, or 'con') and a set of parameters as
    a Python dictionary, return a Python dictionary containing the result.
    Taken verbatim from http://mwa-lfd.haystack.mit.edu/twiki/bin/view/Main/MetaDataWeb
    """

    if params:
        data = urllib.urlencode(params)  # Turn the dictionary into a string with encoded 'name=value' pairs
    else:
        data = ''

    if service.strip().lower() in ['obs', 'find', 'con']:
        service = service.strip().lower()
    else:
        print "invalid service name: %s" % service
        return

    try:
        result = json.load(urllib2.urlopen(BASEURL + service + '?' + data))
    except urllib2.HTTPError as error:
        print "HTTP error from server: code=%d, response:\n %s" % (error.code, error.read())
        return
    except urllib2.URLError as error:
        print "URL or network error: %s" % error.reason
        return

    return result

def get_frequencies(obs_id):
    obsinfo = getmeta(service='obs', params={'obs_id':str(obs_id)})
    freq_array = obsinfo['rfstreams']['0']['frequencies']
    return sfreq(freq_array)






#executables < locale specific should probably clean this up >
import distutils.spawn
voltdownload = distutils.spawn.find_executable("voltdownload.py")
prepare = distutils.spawn.find_executable("prepare.py")
recombine = distutils.spawn.find_executable("recombine.py")
make_beam = distutils.spawn.find_executable("make_beam")
#working dir < ditto >
working_root = "notset"
corrdir = "notset"

#first second
start_time = 1380056664
#stop second
stop_time = 1380060293
#increment
increment = 200
#obsid
obsid = 1064091848
#obsname
obsname="D0004"
#

#PSRFITS MODE

#beam_mode=1


#VDIF MODE

beam_mode=0

# number of "threads per CPU"

jobs_per_node = 8

#chan_list (half BW)

#chan_list_half=["ch01","ch02","ch03","ch04","ch05","ch06","ch07","ch08","ch09","ch10","ch11","ch12"]
chan_list_full=["ch01","ch02","ch03","ch04","ch05","ch06","ch07","ch08","ch09","ch10","ch11","ch12","ch13","ch14","ch15","ch16","ch17","ch18","ch19","ch20","ch21","ch22","ch23","ch24"]
n_coarse = 24
parallel = 3
chan_list = []
# pointing
pointing = " 04:37:15.7 -47:15:08 "
#pointing = " 05:14:06.7 -40:02:48.9 "
# narrow channel count per coarse channel
nchan = 128
#jones
useJones = False
#get fresh  data
getdata = False
#process a single step
single_step = -1
#Run the correlator
runMWAC = False
#Submit the JOB
Go = False
# Run the recombiner
runRECOMBINE = False
#buildPFB
runPFB = True

# set some initial values



def options (opts={}):

    print "\noptions:\n"
    print "-B [1/0]\t Submit download jobs to the copyq - at the moment this mode will only download and will perform <NO> subsequent processing [%d] \n" % (opts['batch_download'])
    print "-b:\t UNIX time of the beginning [%d]]\n" % (opts['begin'])
    print "-c:\t Coarse channel count (how many to process) [%d]\n" % (opts['ncoarse_chan'])
    print "-d:\t Number of parallel downloads to envoke if using '-g' [%d]\n" % (opts['parallel_dl'])
    print "-e:\t UNIX time of the end [%d]\n" % (opts['end'])
    print "-g:\t Get the data? (True/False) add this to get fresh data from the archive [%s]\n" % (opts['get_data'])
    print "-i:\t Increment in seconds (how much we process at once) [%d]\n" % (opts['inc'])
    print "-j:\t [corrdir] Use Jones matrices from the RTS [%s,%s]\n" % (opts['useJones'],opts['corrdir'])
    print "-m:\t Beam forming mode (0 == NO BEAMFORMING 1==PSRFITS, 2==VDIF) [%d]\n" % (opts['mode'])
    print "-n:\t Number of fine channels per coarse channel [%d]\n" % (opts['nchan'])
    print "-o:\t obsid [%s]\n" % opts['obsid']
    print "-p:\t beam pointind [%s]\n" % opts['pointing']
    print "-s:\t single step (only process one increment and this is it (-1 == do them all) [%d]\n" % opts['single_step']
    print "-r:\t [corrdir] Run the offline correlator - this will submit a job to process the .dat files into visibility sets into the specified directory. These are needed if you want an RTS calibration solution [%s]\n" % opts['corrdir']
    print "-G:\t Submit the beamformer/correlator job [Do it = %s]\n" % opts['Go']
    print "-R:\t New VCS mode - requires the recombine operation [runRECOMBINE = %s]\n" % opts['runRECOMBINE']
    print "-w:\t Working root directory [%s]\n" % opts['root']
    print "-z:\t Add to switch off PFB formation/testing [runPFB = %s]\n" % opts['runPFB']


def usage (opts={}):

    print "process_all.py is a script of scripts that downloads prepares and submits jobs to Galaxy. It can be run with just a pointing (-p \"xx:xx:xx xx:xx:xx.x\") and an obsid (\"-o <obsid>\") and it will process all the data in the obs. It will call prepare.py which will attempt to build the phase and calibration information - which will only exist if a calibration obs has already been run. So it will only move past the data prepa stage if the \"-r\" (for run) is used\n"

    options(opts)

    sys.exit()


if __name__ == '__main__':

    the_options = {'begin': start_time, 'ncoarse_chan' : n_coarse, 'end' : stop_time, 'get_data':getdata, 'parallel_dl':parallel, 'inc':increment,'useJones':useJones, 'mode': beam_mode, 'nchan':nchan, 'obsid': obsid, 'pointing' : pointing, 'single_step' : single_step, 'runPFB' : runPFB, 'runMWAC': runMWAC, 'corrdir': corrdir, 'Go':Go, 'runRECOMBINE' : runRECOMBINE, 'root' : working_root, 'batch_download' : 0}

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hB:b:c:e:gd:Gi:j:m:n:o:p:r:Rs:w:z")
    except getopt.GetoptError:
        usage(the_options)
        sys.exit()
    finally:
        if len(sys.argv) < 2:
            usage(the_options)

#print opts

    for opt,arg in opts:

        if (opt == "-h"):
            usage(the_options)
        elif (opt == "-B"):
            the_options['batch_download'] = int(arg)
            the_options['get_data'] = True
            the_options['runRECOMBINE'] = True
        elif (opt == "-b"):
            the_options['begin'] = int(arg)
        elif (opt == "-c"):
            the_options['ncoarse_chan'] = int(arg)
        elif (opt == "-e"):
            the_options['end'] = int(arg)
        elif (opt == "-d"):
            the_options['parallel_dl'] = int(arg)
        elif (opt == "-g"):
            the_options['get_data'] = True
        elif (opt == "-i"):
            the_options['inc'] = int(arg)
        elif (opt == "-j"):
            the_options['useJones'] = True
            the_options['corrdir'] = arg
        elif (opt == "-m"):
            the_options['mode'] = int(arg)
        elif (opt == "-n"):
            the_options['nchan'] = int(arg)
        elif (opt == "-o"):
            the_options['obsid'] = int(arg)
        elif (opt == "-p"):
            the_options['pointing'] = arg
        elif (opt == "-r"):
            the_options['runMWAC'] = True
            the_options['corrdir'] = arg
        elif (opt == "-s"):
            the_options['single_step'] = int(arg)
        elif (opt == "-G"):
            the_options['Go'] = True
        elif (opt == "-R"):
            runRECOMBINE = True
            the_options['runRECOMBINE'] = True
        elif (opt == "-w"):
            the_options['root'] = arg
        elif (opt == "-z"):
            the_options['runPFB'] = False
            runPFB=False


    options (the_options)
    if (the_options['root'] == working_root):
        print "Please set working root with -w\n"
        sys.exit(1)
    if (the_options['runMWAC'] and the_options['corrdir'] == corrdir):
        print "Please set the correlator output dir with -r\n"
        sys.exit(1)

   #    import pdb
#    pdb.set_trace()

    # this is the beamformer mode
    # -f is the psrfits mode and -v is the vdif mode
    if (the_options['mode'] == 1):
        beam_mode_str = "-f"
    elif (the_options['mode'] == 2):
        beam_mode_str = "-v"


    # this is the direction dependent calibration 
    # you do not need this - but you should have it

    if (the_options['useJones'] == True):
        jones = "-j jones.txt"
    else:
        jones = " -i "

    # number of coarse channels to process

    if ((the_options['ncoarse_chan']) == 0):
        skip = "-c"
    else:
        skip = " "

    working_root = the_options['root']
    obsid = the_options['obsid']
    start_time = the_options['begin']
    stop_time = the_options['end']
    increment = the_options['inc']
    getdata = the_options['get_data']
    parallel = the_options['parallel_dl']
    Go = the_options['Go']
    runRECOMBINE = the_options['runRECOMBINE']
    runMWAC = the_options['runMWAC']
    pointing = the_options['pointing']
    beam_mode = the_options['mode']
    nchan = the_options['nchan']

    make_dir = "mkdir %s" % working_root
    subprocess.call(make_dir,shell=True);
    working_dir = "%s/%s" % (working_root,obsid)
    make_dir = "mkdir %s" % working_dir
    subprocess.call(make_dir,shell=True);

    metafits_file = "%s/%d.metafits" % (working_dir,obsid)

    if (os.path.isfile(metafits_file) == False):
        metafile_line = "wget  http://ngas01.ivec.org/metadata/fits?obs_id=%d -O %s\n" % (obsid,metafits_file);
        subprocess.call(metafile_line,shell=True);


    chan_list = get_frequencies(obsid)[0:the_options['ncoarse_chan']]

    print "Will process the following channels:\n"
    print chan_list

    print "Set starting options\n"


    submitted_jobs = []
    submitted_times = []

#get data piecemeal and process it

    step = 0
    last_increment = 0

    print "start: %d Stop: %d Inc: %d \n" % (int(start_time),int(stop_time),int(increment))

    for time_to_get in range(int(start_time),int(stop_time),int(increment)):
        
        # this refers to recombine tasks - I reset it here
        # so that only one increment is popped off the stack for beamforming
        # each time

        a_job_is_done = False;
        print "Time to get is : %s\n" % time_to_get

        if (time_to_get + int(increment) >= int(stop_time)):
            last_increment = 1

        try:
            os.chdir(working_root)
        except:
            print "cannot open working (root) dir:%s" % working_root
            sys.exit()

        print "processing step %d (time %d)\n" % (step, time_to_get)
        step = step + 1

# this little piece of logic lets you advance through the steps w/o processing 
# it is not really required as you can just change the start/stop

        if (the_options['single_step'] > 0):
            if (step != the_options['single_step']):
                continue
    
# Are we downloading data

        if (the_options['get_data'] == True):

# builds the voltdownload command line
            if (runRECOMBINE == False):
                get_data = "%s --obs=%s --type=12 --from=%d --duration=%d --parallel=%d " % (voltdownload,obsid,time_to_get,increment-1,parallel)
            else:
                get_data = "%s --obs=%s --type=11 --from=%d --duration=%d --parallel=%d " % (voltdownload,obsid,time_to_get,increment-1,parallel)

# download everything via the copy queue - do not do anything else
            if (the_options['batch_download'] == 1):
                voltdownload_batch = "%s/volt_%d.batch" % (working_dir,time_to_get)
                secs_to_run = datetime.timedelta(seconds=120*increment)
                with open(voltdownload_batch,'w') as batch_file:

                    batch_line = "#!/bin/bash -l\n\n"
                    batch_file.write(batch_line)
                    batch_line = "%s\n" % (get_data)
                    batch_file.write(batch_line)
            
                submit_line = "sbatch --time=%s --workdir=%s -M zeus --partition=copyq %s\n" % (str(secs_to_run),working_root,voltdownload_batch)
                submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
                continue
            else:
                submit_cmd = subprocess.Popen(get_data,shell=True,stdout=subprocess.PIPE)


            try:
                os.chdir(working_dir)
            except:
                print "cannot open working dir:%s" % working_dir
                sys.exit()

# do I need to recombine this batch

        if (runRECOMBINE == True):



            recombine_batch = "%s/recombine_%d.batch" % (working_dir,time_to_get)




            with open(recombine_batch,'w') as batch_file:
            
                nodes = int(int(increment)/jobs_per_node) + 1


                batch_line = "#!/bin/bash -l\n#SBATCH --time=00:10:00\n#SBATCH \n#SBATCH --export=NONE\n#SBATCH --nodes=%d\n" % (nodes)


                batch_file.write(batch_line)
                batch_line = "module load mpi4py\n"
                batch_file.write(batch_line)
                batch_line = "module load cfitsio\n"
                batch_file.write(batch_line)


                if (jobs_per_node > increment):
                    jobs_per_node = increment

                recombine_line = "aprun -n %d -N %d python %s %s -o %s -s %d -w %s\n" % (increment,jobs_per_node,recombine,skip,obsid,time_to_get,working_dir)

                batch_file.write(recombine_line)



            submit_line = "sbatch --partition=gpuq %s\n" % (recombine_batch)

            submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
            jobid=""
            for line in submit_cmd.stdout:

                if "Submitted" in line:
                    (word1,word2,word3,jobid) = line.split()
                    if (is_number(jobid)): 
                        submitted_jobs.append(jobid)
                        submitted_times.append(time_to_get)
        # we have submitted the recombine job for this increment
        # end if get_data == true
        if (the_options['get_data'] == False):
        # We are not downloading data
        # move into the working directory
            try:
                os.chdir(working_dir)
            except:
                print "cannot open working dir:%s" % working_dir
                sys.exit()

        # this is our current time
            ttg = time_to_get
        # this assumes a recombine job is already completed    
            a_job_is_done = True

        if (runRECOMBINE == True and (the_options['runMWAC'] == True or the_options['mode'] !=0)):
        
        # we are in REcombine mode / we may - or may not have downloaded the data
        # if we have not downloaded the data for recombine to process - some subsequent jobs require
        # these submissions to be complete.
        # We cannot simply make them a dependency as they require the data to be already present before 
        # prepare.py can be run.
        # One option is to run prepare.py as an sbatch submission - which may (or may not) work at the moment

            print submitted_jobs
            for entry,jobid in enumerate(submitted_jobs):
        # but we have submitted a recombine job              
        # now we have to wait until a job is finished before we move on
        # interrogate the queue for the first job
                queue_line = "squeue -j %s\n" % jobid
                queue_cmd = subprocess.Popen(queue_line,shell=True,stdout=subprocess.PIPE)
        # assume the job is finished
                finished = True
                for line in queue_cmd.stdout:

                    if jobid in line:
        # batch job still in the queue so we are not finished
                        finished = False;
        # we test for a_job_is_done because we only one 1 job popped off the stack
                    if (finished == True and a_job_is_done == False):
                        submitted_jobs.pop(entry)
                        # ttg now holds the time of the completed recombine job
                        ttg = submitted_times.pop(entry)
                        a_job_is_done = True
            
            if (last_increment == 1):
                # are we on the last increment - there is no more data to load - so we should just wait till a recombine finishes
                while (len(submitted_jobs) > 0):
                    time.sleep(1)
                    for entry,jobid in enumerate(submitted_jobs):
                        #now we have to wait until this job is finished before we move on
                        queue_line = "squeue -j %s\n" % jobid
                        queue_cmd = subprocess.Popen(queue_line,shell=True,stdout=subprocess.PIPE)
                        finished = True
                        for line in queue_cmd.stdout:

                            if jobid in line:
                            # batch job still in the queue
                                finished = False;
                        if (finished == True):
                            submitted_jobs.pop(entry)
                        # we only want ttg to be the first one - not to get clobbered.
                        # but we want to wait until all the jobs are finished before we move on becuase we may not
                        # come back ...
                        if (finished == True and a_job_is_done == False):
                            ttg = submitted_times.pop(entry)
                            a_job_is_done = True
        else:
            ttg=time_to_get
            a_job_is_done = True


# now process

        if ((Go == True) and (a_job_is_done == True)):
            moved = 0
            pfb_job_list = []
            # are we going to form pfb files - this is a bit of a misnomer as it is actually un-pfbs
            if (runPFB == True):
                # list of submitted pfb batch jobs
                for index,channel in enumerate(chan_list):
                    # pfbfile batch file one for each channel
                    pfb_batch_file = "%s/pfb_build_ch%02d.batch" % (working_dir,index+1)
                    # we need to open the file

                    f=[]
                    with open(pfb_batch_file, 'w') as pfb_build:
                        pfb_build.write("#!/bin/bash -l\n")

                        nodes_line = "#SBATCH --nodes=1\n#SBATCH --export=NONE\n" 
                        pfb_build.write(nodes_line)
                        
                 
                        
                        # we (un)pfb <all> available combined files for this channel
                        # we first move them to the target channel directory
                        # this means the correlator batch does not depend on the PFB ones
                        # but the BF does....

                        files_glob = "%s/combined/*_ch%s*" % (working_dir,channel)
                        for to_convert in sorted(glob.glob(files_glob)):
                            f.append(to_convert)

                        # make the output dir for the channel if not already done

                        make_dir = "mkdir %s/ch%02d" % (working_dir,(index+1))
                        subprocess.call(make_dir,shell=True)
                       
                        to_pfb = 0;
                        for datfile in f:
                            # this is now the input file
                            infile = datfile
                            # the outputfile
                            localdone = "%s.pfb" % datfile
                            # the remote copy
                            donefile = "%s/ch%02d/%s" % (working_dir,(index+1),os.path.basename(localdone))
                            # does the file still exist
                            cp_cmd=""
                            pfb_line=""
                            if (os.path.isfile(infile) == True and os.path.isfile(donefile) == True):
                                # there is an imput file and an output file
                                # is the output file the correct size
                                file_statinfo = os.stat(infile)
                                done_file_statinfo = os.stat(localdone)
                                if (done_file_statinfo.st_size == 2*file_statinfo):
                                    pfb_line = "#read_pfb -i %s -a 128 -n 128  -o %s -4 \n" % (infile,localdone)
                                    cp_cmd = "cp %s %s\n" % (localdone,donefile)
                                elif (os.path.isfile(infile) == True):
                                    # the output file is the wrong size/doesn't exit - try again
                                    pfb_line = "read_pfb -i %s -a 128 -n 128  -o %s -4 \n" % (infile,localdone)
                                    cp_cmd = "cp %s %s\n" % (localdone,donefile)
                                    to_pfb = to_pfb + 1
                            elif (os.path.isfile(infile) == True):
                                # the output file is the wrong size/doesn't exit - try again
                                pfb_line = "read_pfb -i %s -a 128 -n 128  -o %s -4 \n" % (infile,localdone)
                                cp_cmd = "cp %s %s\n" % (localdone,donefile)
                                to_pfb = to_pfb + 1

                            else:
                                print "Cannot find %s" % infile
                                missing_files = "%s/missing.list" % (working_dir)
                            
                                with open(missing_files,'a') as missing:
                                    missing_line = "%s\n" % (infile)
                                    missing.write(missing_line)

                            pfb_build.write(pfb_line)
                            pfb_build.write(cp_cmd)
                    # submit the job

                    secs_to_run = datetime.timedelta(seconds=30*to_pfb)
                    submit_line = "sbatch --time=%s --nodes=1 --workdir=%s --partition=gpuq %s\n" % (str(secs_to_run),working_dir,pfb_batch_file)

                    if (secs_to_run.seconds > 0):
                        print submit_line
                        submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
                        jobid=""
                        for line in submit_cmd.stdout:
                            if "Submitted" in line:
                                (word1,word2,word3,jobid) = line.split()
                                if (is_number(jobid)):
                                    pfb_job_list.append(jobid)
            # all the PFB jobs and associated clean ups have been submitted
            # if this is all we are doing go to the next increment
            # -- but we only want to do this once the 
                
                sys.exit();

            if ((the_options['mode'] == 0) and (runMWAC == False)):
                continue

            channel = 0
            children=[]
            child=0
            
            # this will copy any recombined files to the channel directories
            # if the (un)pfbs are being made this will have happened already

            for index,channel in enumerate(chan_list):
                print "processing %s\n" % channel
                
                channel_dir = "%s/ch%02d" % (working_dir,(index+1))
            
                if (the_options['runMWAC'] == True):
                    make_dir = "mkdir %s" % (channel_dir)
                    subprocess.call(make_dir,shell=True)
                    f=[]
                    files_glob = "%s/combined/*_ch%s*" % (working_dir,channel)
                    for to_move in sorted(glob.glob(files_glob)):
                        f.append(to_move)

                    for file in f:
                        mv_cmd = "mv %s %s/\n" % (file,channel_dir)
                        subprocess.call(mv_cmd,shell=True)

                try:
                    os.chdir(channel_dir)
                except:
                    print "cannot open channel dir:%s" % channel_dir
                    sys.exit()

                (ra,dec) = pointing.split()

                if (runMWAC == True):
                    
                    prepare_line = "%s -r %s -d %s -s %s -f %s" % (prepare,ra,dec,the_options['corrdir'],metafits_file)
                else:
                    prepare_line = "%s -r %s -d %s -g %s -f %s" % (prepare,ra,dec,the_options['corrdir'],metafits_file)

                if (the_options['nchan'] == 88):
                    prepare_line += " -m 0 "
                else:
                    prepare_line += " -m 1 "

                if (runMWAC == True):
                    prepare_line += " -e dat "

                print "Will launch prepare by: %s\n" % prepare_line 
                try:
                    child = os.fork()

                except:
                    print "Error on fork"
                    sys.exit()

                if (child == 0):
                    subprocess.call(prepare_line,shell=True)
                    sys.exit()

                else:
                    time.sleep(1)
                    children.append(child)

                try:
                    os.chdir(working_dir)
                except:
                    print "cannot return to working dir :%s\n" % (working_dir)
                    sys.exit()

            # waiting for prepare to finish
            for i,child in enumerate(children):
                print "Waiting for child %d:%d" % (i,child)
                os.waitpid(child,0)
            
            if (runMWAC == False):
                to_beam = 0
                for index,channel in enumerate(chan_list):
                    print "Checking %s\n" % channel
                    channel_dir = "%s/ch%02d" % (working_dir,(index+1))
                    flags_file = "%s/flags.txt" % channel_dir
                    if (os.path.isfile(flags_file) == True):
                        print "Channel %s passed\n" % channel

                        if (index == 0):
                            files_glob = "%s/*.pfb" % (channel_dir)
                            for entry in sorted(glob.glob(files_glob)):
                                to_beam = to_beam+1


                    else:
                        print "Channel %s failed exiting to avoid confusion" % channel
                        sys.exit()
    # now actually submit the job (5 seconds per second) only if not running the correlator

                
                secs_to_run = datetime.timedelta(seconds=10*to_beam)
                number_of_exe = len(chan_list)
                exe_per_node = 1
                queue = "gpuq"

                batch = "%s_%02d.sh" % (obsid,step)

                with open(batch, 'w') as batch_file:
                    batch_file.write("#!/bin/bash -l\n")

                    nodes_line = "#SBATCH --nodes=%d\n#SBATCH --export=NONE\n" % (number_of_exe/exe_per_node)
                    batch_file.write(nodes_line)
                    time_line = "#SBATCH --time=%s\n" % (str(secs_to_run))
                    batch_file.write(time_line)
                    aprun_line = "aprun -n %d -N %d %s -e pfb -o ch01 -a 128 -n %d -t 1 %s -c phases.txt -w flags.txt -D %s/ch %s psrfits_header.txt\n" % (number_of_exe,exe_per_node,make_beam,nchan,jones,working_dir,beam_mode_str)
                    batch_file.write(aprun_line)

                submit_line = "sbatch --nodes=%d --workdir=%s --time=%s --partition=%s %s\n" % (number_of_exe/exe_per_node,working_dir,str(secs_to_run),queue,batch)
                print submit_line

                submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
                jobid=""
                for line in submit_cmd.stdout:
                    if "Submitted" in line:
                        (word1,word2,word3,jobid) = line.split()


            if (the_options['get_data'] == False):
                print "finished"
                sys.exit()

        # end the beamformer for this iteration   
        # end the "Go" section for this iteration

