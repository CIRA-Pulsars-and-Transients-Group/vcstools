#!/usr/bin/env python
import subprocess
import os
import sys
import glob
import time
import datetime
import distutils.spawn
from astropy.io import fits as pyfits

def getmeta(service='obs', params=None):
    """
    Function to call a JSON web service and return a dictionary:
    Given a JSON web service ('obs', find, or 'con') and a set of parameters as
    a Python dictionary, return a Python dictionary containing the result.
    Taken verbatim from http://mwa-lfd.haystack.mit.edu/twiki/bin/view/Main/MetaDataWeb
    """
    import urllib
    import urllib2
    import json

    # Append the service name to this base URL, eg 'con', 'obs', etc.
    BASEURL = 'http://ngas01.ivec.org/metadata/'


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

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def mdir(path,description):
    try:
        os.mkdir(path)
        os.chmod(path,0761)
    except:
        if (os.path.exists(path)):
            print "{0} Directory Already Exists\n".format(description)
        else:
            sys.exit()

def ensure_metafits(metafits_file):
        if (os.path.isfile(metafits_file) == False):
            metafile_line = "wget  http://ngas01.ivec.org/metadata/fits?obs_id=%d -O %s\n" % (opts.obs,metafits_file)
            subprocess.call(metafile_line,shell=True)


def obs_max_min(obs_id):
    """
    Small function to query the database and returns the times of the first and last file
    :param obs_id:
    :return:
    """
    from file_maxmin import getmeta

    obsinfo = getmeta(service='obs', params={'obs_id':str(obs_id)})

    obs_start = int(min(obsinfo['files'])[11:21])
    obs_end = int(max(obsinfo['files'])[11:21])
    return obs_start, obs_end

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


def get_frequencies(metafits):
    hdulist = pyfits.open(metafits)
    freq_array = hdulist[0].header['CHANNELS']
    return sfreq(freq_array.split(','))


def options (options): # TODO reformat this to print properly

    print "\noptions:\n"
    print "--mode {0}".format(options.mode)
    print "-B [1/0]\t Submit download jobs to the copyq - at the moment this mode will only download and will perform <NO> subsequent processing [%d] \n" % (opts['batch_download'])
    print "-b:\t GPS/UNIX time of the beginning [%d]]\n" % (opts['begin'])
    print "-c:\t Coarse channel count (how many to process) [%d]\n" % (opts['ncoarse_chan'])
    print "-d:\t Number of parallel downloads to envoke if using '-g' [%d]\n" % (opts['parallel_dl'])
    print "-e:\t GPS/UNIX time of the end [%d]\n" % (opts['end'])
 #   print "-g:\t Get the data? (True/False) add this to get fresh data from the archive [%s]\n" % (opts['get_data'])
    print "-i:\t Increment in seconds (how much we process at once) [%d]\n" % (opts['inc'])
    print "-j:\t [corrdir] Use Jones matrices from the RTS [%s,%s]\n" % (opts['useJones'],opts['corrdir'])
    print "-m:\t Beam forming mode (0 == NO BEAMFORMING 1==PSRFITS, 2==VDIF) [%d]\n" % (opts['mode'])
    print "-n:\t Number of fine channels per coarse channel [%d]\n" % (opts['nchan'])
    print "-o:\t obsid [%s]\n" % opts['obsid']
    print "-p:\t beam pointing [%s]\n" % opts['pointing']
    print "-s:\t single step (only process one increment and this is it (-1 == do them all) [%d]\n" % opts['single_step']
#    print "-r:\t [corrdir] Run the offline correlator - this will submit a job to process the .dat files into visibility sets into the specified directory. These are needed if you want an RTS calibration solution [%s]\n" % opts['corrdir']
    print "-G:\t Submit the beamformer/correlator job [Do it = %s]\n" % opts['Go']
#   print "-R:\t New VCS mode - requires the recombine operation [runRECOMBINE = %s]\n" % opts['runRECOMBINE']
    print "-w:\t Working root directory [%s]\n" % opts['root']
#    print "-z:\t Add to switch off PFB formation/testing [runPFB = %s]\n" % opts['runPFB']


def vcs_download(obsid, start_time, stop_time, increment, copyq, format, working_dir, parallel):
    print "Downloading files from archive"
#    voltdownload = distutils.spawn.find_executable("voltdownload.py")
    voltdownload = "/group/mwaops/stremblay/MWA_CoreUtils/voltage/scripts/voltdownload.py"
#   voltdownload = "python /home/fkirsten/software/galaxy-scripts/scripts/voltdownload.py"
    raw_dir = "{0}/raw".format(working_dir)
    mdir(raw_dir, "Raw")

    for time_to_get in range(start_time,stop_time,increment):
        get_data = "{0} --obs={1} --type={2} --from={3} --duration={4} --parallel={5} --dir={6}".format(voltdownload,obsid, format, time_to_get,(increment-1),parallel, raw_dir)
        if copyq:
            voltdownload_batch = "{0}/batch/volt_{1}.batch".format(working_dir,time_to_get)
            secs_to_run = datetime.timedelta(seconds=300*increment)
            with open(voltdownload_batch,'w') as batch_file:

                batch_line = "#!/bin/bash -l\n#SBATCH --export=NONE\n#SBATCH --output={0}/batch/volt_{1}.out\n".format(working_dir,time_to_get)
                batch_file.write(batch_line)
                batch_line = "%s\n" % (get_data)
                batch_file.write(batch_line)

            submit_line = "sbatch --time={0} --workdir={1} -M zeus --partition=copyq {2}\n".format(secs_to_run,raw_dir,voltdownload_batch)
            submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
            continue
        else:
            log_name="{0}/voltdownload_{1}.log".format(working_dir,time_to_get)
            with open(log_name, 'w') as log:
                subprocess.call(get_data, shell=True, stdout=log, stderr=log)


        try:
            os.chdir(working_dir)
        except:
            print "cannot open working dir:{0}".format(working_dir)
            sys.exit()
    check = "checks.py -m download -o {0}".format(obsid)
    

def vcs_recombine(obsid, start_time, stop_time, increment, working_dir):
    print "Running recombine on files"
    jobs_per_node = 8
#    recombine = distutils.spawn.find_executable("recombine.py")
    recombine = "/group/mwaops/stremblay/galaxy-scripts/scripts/recombine.py"
    for time_to_get in range(start_time,stop_time,increment):

        recombine_batch = "{0}/batch/recombine_{1}.batch".format(working_dir,time_to_get)
        with open(recombine_batch,'w') as batch_file:

            nodes = (increment+(-increment%jobs_per_node))//jobs_per_node + 1 # Integer division with ceiling result plus 1 for master node

            batch_line = "#!/bin/bash -l\n#SBATCH --time=06:00:00\n#SBATCH \n#SBATCH --output={0}/batch/recombine_{1}.out\n#SBATCH --export=NONE\n#SBATCH --nodes={2}\n".format(working_dir, time_to_get, nodes)

            batch_file.write(batch_line)
            batch_line = "module load mpi4py\n"
            batch_file.write(batch_line)
            batch_line = "module load cfitsio\n"
            batch_file.write(batch_line)

            if (stop_time - time_to_get) < increment:       # Trying to stop jobs from running over if they aren't perfectly divisible by increment
                increment = stop_time - time_to_get + 1

            if (jobs_per_node > increment):
                jobs_per_node = increment

            recombine_line = "aprun -n {0} -N {1} python {2} -o {3} -s {4} -w {5}\n".format(increment,jobs_per_node,recombine,obsid,time_to_get,working_dir)
            batch_file.write(recombine_line)

        submit_line = "sbatch --partition=gpuq --workdir={0} {1}\n".format(working_dir,recombine_batch)

        submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        jobid=""
        for line in submit_cmd.stdout:

            if "Submitted" in line:
                (word1,word2,word3,jobid) = line.split()
#                if (is_number(jobid)):
#                    submitted_jobs.append(jobid)
#                    submitted_times.append(time_to_get)



def vcs_correlate(obsid,start,stop,increment,working_dir):
    print "Correlating files"
    import os
    import astropy
    from astropy.time import Time
    import datetime
    import calendar

    corr_dir = "{0}/vis".format(working_dir)
    mdir(corr_dir, "Correlator Product")

    chan_list = get_frequencies(metafits_file)

    for time_to_get in range(start,stop,increment):
        inc_start = time_to_get
        inc_stop = time_to_get+increment
        for index,channel in enumerate(chan_list):
            gpubox_label = (index+1)
            f=[]
            for time_to_corr in range(inc_start,inc_stop,1):
                file_to_process = "{0}/combined/{1}_{2}_ch{3:0>2}.dat".format(working_dir,obsid,time_to_corr,channel)
                #check the file exists
                if (os.path.isfile(file_to_process) == True):
                    f.append(file_to_process)

            #now have a full list of files
            #for this increment 
            #and this channel
            if (len(f) > 0):
                corr_batch = "{0}/batch/correlator_{1}_gpubox{2:0>2}.batch".format(working_dir,inc_start,gpubox_label)

                with open(corr_batch, 'w') as batch_file:
                    batch_file.write("#!/bin/bash -l\n#SBATCH --nodes=1\n#SBATCH --export=NONE\n #SBATCH --output={0}.out\n".format(corr_batch))
                    batch_file.write("module load cudatoolkit\nmodule load cfitsio\n")
                
                to_corr = 0
                for file in f:
                    corr_line = ""
                    (current_time,ext) = os.path.splitext(os.path.basename(file))
                    (obsid,gpstime,chan) = current_time.split('_')
                    t = Time(int(gpstime), format='gps', scale='utc')
                    time_str =  t.datetime.strftime('%Y-%m-%d %H:%M:%S')

                    current_time = time.strptime(time_str, "%Y-%m-%d  %H:%M:%S")
                    unix_time = calendar.timegm(current_time)

                    corr_line = " aprun -n 1 -N 1 {0} -o {1} -s {2}/{3} -r 1 -i 100 -f 128 -n 4 -c {4:0>2} -d {4}\n".format("mwac_offline",corr_dir,obsid,unix_time,gpubox_label,file)
                    
                    with open(corr_batch, 'a') as batch_file:
                        batch_file.write(corr_line)
                        to_corr = to_corr+1

                secs_to_run = datetime.timedelta(seconds=10*to_corr)
                batch_submit_line = "sbatch --workdir={0} --time={1} --partition=gpuq {2}\n".format(corr_dir,secs_to_run,corr_batch)
                submit_cmd = subprocess.Popen(batch_submit_line,shell=True,stdout=subprocess.PIPE)
                jobid=""
                for line in submit_cmd.stdout:
                    if "Submitted" in line:
                        (word1,word2,word3,jobid) = line.split()





<<<<<<< HEAD
def coherent_beam(obsid,start,stop,increment):
    print "Forming coherent beam"
    print "Calculating Phase Model (running get_delays)"
=======
def coherent_beam(obs_id, working_dir, metafile, nfine_chan, pointing):
    # Need to run get_delays and then the beamformer on each desired coarse channel
    DI_dir = working_dir+"DIJ"
    RA = pointing[0]
    Dec = pointing[1]

    print "Running get_delays"
    P_dir = working_dir+"pointings"
    mdir(P_dir, "Pointings")
    pointing_dir = "{0}/{1}_{2}".format(P_dir, RA, Dec)
    mdir(pointing_dir, "Pointing {0} {1}".format(RA, Dec))

    for gpubox in ["{0:0>2}".format(i) for i in range(1,25)]:
        #DI_file = "{0}/{1}".format(DI_dir, ?) # Need to finish file path
        pointing_chan_dir = "{0}/{1}".format(pointing_dir,gpubox)
        mdir(pointing_chan_dir, "Pointing {0} {1} gpubox {2}".format(RA, Dec, gpubox))

        DI_file = "{0}/DI_JonesMatrices_node{1}.dat".format(DI_dir, gpubox)
        if (os.path.isfile(DI_file)):
            get_delays_batch = "{0}/batch/gd_{1}.batch".format(working_dir,gpubox)
            with open(get_delays_batch,'w') as batch_file:
                batch_line = "#!/bin/bash -l\n#SBATCH --export=NONE\n#SBATCH --output={0}/batch/gd_{1}.out\n".format(working_dir,gpubox)
                batch_file.write(batch_line)
                #delays_line = "get_delays -a {0} -b {1} -j {2} -m {3} -i -p -z {4} -o {5} -f {6} -n {7} -w 10000 -r {8} -d {9}\n".format(pointing_chan_dir,?,DI_file,metafile,utctime,obs_id,?,nfine_chan,Dec) # need to finish inputs
                batch_file.write(delays_line)
            submit_line = "sbatch --time={0} --workdir={1} --partition=gpuq {2}\n".format(time_to_run, pointing_chan_dir, get_delays_batch)
            submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        else:
            print "WARNING: No Calibration Found for Channel {0}!".format(gpubox)


    # if (the_options['delays'] == True):
    #
    #     DI_file = "%s/DI_JonesMatrices_node0%02d.dat" % (outdir,gpubox_label)
    #     print DI_file
    #
    #     if (old_mode == 1):
    #         if (os.path.isfile(DI_file)):
    #             delays_line = "%s -a ./ -b %d -j %s -m %s -i -p -z %s -o %s -f %s -e %d -n 88 -w 10000 -r %s -d %s" % (get_delays,len(f),DI_file,the_options['metafile'],utctime,obsid,freq_Hz,edge_num,the_options['ra'],the_options['dec'])
    #         else:
    #             delays_line = "%s -a ./ -b %d -i -p -z %s -o %s -f %s -e %d -n 88 -w 10000 -r %s -d %s -m %s" % (get_delays,len(f),utctime,obsid,freq_Hz,edge_num,the_options['ra'],the_options['dec'],the_options['metafile'])
    #     elif (new_mode == 1):
    #         if (os.path.isfile(DI_file)):
    #             delays_line = "%s -a ./ -b %d -j %s -m %s -i -p -z %s -o %s -f %s -n 128 -w 10000 -r %s -d %s" % (get_delays,len(f),DI_file,the_options['metafile'],utctime,obsid,freq_Hz,the_options['ra'],the_options['dec'])
    #         else:
    #             print "WARNING NOT CALIBRATION FOUND\n"
    #             delays_line = "%s -a ./ -b %d -i -p -z %s -o %s -f %s -n 128 -w 10000 -r %s -d %s -m %s" % (get_delays,len(f),utctime,obsid,freq_Hz,the_options['ra'],the_options['dec'],the_options['metafile'])


    print "Forming coherent beam"

    # Run make_beam
    """
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
    """

>>>>>>> 353de5c56b2a0b3e127926d7abaee782d2ceeb0c


if __name__ == '__main__':

    modes=['download','recombine','correlate','beamform']
    jobs_per_node = 8
    chan_list_full=["ch01","ch02","ch03","ch04","ch05","ch06","ch07","ch08","ch09","ch10","ch11","ch12","ch13","ch14","ch15","ch16","ch17","ch18","ch19","ch20","ch21","ch22","ch23","ch24"]
    chan_list = []


    from optparse import OptionParser, OptionGroup

 #   parser=OptionParser(description="process_vcs.py is a script of scripts that downloads prepares and submits jobs to Galaxy. It can be run with just a pointing (-p \"xx:xx:xx xx:xx:xx.x\") and an obsid (\"-o <obsid>\") and it will process all the data in the obs. It will call prepare.py which will attempt to build the phase and calibration information - which will only exist if a calibration obs has already been run. So it will only move past the data prepa stage if the \"-r\" (for run) is used\n"

    parser=OptionParser(description="process_vcs.py is a script for processing the MWA VCS data on Galaxy in steps. It can download data from the archive, call on recombine to form course channels, run the offline correlator, make tile re-ordered and bit promoted PFB files or for a coherent beam for a given pointing.")
    group_download = OptionGroup(parser, 'Download Options')
    group_download.add_option("-B", "--copyq", action="store_true", default=False, help="Submit download jobs to the copyq [default=%default]")
    group_download.add_option("--format", type="choice", choices=['11','12'], default='11', help="Voltage data type (Raw = 11, Recombined Raw = 12) [default=%default]")
    group_download.add_option("-d", "--parallel_dl", type="int", default=3, help="Number of parallel downloads to envoke [default=%default]")

    group_recombine = OptionGroup(parser, 'Recombine Options')

    group_correlate = OptionGroup(parser, 'Correlator Options')
    group_correlate.add_option("--ft_res", metavar="FREQ RES,TIME RES", type="int", nargs=2, default=(40,1), help="Frequency (kHz) and Time (s) resolution to run correlator at. [default=%default]")

    group_pfb = OptionGroup(parser, 'PFB Creation Options')

    group_beamform = OptionGroup(parser, 'Beamforming Options')
    group_beamform.add_option("-p", "--pointing", nargs=2, help="R.A. and Dec. of pointing")
    group_beamform.add_option("--bf_mode", type="choice", choices=['0','1','2'], help="Beam forming mode (0 == NO BEAMFORMING 1==PSRFITS, 2==VDIF)")
    group_beamform.add_option("-j", "--useJones", action="store_true", default=False, help="Use Jones matrices from the RTS [default=%default]")

    parser.add_option("-m", "--mode", type="choice", choices=['download','recombine','correlate','beamform'], help="Mode you want to run. {0}".format(modes))
    parser.add_option("-o", "--obs", metavar="OBS ID", type="int", help="Observation ID you want to process [no default]")
    parser.add_option("-b", "--begin", type="int", help="First GPS time to process [no default]")
    parser.add_option("-e", "--end", type="int", help="Last GPS time to process [no default]")
    parser.add_option("-a", "--all", action="store_true", default=False, help="Perform on entire observation span. Use instead of -b & -e. [default=%default]")
    parser.add_option("-i", "--increment", type="int", default=64, help="Increment in seconds (how much we process at once) [default=%default]")
    parser.add_option("-s", action="store_true", default=False, help="Single step (only process one increment and this is it (False == do them all) [default=%default]")
    parser.add_option("-w", "--work_dir", metavar="DIR", default="/scratch/mwaops/vcs/", help="Base directory you want to run from. This will create a folder for the Obs. ID if it doesn't exist [default=%default]")
    parser.add_option("-c", "--ncoarse_chan", type="int", default=24, help="Coarse channel count (how many to process) [default=%default]")
    parser.add_option("-n", "--nfine_chan", type="int", default=128, help="Number of fine channels per coarse channel [default=%default]")
    parser.add_option_group(group_download)
#    parser.add_option_group(group_recombine)
    parser.add_option_group(group_correlate)
#   parser.add_option_group(group_pfb)
    parser.add_option_group(group_beamform)

    (opts, args) = parser.parse_args()

    if opts.all and (opts.begin or opts.end):
        print "Please specify EITHER (-b,-e) OR -a"
        quit()
    elif opts.all:
        opts.begin, opts.end = obs_max_min(opts.obs)
    if not opts.mode:
      print "Mode required {0}. Please specify with -m or --mode.".format(modes)
      quit()
    if not opts.obs:
        print "Observation ID required, please put in with -o or --obs"
        quit()
    if opts.begin > opts.end:
        print "Starting time is after end time"
        quit()
    if opts.mode == "beamformer":
        if not opts.pointing:
            print "Pointing (-p) required in beamformer mode"
            quit()
        #if opts.pointing[0] or opt.pointing[1]


    mdir(opts.work_dir, "Working")
    obs_dir = "{0}/{1}".format(opts.work_dir,opts.obs)
    mdir(obs_dir, "Observation")
    batch_dir = "{0}/batch".format(obs_dir)
    mdir(batch_dir, "Batch")
    metafits_file = "{0}/{1}.metafits".format(obs_dir,opts.obs)

 #   options(opts)
    print "Processing Obs ID {0} from GPS times {1} till {2}".format(opts.obs, opts.begin, opts.end)

    if opts.mode == 'download':
        print opts.mode
        vcs_download(opts.obs, opts.begin, opts.end, opts.increment, opts.copyq, opts.format, obs_dir, opts.parallel_dl)
    elif opts.mode == 'recombine':
        print opts.mode
        ensure_metafits(metafits_file)
        combined_dir = "{0}/combined".format(obs_dir)
        mdir(combined_dir, "Combined")
        vcs_recombine(opts.obs, opts.begin, opts.end, opts.increment, obs_dir)
    elif opts.mode == 'correlate':
        print opts.mode 
        ensure_metafits(metafits_file)
        vcs_correlate(opts.obs, opts.begin, opts.end, opts.increment, obs_dir)
    elif opts.mode == 'beamformer':
        print opts.mode
        ensure_metafits(metafits_file)
        coherent_beam(opts.obs, obs_dir, metafits_file, opts.nfine_chan, opts.pointing)
    else:
        print "Somehow your non-standard mode snuck through. Try again with one of {0}".format(modes)
        quit()


