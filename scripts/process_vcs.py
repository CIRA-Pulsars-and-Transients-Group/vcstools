#!/usr/bin/env python3


import subprocess
import os
import sys
import tempfile
import atexit
import datetime
from time import sleep, strptime, strftime
import distutils.spawn
import sqlite3 as lite
from astropy.io import fits as pyfits
from astropy.time import Time
from reorder_chans import *
from mdir import mdir
import numpy as np
import logging
import glob

#vcstools functions
from job_submit import submit_slurm
import mwa_metadb_utils as meta

from config_vcs import load_config_file

logger = logging.getLogger(__name__)

try:
    DB_FILE = os.environ['CMD_VCS_DB_FILE']
except KeyError:
    logger.warning("Environmental variable CMD_VCS_DB_FILE is not defined so the "
                   "vcs command database will not be used")
    DB_FILE = None
else:
    import database_vcs


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def gps_to_utc(gps):
    # GPS time as is done in timeconvert.py
    from astropy.utils import iers
    iers.IERS_A_URL = 'https://datacenter.iers.org/data/9/finals2000A.all'
    logger.info(iers.IERS_A_URL)
    utctime = Time(gps, format='gps', scale='utc').fits
    # remove (UTC) that some astropy versions leave on the end
    if utctime.endswith('(UTC)'):
        utctime = strptime(utctime, '%Y-%m-%dT%H:%M:%S.000(UTC)')
        utctime = strftime('%Y-%m-%dT%H:%M:%S', utctime)
    else:
        utctime = strptime(utctime, '%Y-%m-%dT%H:%M:%S.000')
        utctime = strftime('%Y-%m-%dT%H:%M:%S', utctime)
    return utctime


def get_user_email():
    command="echo `ldapsearch -x \"uid=$USER\" mail |grep \"^mail\"|cut -f2 -d' '`"
    email = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True).communicate()[0]
    return email.strip()


def find_combined_beg_end(obsid, base_path="/group/mwavcs/vcs/", channels=None):
    """
    looks through the comined files of the input obsid and returns the max and min in gps time
    Input:
        obsid: The MWA observation ID
    Optional Input:
        base_path: the direct path the base vcs directory. Default: /group/mwavcs/vcs/\
        channels: a list of the frequency channel ids. Default None which then gets the
                  from the mwa metadata
    """
    #TODO have some sort of check to look for gaps
    if glob.glob("{0}/{1}/combined/{1}*_ics.dat".format(base_path, obsid)):
        combined_files = glob.glob("{0}/{1}/combined/{1}*_ics.dat".format(base_path, obsid))
    else:
        channels = meta.get_channels(obsid, channels)
        combined_files = glob.glob("{0}/{1}/combined/{1}*_ch{2}.dat".\
                                   format(base_path, obsid, channels[-1]))
    if len(combined_files) > 0:
        comb_times = []
        for comb in combined_files:
            comb_times.append(int(comb.split("_")[1]))
        beg = min(comb_times)
        end = max(comb_times)
    else:
        logger.warn("No combined files on disk for {0}".format(obsid))
        beg = None
        end = None

    return beg, end


def gps_time_lists(start, stop, chunk):
    time_chunks = []
    while (start + chunk) < stop:
        time_chunks.append((start, start + chunk - 1))
        start += chunk
    time_chunks.append((start, stop))
    return time_chunks


def ensure_metafits(data_dir, obs_id, metafits_file):
    # TODO: To get the actual ppds file should do this with obsdownload -o <obsID> -m
    comp_config = load_config_file()

    if not os.path.exists(metafits_file):
        logger.warning("{0} does not exists".format(metafits_file))
        logger.warning("Will download it from the archive. This can take a "
                      "while so please do not ctrl-C.")
        logger.warning("At the moment, even through the downloaded file is "
                       "labelled as a ppd file this is not true.")
        logger.warning("This is hopefully a temporary measure.")
        #obsdownload = distutils.spawn.find_executable("obsdownload.py")

        get_metafits = "wget http://ws.mwatelescope.org/metadata/fits?obs_id={0} -O {1}".format(obs_id, metafits_file)
        try:
            subprocess.call(get_metafits,shell=True)
        except:
            logger.error("Couldn't download {0}. Aborting.".\
                          format(os.basename(metafits_file)))
            sys.exit(0)
        # clean up
        #os.remove('obscrt.crt')
        #os.remove('obskey.key')
    # make a copy of the file in the product_dir if that directory exists
    # if it doesn't we might have downloaded the metafits file of a calibrator (obs_id only exists on /astro)
    # in case --work_dir was specified in process_vcs call product_dir and data_dir
    # are the same and thus we will not perform the copy
    product_dir = data_dir.replace(comp_config['base_data_dir'], comp_config['base_product_dir']) # being pedantic
    if os.path.exists(product_dir) and not os.path.exists(metafits_file):
        logger.info("Copying {0} to {1}".format(metafits_file, product_dir))
        from shutil import copy2
        copy2("{0}".format(metafits_file), "{0}".format(product_dir))

def create_link(data_dir, target_dir, product_dir, link):
    """
    Creates a symbolic link product_dir/link that points to data_dir/target_dir

    Parameters:
    -----------
    data_dir: string
        The absolute path to the base directory of the true location of the files.
        For our uses this is often a scratch partition like /astro on Galaxy
    target_dir: string
        The folder you would like to be linked to
    product_dir: string
        The absolute path of the link you would like to create
    link: string
        The name of the link you would like to create. Often the same as target_dir
    """
    data_dir = os.path.abspath(data_dir)
    product_dir = os.path.abspath(product_dir)
    if data_dir == product_dir:
        # base directories are the same so doing nothing
        return

    # add product_dir and data_dir to link and target_dir respectively
    link = link.replace(product_dir, '') # just in case...
    link = link.replace('/', '')
    link = os.path.join(product_dir, link)
    target_dir = target_dir.replace(data_dir,'')
    if target_dir.startswith("/"):
        target_dir = target_dir[1:]
    target_dir = os.path.join(data_dir, target_dir)

    # check if link exists and whether it already points to where we'd like it to
    if os.path.exists(link):
        if os.path.islink(link):
            if os.readlink(link) == target_dir:
                return
            else:
                logger.warning("The link {0} already exists but points at {1} while you "
                               "asked it to point at {2}. Deleting the link and creating"
                               "a new one".format(link, os.readlink(link), target_dir))
                os.unlink(link)
                os.symlink(target_dir, link)
        else:
            logger.error("{0} is an existing directory and cannot be turned into a link. Aborting...".format(link))
            sys.exit(0)
    else:
        logger.info("Trying to link {0} against {1}".format(link, target_dir))
        os.symlink(target_dir, link)

def get_frequencies(metafits,resort=False):
    # TODO: for robustness, this should force the entries to be 3-digit numbers
    hdulist    = pyfits.open(metafits)
    freq_str   = hdulist[0].header['CHANNELS']
    freq_array = [int(f) for f in freq_str.split(',')]
    if resort:
        return sfreq(freq_array)
    else:
        return freq_array

def vcs_download(obsid, start_time, stop_time, increment, head, data_dir,
                 product_dir, parallel, vcs_database_id,
                 ics=False, n_untar=2, keep="", vcstools_version="master",
                 nice=0):

    logger.info("Downloading files from archive")
    # voltdownload = distutils.spawn.find_executable("voltdownload.py") #Doesn't seem to be working on zeus for some reason
    voltdownload = "voltdownload.py"
    obsinfo = meta.getmeta(service='obs', params={'obs_id':str(obsid)})
    data_format = obsinfo['dataquality']
    if data_format == 1:
        target_dir = link = '/raw'
        if ics:
            logger.error("Data have not been recombined in the "
                         "archive yet. Exiting")
            sys.exit(0)
        data_type = 11
        dl_dir = "{0}/{1}".format(data_dir, target_dir)
        dir_description = "Raw"
    elif data_format == 6:
        target_dir = link = '/combined'
        if ics:
            data_type = 15
        else:
            data_type = 16
        dl_dir = "{0}/{1}".format(data_dir, target_dir)
        dir_description = "Combined"
    else:
        logger.error("Unable to determine data format from archive. Exiting")
        sys.exit(0)
    mdir(dl_dir, dir_description)
    create_link(data_dir, target_dir, product_dir, link)
    batch_dir = product_dir+"/batch/"

    for time_to_get in range(start_time,stop_time,increment):
        if time_to_get + increment > stop_time:
            increment = stop_time - time_to_get + 1
        get_data = "{0} --obs={1} --type={2} --from={3} --duration={4} --parallel={5} --dir={6}".\
                   format(voltdownload, obsid, data_type, time_to_get,
                          (increment-1), parallel, dl_dir)
        #need to subtract 1 from increment since voltdownload wants how many
        #seconds PAST the first one

        if head:
            log_name="{0}/voltdownload_{1}.log".format(product_dir,time_to_get)
            with open(log_name, 'w') as log:
                subprocess.call(get_data, shell=True, stdout=log, stderr=log)
        else:
            voltdownload_batch = "volt_{0}".format(time_to_get)
            check_batch = "check_volt_{0}".format(time_to_get)
            volt_secs_to_run = datetime.timedelta(seconds=500*increment)
            check_secs_to_run = "15:00"
            if data_type == 16:
                check_secs_to_run = "10:15:00"

            checks = distutils.spawn.find_executable("checks.py")
            # Write out the checks batch file but don't submit it
            commands = []
            if vcs_database_id is not None:
                commands.append(database_vcs.add_database_function())
            commands.append("newcount=0")
            commands.append("let oldcount=$newcount-1")
            commands.append("sed -i -e \"s/oldcount=${{oldcount}}/oldcount=${{newcount}}/\" {0}".\
                            format(batch_dir+voltdownload_batch+".batch"))
            commands.append("oldcount=$newcount; let newcount=$newcount+1")
            commands.append("sed -i -e \"s/_${{oldcount}}.out/_${{newcount}}.out/\" {0}".\
                            format(batch_dir+voltdownload_batch+".batch"))
            checks_command = "-m download -o {0} -w {1} -b {2} -i {3} --data_type {4}".format(obsid,
                              dl_dir, time_to_get, increment, data_type)
            if vcs_database_id is None:
                commands.append('{0} {1}'.format(checks, checks_command))
            else:
                commands.append('run "{0}" "{1}" "{2}"'.format(checks, checks_command, vcs_database_id))
            commands.append("if [ $? -eq 1 ];then")
            commands.append("sbatch {0}".format(batch_dir+voltdownload_batch+".batch"))
            # if we have tarballs we send the untar jobs to the workq
            if data_type == 16:
                commands.append("else")
                untar = distutils.spawn.find_executable('untar.sh')
                untar_command = "-w {0} -o {1} -b {2} -e {3} -j {4} {5}".format(dl_dir,
                                 obsid, time_to_get, time_to_get+increment-1, n_untar,
                                 keep)
                if vcs_database_id is None:
                    commands.append('{0} {1}'.format(untar, untar_command))
                else:
                    commands.append('run "{0}" "{1}" "{2}"'.format(untar,
                                    untar_command, vcs_database_id))

                #commands.append("sbatch {0}.batch".format(batch_dir+tar_batch))
            commands.append("fi")

            # Download and checks should be done on Zeus's cpuq. This will only work
            # on Galaxy as the Ozstar workflow is different
            submit_slurm(check_batch, commands, batch_dir=batch_dir,
                         slurm_kwargs={"time": check_secs_to_run,
                                       "nice": nice},
                         vcstools_version=vcstools_version, submit=False,
                         outfile=batch_dir+check_batch+"_0.out",
                         queue="zcpuq", export="NONE", mem=10240)

            # Write out the tar batch file if in mode 15
            #if format == 16:
            #        body = []
            #        for t in range(time_to_get, time_to_get+increment):
            #                body.append("aprun tar -xf {0}/1149620392_{1}_combined.tar".format(dl_dir,t))
            #        submit_slurm(tar_batch,body,batch_dir=working_dir+"/batch/", slurm_kwargs={"time":"1:00:00", "partition":"gpuq" })


            #module_list=["mwa-voltage/master"]
            #removed the master version load because by default we load the python 3 version
            module_list=[]
            body = []
            if vcs_database_id is not None:
                body.append(database_vcs.add_database_function())
            body.append("oldcount=0")
            body.append("let newcount=$oldcount+1")
            body.append("if [ ${newcount} -gt 10 ]; then")
            body.append("echo \"Tried ten times, this is silly. Aborting here.\";exit")
            body.append("fi")
            body.append("sed -i -e \"s/newcount=${{oldcount}}/newcount=${{newcount}}/\" {0}\n".\
                        format(batch_dir+check_batch+".batch"))
            body.append("sed -i -e \"s/_${{oldcount}}.out/_${{newcount}}.out/\" {0}".\
                        format(batch_dir+check_batch+".batch"))
            body.append("sbatch -d afterany:${{SLURM_JOB_ID}} {0}".\
                        format(batch_dir+check_batch+".batch"))
            voltdownload_command = "--obs={0} --type={1} --from={2} --duration={3} --parallel={4}"\
                                   " --dir={5}".format(obsid, data_type, time_to_get, increment-1,
                                   parallel, dl_dir)
            if vcs_database_id is None:
                body.append("{0} {1}".format(voltdownload, voltdownload_command))
            else:
                body.append('run "{0}" "{1}" "{2}"'.format(voltdownload, voltdownload_command,
                                                           vcs_database_id))
            submit_slurm(voltdownload_batch, body, batch_dir=batch_dir,
                         module_list=module_list,
                         slurm_kwargs={"time": str(volt_secs_to_run),
                                       "nice" : nice},
                         vcstools_version=vcstools_version,
                         outfile=batch_dir+voltdownload_batch+"_1.out",
                         queue="copyq", export="NONE", mem=5120)

            # submit_cmd = subprocess.Popen(volt_submit_line,shell=True,stdout=subprocess.PIPE)
            continue
    # TODO: what is the below doing here???
        try:
            os.chdir(product_dir)
        except:
            logging.error("cannot open working dir:{0}".format(product_dir))
            sys.exit()



def download_cal(obs_id, cal_obs_id, data_dir, product_dir, vcs_database_id, head=False,
                 vcstools_version="master", nice=0):

    batch_dir = product_dir + '/batch/'
    product_dir = '{0}/cal/{1}'.format(product_dir,cal_obs_id)
    mdir(product_dir, 'Calibrator product')
    mdir(batch_dir, 'Batch')
    # obsdownload creates the folder cal_obs_id regardless where it runs
    # this deviates from our initially inteded naming conventions of
    # /astro/mwavcs/vcs/[cal_obs_id]/vis but the renaming and linking is a pain otherwise,
    # hence we'll link vis agains /astro/mwavcs/vcs/[cal_obs_id]/[cal_obs_id]
    #target_dir = '{0}'.format(cal_obs_id)
    link = 'vis'
    csvfile = "{0}{1}_dl.csv".format(batch_dir,cal_obs_id)
    obsdownload = distutils.spawn.find_executable("obsdownload.py")
    #get_data = "{0} -o {1} -d {2}".format(obsdownload,cal_obs_id, data_dir)
    if head:
        logging.error("I'm sorry, this option is no longer supported. Please "
                      "download through the copyq.")
        # print "Will download the data from the archive. This can take a while so please do not ctrl-C."
        # log_name="{0}/caldownload_{1}.log".format(batch_dir,cal_obs_id)
        # with open(log_name, 'w') as log:
        #     subprocess.call(get_data, shell=True, stdout=log, stderr=log)
        # create_link(data_dir, target_dir, product_dir, link)
        # #clean up
        # try:
        #     os.remove('obscrt.crt')
        #     os.remove('obskey.key')
        # except:
        #     pass
    else:
        # we create the link using bash and not our create_link function
        # as we'd like to do this only once the data have arrived,
        # i.e. the download worked.
        make_link = "ln -sfn {0} {1}/{2}".format(data_dir, product_dir, link)
        obsdownload_batch = "caldownload_{0}".format(cal_obs_id)
        secs_to_run = "03:00:00" # sometimes the staging can take a while...
        module_list = ["setuptools"]
        commands = []
        commands.append("module load manta-ray-client")
        if vcs_database_id is not None:
            commands.append(database_vcs.add_database_function())
        commands.append("csvfile={0}".format(csvfile))
        commands.append('cd {0}'.format(data_dir))
        commands.append('if [[ -z ${MWA_ASVO_API_KEY} ]]')
        commands.append('then')
        commands.append('    echo "Error, MWA_ASVO_API_KEY not set"')
        commands.append('    echo "Cannot use client"')
        commands.append('    echo "Please read the MWA ASVO documentation '
                        'about setting this (https://wiki.mwatelescope.org/'
                        'display/MP/MWA+ASVO%3A+Release+Notes)"')
        commands.append('    exit 1')
        commands.append('fi')
        commands.append('echo "obs_id={0}, job_type=d, download_type=vis" > {1}'.\
                        format(cal_obs_id,csvfile))
        commands.append('mwa_client --csv={0} --dir={1}'.format(csvfile,data_dir))
        commands.append(make_link)
        commands.append('unzip *.zip')
        submit_slurm(obsdownload_batch, commands, batch_dir=batch_dir,
                     module_list=module_list,
                     slurm_kwargs={"time": secs_to_run, "nice": nice},
                     vcstools_version=vcstools_version, queue="copyq",
                     export="NONE", mem=4096)


def vcs_recombine(obsid, start_time, stop_time, increment, data_dir, product_dir,
                  vcs_database_id, vcstools_version="master", nice=0):

    #Load computer dependant config file

    logger.info("Running recombine on files")
    jobs_per_node = 8
    target_dir = link = 'combined'
    mdir(data_dir + '/' + target_dir, 'Combined')
    create_link(data_dir, target_dir, product_dir, link)
    batch_dir = product_dir+"/batch/"
    recombine = distutils.spawn.find_executable("recombine.py")
    checks = distutils.spawn.find_executable("checks.py")
    recombine_binary = distutils.spawn.find_executable("recombine")
    for time_to_get in range(start_time,stop_time,increment):
        process_nsecs = increment if (time_to_get + increment <= stop_time) \
                                  else (stop_time - time_to_get + 1)
        if (jobs_per_node > process_nsecs):
                jobs_per_node = process_nsecs
        nodes = (increment+(-increment%jobs_per_node))//jobs_per_node + 1 # Integer division with ceiling result plus 1 for master node
        recombine_batch = "recombine_{0}".format(time_to_get)
        check_batch = "check_recombine_{0}".format(time_to_get)
        module_list = ["module switch PrgEnv-cray PrgEnv-gnu", "python/3.6.3", "numpy/1.13.3", "mwa-voltage/master"]
        commands = []
        if vcs_database_id is not None:
            commands.append(database_vcs.add_database_function())
        commands.append("newcount=0")
        commands.append("let oldcount=$newcount-1")
        commands.append("sed -i -e \"s/oldcount=${{oldcount}}/oldcount=${{newcount}}/\" {0}".\
                        format(batch_dir+recombine_batch+".batch"))
        commands.append("oldcount=$newcount; let newcount=$newcount+1")
        commands.append("sed -i -e \"s/_${{oldcount}}.out/_${{newcount}}.out/\" {0}".\
                        format(batch_dir+recombine_batch+".batch"))
        checks_command = "-m recombine -o {0} -w {1}/combined/ -b {2} -i {3}".format(obsid,
                          data_dir, time_to_get, process_nsecs)
        if vcs_database_id is None:
            commands.append("{0} {1}".format(checks, checks_command))
        else:
            commands.append('run "{0}" "{1}" "{2}"'.format(checks, checks_command,
                                                           vcs_database_id))
        commands.append("if [ $? -eq 1 ];then")
        commands.append("sbatch {0}".format(batch_dir+recombine_batch+".batch"))
        commands.append("fi")
        submit_slurm(check_batch, commands, batch_dir=batch_dir,
                     module_list=module_list,
                     slurm_kwargs={"time": "15:00", "nice": nice},
                     vcstools_version=vcstools_version, submit=False,
                     outfile=batch_dir+check_batch+"_0.out",
                     queue='gpuq', export="NONE")

        module_list = ["module switch PrgEnv-cray PrgEnv-gnu", "python/3.6.3",
                       "numpy/1.13.3", "mwa-voltage/master", "mpi4py", "cfitsio"]
        commands = []
        if vcs_database_id is not None:
            commands.append(database_vcs.add_database_function())
        commands.append("oldcount=0")
        commands.append("let newcount=$oldcount+1")
        commands.append("if [ ${newcount} -gt 10 ]; then")
        commands.append("echo \"Tried ten times, this is silly. Aborting here.\";exit")
        commands.append("fi")
        commands.append("sed -i -e \"s/newcount=${{oldcount}}/newcount=${{newcount}}/\" {0}".\
                        format(batch_dir+check_batch+".batch"))
        commands.append("sed -i -e \"s/_${{oldcount}}.out/_${{newcount}}.out/\" {0}".\
                        format(batch_dir+check_batch+".batch"))
        commands.append("sbatch -d afterany:${{SLURM_JOB_ID}} {0}".\
                        format(batch_dir+check_batch+".batch"))
        recombine_command = "-o {0} -s {1} -w {2} -e {3}".format(obsid, time_to_get,
                            data_dir, recombine_binary)
        if vcs_database_id is None:
            commands.append("srun --export=all python3 {0} {1}".format(recombine, recombine_command))
        else:
            commands.append('run "srun --export=all python3 {0}" "{1}" "{2}"'.format(recombine,
                            recombine_command, vcs_database_id))

        submit_slurm(recombine_batch, commands, batch_dir=batch_dir,
                     module_list=module_list,
                     slurm_kwargs={"time": "06:00:00", "nodes": str(nodes),
                                   "ntasks-per-node": jobs_per_node,
                                   "nice": nice},
                     vcstools_version=vcstools_version,
                     outfile=batch_dir+recombine_batch+"_1.out",
                     queue='gpuq', export="NONE")




def vcs_correlate(obsid,start,stop,increment, data_dir, product_dir, ft_res,
                  vcs_database_id, metafits, vcstools_version="master", nice=0):

    logger.info("Correlating files at {0} kHz and {1} milliseconds".\
                format(ft_res[0], ft_res[1]))

    batch_dir = product_dir+"/batch/"
    target_dir = link = 'vis'

    if data_dir == product_dir:
        corr_dir = "{0}/cal/{1}/{2}".format(product_dir, obsid, target_dir)
    else:
        corr_dir = "{0}/{1}".format(data_dir, target_dir)
        product_dir = "{0}/cal/{1}/".format(product_dir, obsid)
        mdir(product_dir, "Correlator")
    mdir(corr_dir, "Correlator Product")
    create_link(data_dir, target_dir, product_dir, link)

    chan_list = get_frequencies(metafits_file, resort=True)
    #gpu_int = 0.01 # Code was compiled with a hard-coded 100 sample minimum intigration. For 'normal' data this means 0.01 seconds
    gpu_int = 10 # Code was compiled with a hard-coded 100 sample minimum integration. For 'normal' data this means 10 milliseconds.
    integrations=int(ft_res[1]/gpu_int)
    #num_frames=int(1.0/ft_res[1])
    num_frames=int(1000/ft_res[1])

    logger.info("Input chan list is {0}".format(chan_list))

    for time_to_get in range(start,stop,increment):
        inc_start = time_to_get
        inc_stop = time_to_get+increment
        for index,channel in enumerate(chan_list):
            gpubox_label = (index+1)
            f=[]
            for time_to_corr in range(inc_start,inc_stop,1):
                file_to_process = "{0}/combined/{1}_{2}_ch{3:0>2}.dat".\
                                  format(data_dir,obsid,time_to_corr,channel)
                #check the file exists
                if (os.path.isfile(file_to_process) == True):
                    f.append(file_to_process)

            #now have a full list of files
            #for this increment
            #and this channel
            if (len(f) > 0):
                corr_batch = "correlator_{0}_gpubox{1:0>2}".format(inc_start,gpubox_label)
                body = []
                if vcs_database_id is not None:
                    body.append(database_vcs.add_database_function())
                to_corr = 0
                for file in f:
                    (current_time,_) = os.path.splitext(os.path.basename(file))
                    (obsid,gpstime,_) = current_time.split('_')
                    t = Time(int(gpstime), format='gps', scale='utc')
                    unix_time = int(t.unix)

                    offline_correlator_command = "-o {0}/{1} -s {2} -r {3} -i {4} -f 128 -n {5} "\
                                                 "-c {6:0>2} -d {7}".format(corr_dir, obsid,
                                                 unix_time, num_frames, integrations,
                                                 int(ft_res[0]/10), gpubox_label, file)
                    if vcs_database_id is None:
                        body.append("{0} {1}".format("offline_correlator",
                                    offline_correlator_command))
                    else:
                        body.append('run "{0}" "{1}" "{2}"'.format("offline_correlator",
                                    offline_correlator_command, vcs_database_id))
                    to_corr += 1

                module_list = ["module switch PrgEnv-cray PrgEnv-gnu"]
                secs_to_run = str(datetime.timedelta(seconds=2*12*num_frames*to_corr))
                # added factor two on 10 April 2017 as galaxy seemed really slow...
                submit_slurm(corr_batch, body, module_list=module_list,
                             slurm_kwargs={"time": secs_to_run, "nice": nice,
                                           "gres": "gpu:1"},
                             queue='gpuq', vcstools_version=vcstools_version,
                             batch_dir=batch_dir, export="NONE")
            else:
                logger.error("Couldn't find any recombine files. Aborting here.")


def coherent_beam(obs_id, start, stop, data_dir, product_dir, batch_dir,
                  metafits_file, nfine_chan, pointing_list,
                  rts_flag_file=None, bf_formats=None, DI_dir=None,
                  execpath=None, calibration_type='rts', ipfb_filter="LSQ12",
                  vcstools_version="master", nice=0, channels_to_beamform=None,
                  dpp=False):
    """
    This function runs the new version of the beamformer. It is modelled after
    the old function above and will likely be able to be streamlined after
    working implementation (SET)

    Streamlining underway, as well as full replacement of the old function (SET March 28, 2018)
    """

    #Load computer dependant config file
    comp_config = load_config_file()

    # If execpath is given, change the make_beam executable command
    # otherwise, it should be on the PATH if vcstools has been installed
    if execpath:
        make_beam_cmd = "{0}/make_beam".format(execpath)
        make_beam_version_cmd = "{0}/make_beam -V".format(execpath)
    else:
        make_beam_cmd = "make_beam"
        make_beam_version_cmd = "make_beam -V"

    make_beam_version = subprocess.Popen(make_beam_version_cmd,
                           stdout=subprocess.PIPE, shell=True).communicate()[0]
    #tested_version = "?.?.?"
    logger.info("Current version of make_beam = {0}".format(make_beam_version.strip()))
    #logger.info("Tested version of make_beam = {0}".format(tested_version.strip()))

    metafile = "{0}/{1}.meta".format(product_dir, obs_id)
    channels = None
    # No channels given so first check for a metafile
    if os.path.isfile(metafile):
        logger.info("Found observation metafile: {0}".format(metafile))
        with open(metafile, 'r') as m:
            for line in m.readlines():
                if line.startswith("channels"):
                    channels = line.split(",")[1:]
                    channels = np.array(channels, dtype=np.int)
    else:
        logger.debug("No metafile in {0}".format(metafile))
    logger.debug("Channels before meta.get_channels: {0}".format(channels))
    # If channels is still None get_channels will get it from the metadata
    channels = meta.get_channels(obs_id, channels=channels)

    # Make a metafile containing the channels so no future metadata calls are required
    if not os.path.isfile(metafile):
        with open(metafile, "w") as m:
            m.write("#Metadata for obs ID {0} required to determine if: normal or "
                    "picket-fence\n".format(obs_id))
            m.write("channels,{0}".format(",".join([str(c) for c in channels])))
    channels = np.array(channels, dtype=np.int)
    hichans = [c for c in channels if c>128]
    lochans = [c for c in channels if c<=128]
    lochans.extend(list(reversed(hichans)))
    ordered_channels = lochans

    if channels_to_beamform is None:
        # If no channels_to_beamform given fold on everything
        channels_to_beamform = ordered_channels

    # Run for each coarse channel. Calculates delays and makes beam

    if not DI_dir:
        logger.error("You need to specify the path to the calibrator files, "
                     "either where the DIJs are or where the Offringa "
                     "calibration_solutions.bin file is. Aborting here")
        sys.exit(0)
    DI_dir = os.path.abspath(DI_dir)

    # make_beam_small requires the start time in UTC, get it from the start
    utctime = gps_to_utc(start)

    if dpp:
        # running the Data Processing Pipeline so make a symlink to /astro so we
        # don't run out of space on /group
        target_dir = link = "dpp_pointings"
        P_dir = os.path.join(product_dir, link)
        mdir(os.path.join(data_dir, target_dir), "DPP Pointings")
        create_link(data_dir, target_dir, product_dir, link)
    else:
        P_dir = os.path.join(product_dir, "pointings")
        mdir(P_dir, "Pointings")
    # startjobs = True

    # Set up supercomputer dependant parameters
    import socket
    hostname = socket.gethostname()
    if hostname.startswith('john') or hostname.startswith('farnarkle'):
        max_pointing = 120
        #Work out required SSD size
        temp_mem = int(0.0012 * (float(stop) - float(start) + 1.) * \
                       float(len(pointing_list)) ) + 1
        # Split it up into 400 chuncks to not use more than 60BG
        #time_chunks = gps_time_lists(start, stop, 400)
        time_chunks = gps_time_lists(start, stop, 10000)
    else:
        max_pointing = 15
        temp_mem = None
        # Do it all at once on Galaxy
        time_chunks = gps_time_lists(start, stop, 10000)

    # set up SLURM requirements
    if len(pointing_list) > max_pointing:
        seconds_to_run = 8 * (stop - start + 1) * max_pointing
    else:
        seconds_to_run = 8 * (stop - start + 1) * len(pointing_list)

    if seconds_to_run > 86399.:
        secs_to_run = datetime.timedelta(seconds=86399)
    else:
        secs_to_run = datetime.timedelta(seconds=seconds_to_run)

    # Get the project id (eg G0057) from the metafits file
    with pyfits.open(metafits_file) as hdul:
        project_id = hdul[0].header['project']

    # splits the pointing list into lists of length max_pointing
    pointing_list_list = list(chunks(pointing_list, max_pointing))
    time_now = str(datetime.datetime.now()).replace(" ", "_")

    logging.info("Running make_beam")
    job_id_list_list = []
    for pl, pointing_list in enumerate(pointing_list_list):
        pointing_str = ",".join(pointing_list)
        # Run one coarse channel per node
        job_id_list = []
        for gpubox, coarse_chan in enumerate(ordered_channels, 1):
            if coarse_chan not in channels_to_beamform:
                continue
            if calibration_type == 'rts':
                #chan_list = get_frequencies(metafits_file, resort=True)
                DI_file = "{0}/DI_JonesMatrices_node{1:0>3}.dat".format(DI_dir, gpubox)
                jones_option = "-J {0}".format(DI_file)
            elif calibration_type == 'offringa':
                #chan_list = get_frequencies(metafits_file, resort=False)
                DI_file = "{0}/calibration_solution.bin".format(DI_dir)
                jones_option = "-O {0} -C {1}".format(DI_file, int(gpubox) - 1)
            else:
                logger.info("Please an accepted calibratin type. Aborting here.")
                sys.exit(0)

            # Making pointing directories
            for pointing in pointing_list:
                mdir("{0}/{1}".format(P_dir, pointing), "Pointing {0}".format(pointing))

            n_omp_threads = 1
            if "v" in bf_formats:
                for pointing in pointing_list:
                    make_beam_small_batch = "mb_{0}_ch{1}".format(pointing, coarse_chan)
                    module_list = [comp_config['container_module']]
                    commands = []
                    commands.append("cd {0}/{1}".format(P_dir,pointing))
                    runline = "srun --export=all -n 1"
                    runline += " -c {}".format(n_omp_threads)
                    runline += " {}".format(comp_config['container_command'])
                    runline += " {}".format(make_beam_cmd)
                    runline += " -o {}".format(obs_id)
                    runline += " -b {}".format(start)
                    runline += " -e {}".format(stop)
                    runline += " -a 128"
                    runline += " -n 128"
                    runline += " -f {}".format(coarse_chan)
                    runline += " {}".format(jones_option)
                    runline += " -d {}/combined".format(data_dir)
                    runline += " -P {}".format(pointing)
                    runline += " -r 10000"
                    runline += " -m {}".format(metafits_file)
                    runline += " -z {}".format(utctime)
                    runline += " {}".format(bf_formats)
                    runline += " -F {}".format(rts_flag_file)
                    runline += " -S {}".format(ipfb_filter)
                    commands.append(runline)
        
                    job_id = submit_slurm(make_beam_small_batch, commands,
                                batch_dir=batch_dir, module_list=module_list,
                                slurm_kwargs={"time":secs_to_run, "nice":nice},
                                queue='gpuq', vcstools_version=vcstools_version,#forces olf version with vdif
                                submit=True, export="NONE", gpu_res=1,
                                cpu_threads=n_omp_threads,
                                mem=comp_config['gpu_beamform_mem'])
                    job_id_list.append(job_id)

            else:
                make_beam_small_batch = "mb_{0}_{1}_ch{2}".format(pl, time_now, coarse_chan)
                module_list = [comp_config['container_module']]
                commands = []
                if hostname.startswith('john') or hostname.startswith('farnarkle'):
                    # Write outputs to SSDs if on Ozstar
                    commands.append("cd $JOBFS")
                else:
                    commands.append("cd {0}".format(P_dir))

                # Loop over each GPS time chunk. Will only be one on Galaxy
                for tci, (start, stop) in enumerate(time_chunks):
                    utctime = gps_to_utc(start)
                    commands.append("srun --export=all -n 1 -c {0} {1} {2} -o {3} -b {4} "
                                    "-e {5} -a 128 -n 128 -f {6} {7} -d {8}/combined "
                                    "-P {9} -r 10000 -m {10} -z {11} {12} -F {13}".format(
                                    n_omp_threads, comp_config['container_command'],
                                    make_beam_cmd, obs_id, start, stop, coarse_chan,
                                    jones_option, data_dir, pointing_str, metafits_file,
                                    utctime, bf_formats, rts_flag_file))
                    commands.append("")
                    if hostname.startswith('john') or hostname.startswith('farnarkle'):
                        for pointing in pointing_list:
                            """
                            # Move outputs off the SSD and renames them as if makebeam
                            # was only run once
                            # G0024_1166459712_07:42:49.00_-28:21:43.00_ch132_0001.fits
                            commands.append("cp $JOBFS/{0}/{1}_{2}_{0}_ch{3}_0001.fits "
                                    "{4}/{0}/{1}_{2}_{0}_ch{3}_00{5:02}.fits".format(pointing,
                                    project_id, obs_id, coarse_chan, P_dir, (tci+1)*2-1))
                            commands.append("cp $JOBFS/{0}/{1}_{2}_{0}_ch{3}_0002.fits "
                                    "{4}/{0}/{1}_{2}_{0}_ch{3}_00{5:02}.fits".format(pointing,
                                    project_id, obs_id, coarse_chan, P_dir, (tci+1)*2))
                            """
                            commands.append("cp $JOBFS/{0}/{1}_{2}_{0}_ch{3}_00*.fits "
                                            "{4}/{0}/".format(pointing, project_id,
                                                              obs_id, coarse_chan, P_dir))
                    commands.append("")

                job_id = submit_slurm(make_beam_small_batch, commands,
                            batch_dir=batch_dir, module_list=module_list,
                            slurm_kwargs={"time":secs_to_run, "nice":nice},
                            queue='gpuq', vcstools_version=vcstools_version,
                            submit=True, export="NONE", gpu_res=1,
                            cpu_threads=n_omp_threads,
                            mem=comp_config['gpu_beamform_mem'], temp_mem=temp_mem)
                job_id_list.append(job_id)
        job_id_list_list.append(job_id_list)

    return job_id_list_list, make_beam_small_batch.split('ch')[0]

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

if __name__ == '__main__':

    # Dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING)

    modes=['download', 'download_ics', 'download_cal', 'recombine','correlate', 'beamform']
    bf_out_modes=['psrfits', 'vdif', 'both']
    jobs_per_node = 8
    chan_list_full=["ch01","ch02","ch03","ch04","ch05","ch06","ch07","ch08",
                    "ch09","ch10","ch11","ch12","ch13","ch14","ch15","ch16",
                    "ch17","ch18","ch19","ch20","ch21","ch22","ch23","ch24"]
    chan_list = []
    jones = "-j jones.txt"

    import argparse

 #   parser=OptionParser(description="process_vcs.py is a script of scripts that downloads prepares and submits jobs to Galaxy. It can be run with just a pointing (-p \"xx:xx:xx xx:xx:xx.x\") and an obsid (\"-o <obsid>\") and it will process all the data in the obs. It will call prepare.py which will attempt to build the phase and calibration information - which will only exist if a calibration obs has already been run. So it will only move past the data prepa stage if the \"-r\" (for run) is used\n"

    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            description="process_vcs.py is a script for processing the MWA VCS data on Galaxy in steps. It can download data from the archive, call on recombine to form course channels, run the offline correlator, make tile re-ordered and bit promoted PFB files or for a coherent beam for a given pointing.")
    group_download = parser.add_argument_group('Download Options')
    group_download.add_argument("--head", action="store_true", default=False, help="Submit download jobs to the headnode instead of the copyqueue ")
    #group_download.add_argument("--format", type="choice", choices=['11','15','16'], default='11', help="Voltage data type (Raw = 11, ICS Only = 15, Recombined and ICS = 16) ")
    group_download.add_argument("-d", "--parallel_dl", type=int, default=3, help="Number of parallel downloads to envoke ")
    group_download.add_argument("-j", "--untar_jobs", type=int, default=2, help="Number of parallel jobs when untaring downloaded tarballs. ")
    group_download.add_argument("-k", "--keep_tarball", action="store_true", default=False, help="Keep the tarballs after unpacking. ")
    group_correlate = parser.add_argument_group('Correlator Options')
    group_correlate.add_argument("--ft_res", metavar="FREQ RES,TIME RES", type=int, nargs=2, default=(10,1000), help="Frequency (kHz) and Time (ms) resolution for running the correlator. Please make divisible by 10 kHz and 10 ms respectively. ")

    group_beamform = parser.add_argument_group('Beamforming Options')
    group_beamform.add_argument("-p", "--pointings", type=str, nargs='*', help="A space sepertated list of pointings with the RA and Dec seperated by _ in the format HH:MM:SS_+DD:MM:SS, e.g. \"19:23:48.53_-20:31:52.95 19:23:40.00_-20:31:50.00\"")
    group_beamform.add_argument("--pointing_file", help="A file containing pointings with the RA and Dec seperated by _ in the format HH:MM:SS_+DD:MM:SS on each line, e.g. \"19:23:48.53_-20:31:52.95\n19:23:40.00_-20:31:50.00\"")
    group_beamform.add_argument("--DI_dir", default=None, help="Directory containing either Direction Independent Jones Matrices (as created by the RTS) or calibration_solution.bin as created by Andre Offringa's tools.[no default]")
    group_beamform.add_argument("--bf_out_format", type=str, choices=['psrfits','vdif','both'], help="Beam former output format. Choices are {0}. ".format(bf_out_modes), default='psrfits')
    group_beamform.add_argument("--incoh", action="store_true", default=False, help="Add this flag if you want to form an incoherent sum as well. ")
    group_beamform.add_argument("--sum", action="store_true", default=False, help="Add this flag if you the beamformer output as summed polarisations (only Stokes I). This reduces the data size by a factor of 4.")
    group_beamform.add_argument("--flagged_tiles", type=str, default=None, help="Path (including file name) to file containing the flagged tiles as used in the RTS, will be used by get_delays. ")
    group_beamform.add_argument('--cal_type', type=str, help="Use either RTS (\"rts\") solutions or Andre-Offringa-style (\"offringa\") solutions. Default is \"rts\". If using Offringa's tools, the filename of calibration solution must be \"calibration_solution.bin\".", default="rts")
    group_beamform.add_argument("-E", "--execpath", type=str, default=None, help="Supply a path into this option if you explicitly want to run files from a different location for testing. Default is None (i.e. whatever is on your PATH).")
    group_beamform.add_argument("--ipfb_filter", type=str, choices=['LSQ12','MIRROR'], help="The filter to use when performing the inverse PFB", default='LSQ12')
   

    parser.add_argument("-m", "--mode", type=str, choices=['download','download_ics', 'download_cal', 'recombine','correlate', 'beamform'], help="Mode you want to run. {0}".format(modes))
    parser.add_argument("-o", "--obs", metavar="OBS ID", type=int, help="Observation ID you want to process [no default]")
    parser.add_argument('--cal_obs', '-O', metavar="CALIBRATOR OBS ID", type=int, help="Only required in 'download_cal' mode."
                          "Observation ID of calibrator you want to process. In case of "
                          "in-beam calibration should be the same as input to -o (obsID). [no default]", default=None)
    parser.add_argument("-b", "--begin", type=int, help="First GPS time to process [no default]")
    parser.add_argument("-e", "--end", type=int, help="Last GPS time to process [no default]")
    parser.add_argument("-a", "--all", action="store_true", default=False, help="Perform on entire observation span. Use instead of -b & -e. ")
    parser.add_argument("--all_avail", action="store_true", default=False,
                      help="Uses all of the available combined files available (does not check for gaps). [default=False]")
    parser.add_argument("-i", "--increment", type=int, default=64, help="Increment in seconds (how much we process at once) ")
    parser.add_argument("-s", action="store_true", default=False, help="Single step (only process one increment and this is it (False == do them all) ")
    parser.add_argument("-w", "--work_dir", metavar="DIR", default=None, help="Base directory you want run things in. USE WITH CAUTION! Per default "
                          "raw data will will be downloaded into /astro/mwavcs/vcs/[obsID] and data products will be in /group/mwavcs/vcs/[obsID]."
                          " If set, this will create a folder for the Obs. ID if it doesn't exist ")
    parser.add_argument("-c", "--ncoarse_chan", type=int, default=24, help="Coarse channel count (how many to process) ")
    parser.add_argument("-n", "--nfine_chan", type=int, default=128, help="Number of fine channels per coarse channel ")
    parser.add_argument("--mail",action="store_true", default=False, help="Enables e-mail notification about start, end, and fail of jobs. Currently only implemented for beamformer mode.")
    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO",
                                        default="INFO")
    parser.add_argument("-V", "--version", action="store_true", help="Print version and quit")
    parser.add_argument("--vcstools_version", type=str, default="master", help="VCSTools version to load in jobs (i.e. on the queues) ")
    parser.add_argument("--nice", type=int, default=0, help="Reduces your priority of Slurm Jobs. ")

    args = parser.parse_args()

    # set up the logger for stand-alone execution
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.propagate = False

    logger.info("Using vcstools/{0}".format(args.vcstools_version))
    if args.version:
        try:
            import version
            logger.info(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            logger.error("Couldn't import version.py - have you installed vcstools?")
            logger.error("ImportError: {0}".format(ie))
            sys.exit(0)

    #Option parsing
    if not args.obs:
        logger.error("Observation ID required, please put in with -o or --obs")
        sys.exit(0)
    if args.all and (args.begin or args.end):
        logger.error("Please specify EITHER (-b,-e) OR -a")
        sys.exit(0)
    elif args.all:
        args.begin, args.end = meta.obs_max_min(args.cal_obs\
                               if args.mode == 'download_cal' else args.obs)
    elif args.mode != 'download_cal' and (not args.begin or not args.end):
        logger.debug(args.mode)
        logger.error("Please specify EITHER (-b,-e) OR -a")
        sys.exit(0)
    # make sure we can process increments smaller than 64 seconds when not in calibration related mode
    if args.mode != 'download_cal':
        if args.end - args.begin +1 < args.increment:
            args.increment = args.end - args.begin + 1
    e_mail = ""
    if args.mail:
        e_mail = get_user_email()
        logger.info("Sending info to {0}".format(e_mail))
    if not args.mode:
      logger.error("Mode required {0}. Please specify with -m or --mode.".format(modes))

      sys.exit(0)
    if args.begin and args.end:
        if args.begin > args.end:
            logger.error("Starting time is after end time")
            sys.exit(0)
    if (args.mode == "beamform" or args.incoh):
        bf_format = ""
        if not (args.pointings or args.pointing_file):
            logger.info("Beamformer mode required that you specify pointings "
                        "using either -p or --pointing_file.")
            sys.exit(0)
        #check if they're using more than one of the pointing options
        if (args.pointings and args.pointing_file):
            logger.info("Beamformer mode requires only one pointing option. "
                        "Please use either -p or --pointing_file.")
            sys.exit(0)
        if (args.bf_out_format == 'psrfits' or args.bf_out_format == 'both'):
            bf_format +=" -p"
            logger.info("Writing out PSRFITS.")
        if  (args.bf_out_format == 'vdif' or args.bf_out_format == 'both'):
            bf_format += " -v"
            logger.info("Writing out upsampled VDIF.")
        if (args.incoh):
            bf_format += " -i"
            logger.info("Writing out incoherent sum.")
        if (args.sum):
            bf_format += " -s"

        # This isn't necessary as checks for execpath are done in beamforming function (BWM 6/4/18)
        #if args.execpath:
        #    execpath = args.execpath

    #Load computer dependant config file
    comp_config = load_config_file()

    if args.work_dir:
        logger.warning("YOU ARE MESSING WITH THE DEFAULT DIRECTORY STRUCTURE "
                       "FOR PROCESSING -- BE SURE YOU KNOW WHAT YOU ARE DOING!")
        sleep(5)
        data_dir = product_dir = "{0}/{1}".format(args.work_dir, args.obs)
    else:
        data_dir = '{0}{1}'.format(comp_config['base_data_dir'],args.obs)
        product_dir = '{0}{1}'.format(comp_config['base_product_dir'],args.obs)
    batch_dir = "{0}/batch".format(product_dir)
    mdir(data_dir, "Data")
    mdir(product_dir, "Products")
    mdir(batch_dir, "Batch")
    metafits_file = "{0}/{1}_metafits_ppds.fits".format(data_dir, args.obs)
    # TODO: modify metafits downloader to not just do a trivial wget

    logger.info("Processing Obs ID {0} from GPS times {1} till {2}".\
                format(args.obs, args.begin, args.end))

    # Record command in database
    if DB_FILE is not None:
        vcs_database_id = database_vcs.database_command(sys.argv, args.obs)
    else:
        vcs_database_id = None

    if args.mode == 'download_ics':
        logger.info("Mode: {0}".format(args.mode))
        vcs_download(args.obs, args.begin, args.end, args.increment, args.head,
                     data_dir, product_dir, args.parallel_dl, vcs_database_id,
                     ics=True, vcstools_version=args.vcstools_version,
                     nice=args.nice)
    elif args.mode == 'download':
        logger.info("Mode: {0}".format(args.mode))
        vcs_download(args.obs, args.begin, args.end, args.increment, args.head,
                     data_dir, product_dir, args.parallel_dl, vcs_database_id,
                     n_untar=args.untar_jobs,
                     keep='-k' if args.keep_tarball else "",
                     vcstools_version=args.vcstools_version, nice=args.nice)
    elif args.mode == 'recombine':
        logger.info("Mode: {0}".format(args.mode))
        ensure_metafits(data_dir, args.obs, metafits_file)
        vcs_recombine(args.obs, args.begin, args.end, args.increment, data_dir,
                      product_dir, vcs_database_id,
                      vcstools_version=args.vcstools_version, nice=args.nice)
    elif args.mode == 'correlate':
        logger.info("Mode: {0}".format(args.mode))
        ensure_metafits(data_dir, args.obs, metafits_file)
        vcs_correlate(args.obs, args.begin, args.end, args.increment, data_dir,
                      product_dir, args.ft_res, vcs_database_id, metafits_file,
                      vcstools_version=args.vcstools_version, nice=args.nice)
    elif args.mode == 'download_cal':
        logger.info("Mode: {0}".format(args.mode))
        if not args.cal_obs:
            logger.error("You need to also pass the calibrator observation ID."
                         " Aborting here.")
            sys.exit(0)
        if args.cal_obs == args.obs:
            logging.error("The calibrator obsID cannot be the same as the "
                          "target obsID -- there are not gpubox files for "
                          "VCS data on the archive.")
            sys.exit(0)
        data_dir = data_dir.replace(str(args.obs), str(args.cal_obs))
        mdir(data_dir, "Calibrator Data")
        download_cal(args.obs, args.cal_obs, data_dir, product_dir, vcs_database_id,
                     args.head, vcstools_version=args.vcstools_version,
                     nice=args.nice)
    elif args.mode == ('beamform' or 'incoh'):
        logger.info("Mode: {0}".format(args.mode))
        if not args.DI_dir:
            logger.error("You need to specify the path to either where the "
                         "DIJs are or where the offringe calibration_solution.bin"
                         "file is. Aborting here.")
            sys.exit(0)
        if args.flagged_tiles:
            flagged_tiles_file = os.path.abspath(args.flagged_tiles)
            if not os.path.isfile(args.flagged_tiles):
                logger.error("Your are not pointing at a file with your input "
                             "to --flagged_tiles. Aborting here as the "
                             "beamformer will not run...")
                sys.exit(0)
        else:
            if os.path.isfile("{0}/flagged_tiles.txt".format(args.DI_dir)):
                flagged_tiles_file = "{0}/flagged_tiles.txt".format(args.DI_dir)
                logger.info("Found tiles flags in {0}/flagged_tiles.txt. "
                            "Using it by default".format(args.DI_dir))
            else:
                flagged_tiles_file = None
        ensure_metafits(data_dir, args.obs, metafits_file)
        #Turn the pointings into a list
        if args.pointings:
            pointing_list = args.pointings
        elif args.pointing_file:
            with open(args.pointing_file) as f:
                pointing_list = f.readlines()
        else:
            logger.error("Please use either --pointing, --pointing_list or "
                         "--pointing_file when beamforming. Exiting here.")
            sys.exit(0)
        logger.debug(pointing_list)

        coherent_beam(args.obs, args.begin, args.end,
                      data_dir, product_dir, batch_dir,
                      metafits_file, args.nfine_chan, pointing_list,
                      rts_flag_file=flagged_tiles_file,
                      bf_formats=bf_format, DI_dir=args.DI_dir,
                      calibration_type=args.cal_type,
                      vcstools_version=args.vcstools_version, nice=args.nice,
                      execpath=args.execpath, ipfb_filter=args.ipfb_filter)
    else:
        logger.error("Somehow your non-standard mode snuck through. "
                     "Try again with one of {0}".format(modes))
        sys.exit(0)
