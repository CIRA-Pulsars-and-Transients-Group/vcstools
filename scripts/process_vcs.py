#!/usr/bin/env python3


import subprocess
import os
import sys
import glob
import tempfile
import atexit
import hashlib
import datetime
from time import sleep, strptime, strftime
import distutils.spawn
import sqlite3 as lite
from astropy.io import fits as pyfits
from astropy.time import Time
from reorder_chans import *
from mdir import mdir
import logging
import numpy as np

#vcstools functions
from job_submit import submit_slurm
import mwa_metadb_utils as meta
import database_vcs
import config

logger = logging.getLogger(__name__)

def tmp(suffix=".sh"):
    t = tempfile.mktemp(suffix=suffix)
    atexit.register(os.unlink, t)
    return t


def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def get_user_email():
    command="echo `ldapsearch -x \"uid=$USER\" mail |grep \"^mail\"|cut -f2 -d' '`"
    email = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True).communicate()[0]
    return email.strip()


def ensure_metafits(data_dir, obs_id, metafits_file):
    # TODO: To get the actual ppds file should do this with obsdownload -o <obsID> -m

    if not os.path.exists(metafits_file):
        logger.warning("{0} does not exists".format(metafits_file))
        logger.warning("Will download it from the archive. This can take a "+\
                      "while so please do not ctrl-C.")
        logger.warning("At the moment, even through the downloaded file is "+\
                       "labelled as a ppd file this is not true.")
        logger.warning("This is hopefully a temporary measure.")
        #obsdownload = distutils.spawn.find_executable("obsdownload.py")

        get_metafits = "wget http://ws.mwatelescope.org/metadata/fits?obs_id={0} -O {1}".format(obs_id, metafits_file)
        try:
            subprocess.call(get_metafits,shell=True)
        except:
            logger.error("Couldn't download {0}. Aborting.".\
                          format(os.basename(metafits_file)))
            quit()
        # clean up
        #os.remove('obscrt.crt')
        #os.remove('obskey.key')
    # make a copy of the file in the product_dir if that directory exists
    # if it doesn't we might have downloaded the metafits file of a calibrator (obs_id only exists on /astro)
    # in case --work_dir was specified in process_vcs call product_dir and data_dir
    # are the same and thus we will not perform the copy
    product_dir = data_dir.replace('/astro/mwaops/vcs/', '/group/mwaops/vcs/') # being pedantic
    if os.path.exists(product_dir) and not os.path.exists(metafits_file):
        logger.info("Copying {0} to {1}".format(metafits_file, product_dir))
        from shutil import copy2
        copy2("{0}".format(metafits_file), "{0}".format(product_dir))

def create_link(data_dir, target_dir, product_dir, link):
    data_dir = os.path.abspath(data_dir)
    product_dir = os.path.abspath(product_dir)
    if data_dir == product_dir:
        return
    link = link.replace(product_dir, '') # just in case...
    link = product_dir + '/' + link
    target_dir = target_dir.replace(data_dir,'')
    if target_dir.startswith("/"):
        target_dir = target_dir[1:]
    if data_dir.endswith("/"):
        data_dir = data_dir[:-1]
    target_dir = data_dir + '/' + target_dir
    # check if link exists and whether it already points to where we'd like it to
    if os.path.exists(link):
        if os.path.islink(link):
            if os.readlink(link) == target_dir:
                return
            else:
                logger.error("The link {0} already exists but points at {1} while you asked it to point at {2}. Aborting...".format(link, os.readlink(link), target_dir))
                quit()
        else:
            logger.error("{0} is an existing directory and cannot be turned into a link. Aborting...".format(link))
            quit()
    else:
        #needs to point at vis on scratch for gpubox files
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
                 product_dir, parallel, args,
                 ics=False, n_untar=2, keep="", vcstools_version="master",
                 nice=0):

    vcs_database_id = database_vcs.database_command(args, obsid)
    #Load computer dependant config file
    comp_config = config.load_config_file()

    logger.info("Downloading files from archive")
    # voltdownload = distutils.spawn.find_executable("voltdownload.py") #Doesn't seem to be working on zeus for some reason
    voltdownload = "voltdownload.py"
    obsinfo = meta.getmeta(service='obs', params={'obs_id':str(obsid)})
    data_format = obsinfo['dataquality']
    if data_format == 1:
        target_dir = link = '/raw'
        if ics:
            logger.error("Data have not been recombined in the "+\
                         "archive yet. Exiting")
            quit()
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
        quit()
    mdir(dl_dir, dir_description)
    create_link(data_dir, target_dir, product_dir, link)
    batch_dir = product_dir+"/batch/"

    for time_to_get in range(start_time,stop_time,increment):
        if time_to_get + increment > stop_time:
            increment = stop_time - time_to_get + 1
        get_data = "{0} --obs={1} --type={2} --from={3} --duration={4} --parallel={5} --dir={6}".\
                   format(voltdownload,obsid, data_type, time_to_get,
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
            volt_secs_to_run = datetime.timedelta(seconds=300*increment)
            check_secs_to_run = "15:00"
            if data_type == 16:
                check_secs_to_run = "10:15:00"
                #tar_batch = "untar_{0}".format(time_to_get)
                #tar_secs_to_run = "10:00:00"
                #body = []
                #untar = distutils.spawn.find_executable('untar.sh')
                #body.append("export CMD_VCS_DB_FILE={0}".format(os.environ['CMD_VCS_DB_FILE']])
                #body.append(database_vcs.add_database_function())
                #body.append('run "{0}"  "-w {1} -o {2} -b {3} -e {4} -j {5} {6}" "{7}"'.format(
                #            untar, dl_dir, obsid, time_to_get, 
                #            time_to_get+increment-1, n_untar, keep, vcs_database_id))

                #submit_slurm(tar_batch, body, batch_dir=batch_dir,
                #             slurm_kwargs={"time": str(tar_secs_to_run), "partition": "workq",
                #                           "nice": nice},
                #             vcstools_version=vcstools_version, submit=False,
                #             outfile=batch_dir+tar_batch+".out", cluster="zeus", export="NONE")

            checks = distutils.spawn.find_executable("checks.py")
            # Write out the checks batch file but don't submit it
            commands = []
            #commands.append("module load numpy")
            commands.append("export CMD_VCS_DB_FILE={0}".format(os.environ['CMD_VCS_DB_FILE']))
            commands.append(database_vcs.add_database_function())
            commands.append("newcount=0")
            commands.append("let oldcount=$newcount-1")
            commands.append("sed -i -e \"s/oldcount=${{oldcount}}/oldcount=${{newcount}}/\" {0}".\
                            format(batch_dir+voltdownload_batch+".batch"))
            commands.append("oldcount=$newcount; let newcount=$newcount+1")
            commands.append("sed -i -e \"s/_${{oldcount}}.out/_${{newcount}}.out/\" {0}".\
                            format(batch_dir+voltdownload_batch+".batch"))
            commands.append('run "{0}" "-m download -o {1} -w {2} -b {3} -i {4} --data_type {5}" "{6}"'.format(checks, obsid, dl_dir, time_to_get, increment, str(data_type),vcs_database_id))
            commands.append("if [ $? -eq 1 ];then")
            commands.append("sbatch {0}".format(batch_dir+voltdownload_batch+".batch"))
            # if we have tarballs we send the untar jobs to the workq
            if data_type == 16:
                commands.append("else")
                untar = distutils.spawn.find_executable('untar.sh')
                commands.append('run "{0}" "-w {1} -o {2} '
                                '-b {3} -e {4} -j {5} {6}" "{7}"'.format(untar, dl_dir,
                                       obsid, time_to_get, time_to_get+increment-1, n_untar,
                                       keep, vcs_database_id))

                #commands.append("sbatch {0}.batch".format(batch_dir+tar_batch))
            commands.append("fi")

            # Download and checks should be done on Zeus's cpuq. This will only work
            # on Galaxy as the Ozstar workflow is different
            submit_slurm(check_batch, commands, batch_dir=batch_dir, 
                         slurm_kwargs={"time": check_secs_to_run, 
                                       "nice": nice, 
                                       "mem-per-cpu": "10GB"},
                         vcstools_version=vcstools_version, submit=False,
                         outfile=batch_dir+check_batch+"_0.out", 
                         queue="zcpuq", export="NONE")

            # Write out the tar batch file if in mode 15
            #if format == 16:
            #        body = []
            #        for t in range(time_to_get, time_to_get+increment):
            #                body.append("aprun tar -xf {0}/1149620392_{1}_combined.tar".format(dl_dir,t))
            #        submit_slurm(tar_batch,body,batch_dir=working_dir+"/batch/", slurm_kwargs={"time":"1:00:00", "partition":"gpuq" })



            module_list=["mwa-voltage/master"]
            body = []
            body.append("export CMD_VCS_DB_FILE={0}".format(os.environ['CMD_VCS_DB_FILE']))
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
            body.append('run "{0}" "--obs={1} --type={2} --from={3} --duration={4} --parallel={5} --dir={6}" "{7}"'.\
                        format(voltdownload,obsid, data_type, time_to_get,
                               (increment-1),parallel, dl_dir, vcs_database_id))
            submit_slurm(voltdownload_batch, body, batch_dir=batch_dir, 
                         module_list=module_list, 
                         slurm_kwargs={"time": str(volt_secs_to_run), 
                                       "nice" : nice},
                         vcstools_version=vcstools_version, 
                         outfile=batch_dir+voltdownload_batch+"_1.out",
                         queue="copyq", export="NONE")

            # submit_cmd = subprocess.Popen(volt_submit_line,shell=True,stdout=subprocess.PIPE)
            continue
    # TODO: what is the below doing here???
        try:
            os.chdir(product_dir)
        except:
            logging.error("cannot open working dir:{0}".format(product_dir))
            sys.exit()



def download_cal(obs_id, cal_obs_id, data_dir, product_dir, args, head=False,
                 vcstools_version="master", nice=0):

    vcs_database_id = database_vcs.database_command(args, obs_id)
    #Load computer dependant config file
    comp_config = config.load_config_file()

    batch_dir = product_dir + '/batch/'
    product_dir = '{0}/cal/{1}'.format(product_dir,cal_obs_id)
    mdir(product_dir, 'Calibrator product')
    mdir(batch_dir, 'Batch')
    # obsdownload creates the folder cal_obs_id regardless where it runs
    # this deviates from our initially inteded naming conventions of
    # /astro/mwaopos/vcs/[cal_obs_id]/vis but the renaming and linking is a pain otherwise,
    # hence we'll link vis agains /astro/mwaopos/vcs/[cal_obs_id]/[cal_obs_id]
    target_dir = '{0}'.format(cal_obs_id)
    link = 'vis'
    csvfile = "{0}{1}_dl.csv".format(batch_dir,cal_obs_id)
    obsdownload = distutils.spawn.find_executable("obsdownload.py")
    get_data = "{0} -o {1} -d {2}".format(obsdownload,cal_obs_id, data_dir)
    if head:
        logging.error("I'm sorry, this option is no longer supported. Please "+\
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
        make_link = "ln -s {0} {1}/{2}".format(data_dir, product_dir, link)
        obsdownload_batch = "caldownload_{0}".format(cal_obs_id)
        secs_to_run = "03:00:00" # sometimes the staging can take a while...
        module_list = ["setuptools"]
        commands = []
        commands.append("module load manta-ray-client")
        commands.append("export CMD_VCS_DB_FILE={0}".format(os.environ['CMD_VCS_DB_FILE']))
        commands.append(database_vcs.add_database_function())
        commands.append("csvfile={0}".format(csvfile))
        # commands.append('source /group/mwaops/PULSAR/psrBash.profile')
        # commands.append('module load setuptools')
        commands.append('cd {0}'.format(data_dir))
        commands.append('if [[ -z ${MWA_ASVO_API_KEY} ]]')
        commands.append('then')
        commands.append('    echo "Error, MWA_ASVO_API_KEY not set"')
        commands.append('    echo "Cannot use client"')
        commands.append('    echo "Please read the MWA ASVO documentation '+\
                        'about setting this (https://wiki.mwatelescope.org/'+\
                        'display/MP/MWA+ASVO%3A+Release+Notes)"')
        commands.append('    exit 1')
        commands.append('fi')
        commands.append('echo "obs_id={0}, job_type=d, download_type=vis" > {1}'.\
                        format(cal_obs_id,csvfile))
        commands.append('mwa_client --csv={0} --dir={1}'.format(csvfile,data_dir))
        # commands.append('run "{0}" "-o {1} -d {2}" "{3}"'.format(obsdownload,cal_obs_id, data_dir,vcs_database_id))
        commands.append(make_link)
        commands.append('unzip *.zip')
        submit_slurm(obsdownload_batch, commands, batch_dir=batch_dir, 
                     module_list=module_list,
                     slurm_kwargs={"time": secs_to_run, "nice": nice},
                     vcstools_version=vcstools_version, queue="copyq", 
                     export="NONE", mem=2048)


def vcs_recombine(obsid, start_time, stop_time, increment, data_dir, product_dir, args,
                  vcstools_version="master", nice=0):

    vcs_database_id = database_vcs.database_command(args, obsid)
    #Load computer dependant config file
    comp_config = config.load_config_file()
    
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
        commands.append("export CMD_VCS_DB_FILE={0}".format(os.environ['CMD_VCS_DB_FILE']))
        commands.append(database_vcs.add_database_function())
        commands.append("newcount=0")
        commands.append("let oldcount=$newcount-1")
        commands.append("sed -i -e \"s/oldcount=${{oldcount}}/oldcount=${{newcount}}/\" {0}".\
                        format(batch_dir+recombine_batch+".batch"))
        commands.append("oldcount=$newcount; let newcount=$newcount+1")
        commands.append("sed -i -e \"s/_${{oldcount}}.out/_${{newcount}}.out/\" {0}".\
                        format(batch_dir+recombine_batch+".batch"))
        commands.append('run "{0}" "-m recombine -o {1} -w {2}/combined/ -b {3} -i {4}" "{5}"'.\
                        format(checks, obsid, data_dir, time_to_get, 
                               process_nsecs, vcs_database_id)) 
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
        commands.append("export CMD_VCS_DB_FILE={0}".format(os.environ['CMD_VCS_DB_FILE']))
        commands.append(database_vcs.add_database_function())
        #commands.append("module switch PrgEnv-cray PrgEnv-gnu")
        #commands.append("module load mpi4py")
        #commands.append("module load cfitsio")
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
        commands.append('run "srun --export=all python3 {0}" "-o {1} -s {2} -w {3} -e {4}" "{5}"'.\
                        format(recombine, obsid, time_to_get, data_dir,
                               recombine_binary, vcs_database_id))

        submit_slurm(recombine_batch, commands, batch_dir=batch_dir, 
                     module_list=module_list,
                     slurm_kwargs={"time": "06:00:00", "nodes": str(nodes), 
                                   "ntasks-per-node": jobs_per_node, 
                                   "nice": nice},
                     vcstools_version=vcstools_version, 
                     outfile=batch_dir+recombine_batch+"_1.out",
                     queue='gpuq', export="NONE")




def vcs_correlate(obsid,start,stop,increment, data_dir, product_dir, ft_res, args, metafits,
                  vcstools_version="master", nice=0):

    vcs_database_id = database_vcs.database_command(args, obsid)
    logger.info("Correlating files at {0} kHz and {1} milliseconds".\
                format(ft_res[0], ft_res[1]))
    from astropy.time import Time
    import calendar

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
                body.append("export CMD_VCS_DB_FILE={0}".format(os.environ['CMD_VCS_DB_FILE']))
                body.append(database_vcs.add_database_function())
                #body.append("source /group/mwaops/PULSAR/psrBash.profile")
                #body.append("module swap craype-ivybridge craype-sandybridge")

                # with open(corr_batch, 'w') as batch_file:
                #     batch_file.write("#!/bin/bash -l\n#SBATCH --nodes=1\n#SBATCH --account=mwaops\n#SBATCH --export=NONE\n#SBATCH --output={0}.out\n".format(corr_batch[:-6]))
                #     batch_file.write('source /group/mwaops/PULSAR/psrBash.profile\n')
                #     batch_file.write('module swap craype-ivybridge craype-sandybridge\n')
                to_corr = 0
                for file in f:
                    corr_line = ""
                    (current_time,ext) = os.path.splitext(os.path.basename(file))
                    (obsid,gpstime,chan) = current_time.split('_')
                    t = Time(int(gpstime), format='gps', scale='utc')
                    unix_time = int(t.unix)

                    body.append('run "srun --export=all {0}" "-o {1}/{2} -s {3} -r {4} -i {5} -f 128 -n {6} -c {7:0>2} -d {8}" "{9}"'.\
                                format("offline_correlator", corr_dir, obsid, 
                                       unix_time, num_frames, integrations, 
                                       int(ft_res[0]/10), gpubox_label, file,
                                       vcs_database_id))
                    to_corr += 1
                    # with open(corr_batch, 'a') as batch_file:
                    #     batch_file.write(corr_line)
                    #     to_corr = to_corr+1

                module_list = ["module switch PrgEnv-cray PrgEnv-gnu"]
                secs_to_run = str(datetime.timedelta(seconds=2*12*num_frames*to_corr)) 
                # added factor two on 10 April 2017 as galaxy seemed really slow...
                submit_slurm(corr_batch, body, module_list=module_list,
                             slurm_kwargs={"time": secs_to_run, "nice": nice,
                                           "gres": "gpu:1"},
                             queue='gpuq', vcstools_version=vcstools_version, 
                             batch_dir=batch_dir, export="NONE")
                # batch_submit_line = "sbatch --workdir={0} --time={1} --partition=gpuq --gid=mwaops {2} \n".format(corr_dir,secs_to_run,corr_batch)
                # submit_cmd = subprocess.Popen(batch_submit_line,shell=True,stdout=subprocess.PIPE)
                # jobid=""
                # for line in submit_cmd.stdout:
                #     if "Submitted" in line:
                #         (word1,word2,word3,jobid) = line.split()
            else:
                logger.error("Couldn't find any recombine files. Aborting here.")


def coherent_beam(obs_id, start, stop, data_dir, product_dir, batch_dir, 
                  metafits_file, nfine_chan, pointing, args, 
                  rts_flag_file=None, bf_formats=None, DI_dir=None, 
                  execpath=None, calibration_type='rts', 
                  vcstools_version="master", nice=0):
    """
    This function runs the new version of the beamformer. It is modelled after 
    the old function above and will likely be able to be streamlined after 
    working implementation (SET)

    Streamlining underway, as well as full replacement of the old function (SET March 28, 2018)
    """
    vcs_database_id = database_vcs.database_command(args, obs_id)  
    
    #Load computer dependant config file
    comp_config = config.load_config_file()

    # If execpath is given, change the make_beam executable command
    # otherwise, it should be on the PATH if vcstools has been installed
    if execpath:
        make_beam_cmd = "{0}/make_beam".format(execpath)
        make_beam_version_cmd = "{0}/make_beam -V".format(execpath)
    else:
        make_beam_cmd = "make_beam"
        make_beam_version_cmd = "make_beam -V"

    #make_beam_version_cmd = "make_beam -V"
    make_beam_version = subprocess.Popen(make_beam_version_cmd, 
                           stdout=subprocess.PIPE, shell=True).communicate()[0]
    tested_version = "?.?.?"
    logger.info("Current version of make_beam = {0}".format(make_beam_version.strip()))
    logger.info("Tested version of make_beam = {0}".format(tested_version.strip()))

    metafile = "{0}/{1}.meta".format(product_dir, obs_id)
    metafile_exists = False
    if os.path.isfile(metafile):
        logger.info("Found observation metafile: {0}".format(metafile))
        channels = None
        with open(metafile, 'r') as m:
            for chan_line in m.readlines():
                if chan_line.startswith("channels"):
                    channels = chan_line.split(",")[1:]
        if channels == None:
            logger.info("Channels keyword not found in metafile. Re-querying "+\
                        "the database.")
        else:
            metafile_exists = True
    # Grabbing this from the calibrate section for now. This should be streamlined to call an external function (SET)
    if metafile_exists == False:
        logger.info("Querying the database for calibrator obs ID {0}...".\
                    format(obs_id))
        obs_info = meta.getmeta(service='obs', params={'obs_id': str(obs_id)})
        channels = obs_info[u'rfstreams'][u"0"][u'frequencies']
        with open(metafile, "w") as m:
            m.write("#Metadata for obs ID {0} required to determine if: normal or picket-fence\n".format(obs_id))
            m.write("channels,{0}".format(",".join([str(c) for c in channels])))
    channels = np.array(channels, dtype=np.int)
    hichans = [c for c in channels if c>128]
    lochans = [c for c in channels if c<=128]
    lochans.extend(list(reversed(hichans)))
    ordered_channels = lochans

    # Run for each coarse channel. Calculates delays and makes beam

    if not DI_dir:
        logger.error("You need to specify the path to the calibrator files, "+\
                     "either where the DIJs are or where the Offringa "+\
                     "calibration_solutions.bin file is. Aborting here")
        quit()
    DI_dir = os.path.abspath(DI_dir)
    RA = pointing[0]
    Dec = pointing[1]

    # make_beam_small requires the start time in UTC, get it from the start
    # GPS time as is done in timeconvert.py
    utctime = Time(start, format='gps', scale='utc').fits
    # remove (UTC) that some astropy versions leave on the end
    if utctime.endswith('(UTC)'):
        utctime = strptime(utctime, '%Y-%m-%dT%H:%M:%S.000(UTC)')
        utctime = strftime('%Y-%m-%dT%H:%M:%S', utctime)
    else:
        utctime = strptime(utctime, '%Y-%m-%dT%H:%M:%S.000')
        utctime = strftime('%Y-%m-%dT%H:%M:%S', utctime)


    logging.info("Running make_beam")
    P_dir = product_dir+"/pointings"
    mdir(P_dir, "Pointings")
    pointing_dir = "{0}/{1}_{2}".format(P_dir, RA, Dec)
    mdir(pointing_dir, "Pointing {0} {1}".format(RA, Dec))
    # startjobs = True


    # set up SLURM requirements
    seconds_to_run = 15 * (stop - start + 1)  # This should be able to be reduced after functional testing
    if seconds_to_run > 86399.:
        secs_to_run = datetime.timedelta(seconds=86399)
    else:
        secs_to_run = datetime.timedelta(seconds=seconds_to_run)

    # VDIF will need gpuq if inverting pfb with '-m' option, otherwise cpuq is fine
    # In general this needs to be cleaned up, prefferably to be able to intelligently select a
    # queue and a maximum wall time (SET)
    n_omp_threads = 1
    openmp_line = "export OMP_NUM_THREADS={0}".format(n_omp_threads)

    # Run one coarse channel per node
    #for coarse_chan in range(24):
    job_id_list = []
    for gpubox, coarse_chan in enumerate(ordered_channels, 1):
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
            quit()

        make_beam_small_batch = "mb_{0}_{1}_ch{2}".format(RA, Dec, coarse_chan)
        module_list = [comp_config['container_module']]
        commands = []
        #commands.append("source /group/mwaops/PULSAR/psrBash.profile")
        #commands.append("module swap craype-ivybridge craype-sandybridge")
        commands.append(openmp_line)
        commands.append("cd {0}".format(pointing_dir))
        commands.append("srun --export=all -n 1 -c {0} {1} {2} -o {3} -b {4} -e {5} -a 128 -n 128 -f {6} {7} -d "
                        "{8}/combined -R {9} -D {10} -r 10000 -m {11} -z {12} {13} -F {14}".format(n_omp_threads, comp_config['container_command'], make_beam_cmd, obs_id, start,
                        stop, coarse_chan, jones_option, data_dir, RA, Dec, metafits_file, utctime, bf_formats, rts_flag_file))  # check these

        job_id = submit_slurm(make_beam_small_batch, commands,
                    batch_dir=batch_dir, module_list=module_list,
                    slurm_kwargs={"time":secs_to_run, "nice":nice},
                    queue='gpuq', vcstools_version=vcstools_version, 
                    submit=True, export="NONE", gpu_res=1)
        job_id_list.append(job_id)
    #TODO This can be returned as a job id string that can be slapped right on for dependancies
    #job_id_str = ""
    #for i in job_id_list:
    #    job_id_str += ":" + str(i)
    #return job_id_str
    return job_id_list

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

    from optparse import OptionParser, OptionGroup, SUPPRESS_HELP

 #   parser=OptionParser(description="process_vcs.py is a script of scripts that downloads prepares and submits jobs to Galaxy. It can be run with just a pointing (-p \"xx:xx:xx xx:xx:xx.x\") and an obsid (\"-o <obsid>\") and it will process all the data in the obs. It will call prepare.py which will attempt to build the phase and calibration information - which will only exist if a calibration obs has already been run. So it will only move past the data prepa stage if the \"-r\" (for run) is used\n"

    parser=OptionParser(description="process_vcs.py is a script for processing the MWA VCS data on Galaxy in steps. It can download data from the archive, call on recombine to form course channels, run the offline correlator, make tile re-ordered and bit promoted PFB files or for a coherent beam for a given pointing.")
    group_download = OptionGroup(parser, 'Download Options')
    group_download.add_option("--head", action="store_true", default=False, help="Submit download jobs to the headnode instead of the copyqueue [default=%default]")
    #group_download.add_option("--format", type="choice", choices=['11','15','16'], default='11', help="Voltage data type (Raw = 11, ICS Only = 15, Recombined and ICS = 16) [default=%default]")
    group_download.add_option("-d", "--parallel_dl", type="int", default=3, help="Number of parallel downloads to envoke [default=%default]")
    group_download.add_option("-j", "--untar_jobs", type='int', default=2, help="Number of parallel jobs when untaring downloaded tarballs. [default=%default]")
    group_download.add_option("-k", "--keep_tarball", action="store_true", default=False, help="Keep the tarballs after unpacking. [default=%default]")
    group_correlate = OptionGroup(parser, 'Correlator Options')
    group_correlate.add_option("--ft_res", metavar="FREQ RES,TIME RES", type="int", nargs=2, default=(10,1000), help="Frequency (kHz) and Time (ms) resolution for running the correlator. Please make divisible by 10 kHz and 10 ms respectively. [default=%default]")

    group_beamform = OptionGroup(parser, 'Beamforming Options')
    group_beamform.add_option("-p", "--pointing", nargs=2, help="required, R.A. and Dec. of pointing, e.g. \"19:23:48.53\" \"-20:31:52.95\"")
    group_beamform.add_option("--DI_dir", default=None, help="Directory containing either Direction Independent Jones Matrices (as created by the RTS) " +\
                                  "or calibration_solution.bin as created by Andre Offringa's tools.[no default]")
    group_beamform.add_option("--bf_out_format", type="choice", choices=['psrfits','vdif','both'], help="Beam former output format. Choices are {0}. [default=%default]".format(bf_out_modes), default='psrfits')
    group_beamform.add_option("--incoh", action="store_true", default=False, help="Add this flag if you want to form an incoherent sum as well. [default=%default]")
    group_beamform.add_option("--flagged_tiles", type="string", default=None, help="Path (including file name) to file containing the flagged tiles as used in the RTS, will be used by get_delays. [default=%default]")
    group_beamform.add_option('--cal_type', type='string', help="Use either RTS (\"rts\") solutions or Andre-Offringa-style (\"offringa\") solutions. Default is \"rts\". If using Offringa's tools, the filename of calibration solution must be \"calibration_solution.bin\".", default="rts")
    group_beamform.add_option("-E", "--execpath", type="string", default=None, help="Supply a path into this option if you explicitly want to run files from a different location for testing. Default is None (i.e. whatever is on your PATH).")

    parser.add_option("-m", "--mode", type="choice", choices=['download','download_ics', 'download_cal', 'recombine','correlate', 'calibrate', 'beamform'], help="Mode you want to run. {0}".format(modes))
    parser.add_option("-o", "--obs", metavar="OBS ID", type="int", help="Observation ID you want to process [no default]")
    parser.add_option('--cal_obs', '-O', metavar="CALIBRATOR OBS ID", type="int", help="Only required in 'calibrate' and 'download_cal' mode."+\
                          "Observation ID of calibrator you want to process. In case of " + \
                          "in-beam calibration should be the same as input to -o (obsID). [no default]", default=None)
    parser.add_option("-b", "--begin", type="int", help="First GPS time to process [no default]")
    parser.add_option("-e", "--end", type="int", help="Last GPS time to process [no default]")
    parser.add_option("-a", "--all", action="store_true", default=False, help="Perform on entire observation span. Use instead of -b & -e. [default=%default]")
    parser.add_option("-i", "--increment", type="int", default=64, help="Increment in seconds (how much we process at once) [default=%default]")
    parser.add_option("-s", action="store_true", default=False, help="Single step (only process one increment and this is it (False == do them all) [default=%default]")
    parser.add_option("-w", "--work_dir", metavar="DIR", default=None, help="Base directory you want run things in. USE WITH CAUTION! Per default " + \
                          "raw data will will be downloaded into /astro/mwaops/vcs/[obsID] and data products will be in /group/mwaops/vcs/[obsID]."+ \
                          " If set, this will create a folder for the Obs. ID if it doesn't exist [default=%default]")
    parser.add_option("-c", "--ncoarse_chan", type="int", default=24, help="Coarse channel count (how many to process) [default=%default]")
    parser.add_option("-n", "--nfine_chan", type="int", default=128, help="Number of fine channels per coarse channel [default=%default]")
    parser.add_option("--mail",action="store_true", default=False, help="Enables e-mail notification about start, end, and fail of jobs. Currently only implemented for beamformer mode.[default=%default]")
    parser.add_option("-L", "--loglvl", type="string", help="Logger verbosity level. Default: INFO", 
                                        default="INFO")
    parser.add_option("-V", "--version", action="store_true", help="Print version and quit")
    parser.add_option("--vcstools_version", type="string", default="master", help="VCSTools version to load in jobs (i.e. on the queues) [default=%default]")
    parser.add_option("--nice", type="int", default=0, help="Reduces your priority of Slurm Jobs. [default=%default]")
    parser.add_option_group(group_download)
    parser.add_option_group(group_correlate)
    parser.add_option_group(group_beamform)

    (opts, args) = parser.parse_args()

    # set up the logger for stand-alone execution
    logger.setLevel(loglevels[opts.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[opts.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.propagate = False

    logger.info("Using vcstools/{0}".format(opts.vcstools_version))
    if opts.version:
        try:
            import version
            logger.info(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            logger.error("Couldn't import version.py - have you installed vcstools?")
            logger.error("ImportError: {0}".format(ie))
            sys.exit(0)

    #Option parsing
    if opts.all and (opts.begin or opts.end):
        logger.error("Please specify EITHER (-b,-e) OR -a")
        quit()
    elif opts.all:
        opts.begin, opts.end = meta.obs_max_min(opts.cal_obs\
                               if opts.mode == 'download_cal' else opts.obs)
    # make sure we can process increments smaller than 64 seconds when not in calibration related mode
    if opts.mode not in ['download_cal','calibrate']:
        if opts.end - opts.begin +1 < opts.increment:
            opts.increment = opts.end - opts.begin + 1
    e_mail = ""
    if opts.mail:
        e_mail = get_user_email()
        logger.info("Sending info to {0}".format(e_mail))
    if not opts.mode:
      logger.error("Mode required {0}. Please specify with -m or --mode.".format(modes))

      quit()
    if not opts.obs:
        logger.error("Observation ID required, please put in with -o or --obs")
        quit()
    if opts.begin and opts.end:
        if opts.begin > opts.end:
            logger.error("Starting time is after end time")
            quit()
    if (opts.mode == "beamform" or opts.incoh):
        bf_format = ""
        if not opts.pointing:
            logger.error("Pointing (-p) required in beamformer mode")
            quit()
        if (opts.bf_out_format == 'psrfits' or opts.bf_out_format == 'both'):
            bf_format +=" -p"
            logger.info("Writing out PSRFITS.")
        if  (opts.bf_out_format == 'vdif' or opts.bf_out_format == 'both'):
            bf_format += " -u"
            logger.info("Writing out upsampled VDIF.")
        if (opts.incoh):
            bf_format += " -i"
            logger.info("Writing out incoherent sum.")

        # This isn't necessary as checks for execpath are done in beamforming function (BWM 6/4/18)
        #if opts.execpath:
        #    execpath = opts.execpath

    #Load computer dependant config file
    comp_config = config.load_config_file()
    
    if opts.work_dir:
        logger.warning("YOU ARE MESSING WITH THE DEFAULT DIRECTORY STRUCTURE "+\
                       "FOR PROCESSING -- BE SURE YOU KNOW WHAT YOU ARE DOING!")
        sleep(5)
        data_dir = product_dir = "{0}/{1}".format(opts.work_dir, opts.obs)
    else:
        data_dir = '{0}{1}'.format(comp_config['base_data_dir'],opts.obs)
        product_dir = '{0}{1}'.format(comp_config['base_product_dir'],opts.obs)
    batch_dir = "{0}/batch".format(product_dir)
    mdir(data_dir, "Data")
    mdir(product_dir, "Products")
    mdir(batch_dir, "Batch")
    metafits_file = "{0}/{1}_metafits_ppds.fits".format(data_dir, opts.obs)
    # TODO: modify metafits downloader to not just do a trivial wget

    logger.info("Processing Obs ID {0} from GPS times {1} till {2}".\
                format(opts.obs, opts.begin, opts.end))

    if opts.mode == 'download_ics':
        logger.info("Mode: {0}".format(opts.mode))
        vcs_download(opts.obs, opts.begin, opts.end, opts.increment, opts.head,
                     data_dir, product_dir, opts.parallel_dl, sys.argv,
                     ics=True, vcstools_version=opts.vcstools_version,
                     nice=opts.nice)
    elif opts.mode == 'download':
        logger.info("Mode: {0}".format(opts.mode))
        vcs_download(opts.obs, opts.begin, opts.end, opts.increment, opts.head,
                     data_dir, product_dir, opts.parallel_dl, sys.argv,
                     n_untar=opts.untar_jobs,
                     keep='-k' if opts.keep_tarball else "",
                     vcstools_version=opts.vcstools_version, nice=opts.nice)
    elif opts.mode == 'recombine':
        logger.info("Mode: {0}".format(opts.mode))
        ensure_metafits(data_dir, opts.obs, metafits_file)
        vcs_recombine(opts.obs, opts.begin, opts.end, opts.increment, data_dir,
                      product_dir, sys.argv, 
                      vcstools_version=opts.vcstools_version, nice=opts.nice)
    elif opts.mode == 'correlate':
        logger.info("Mode: {0}".format(opts.mode))
        ensure_metafits(data_dir, opts.obs, metafits_file)
        vcs_correlate(opts.obs, opts.begin, opts.end, opts.increment, data_dir,
                      product_dir, opts.ft_res, sys.argv, metafits_file, 
                      vcstools_version=opts.vcstools_version, nice=opts.nice)
    elif opts.mode == 'download_cal':
        logger.info("Mode: {0}".format(opts.mode))
        if not opts.cal_obs:
            logger.error("You need to also pass the calibrator observation ID."+\
                         " Aborting here.")
            quit()
        if opts.cal_obs == opts.obs:
            logging.error("The calibrator obsID cannot be the same as the "+\
                          "target obsID -- there are not gpubox files for "+\
                          "VCS data on the archive.")
            quit()
        data_dir = data_dir.replace(str(opts.obs), str(opts.cal_obs))
        mdir(data_dir, "Calibrator Data")
        download_cal(opts.obs, opts.cal_obs, data_dir, product_dir, sys.argv, 
                     opts.head, vcstools_version=opts.vcstools_version, 
                     nice=opts.nice)
    elif opts.mode == ('beamform' or 'incoh'):
        logger.info("Mode: {0}".format(opts.mode))
        if not opts.DI_dir:
            logger.error("You need to specify the path to either where the "+\
                         "DIJs are or where the offringe calibration_solution.bin"+\
                         "file is. Aborting here.")
            quit()
        if opts.flagged_tiles:
            flagged_tiles_file = os.path.abspath(opts.flagged_tiles)
            if not os.path.isfile(opts.flagged_tiles):
                logger.error("Your are not pointing at a file with your input "+\
                             "to --flagged_tiles. Aborting here as the "+\
                             "beamformer will not run...")
                quit()
        else:
            flagged_tiles_file = None
        ensure_metafits(data_dir, opts.obs, metafits_file)
        coherent_beam(opts.obs, opts.begin, opts.end, data_dir, product_dir, 
                      batch_dir, metafits_file, opts.nfine_chan, opts.pointing,
                      sys.argv, rts_flag_file=flagged_tiles_file, 
                      bf_formats=bf_format, DI_dir=opts.DI_dir, 
                      calibration_type=opts.cal_type,
                      vcstools_version=opts.vcstools_version, nice=opts.nice,
                      execpath=opts.execpath)
    else:
        logger.error("Somehow your non-standard mode snuck through. "+\
                     "Try again with one of {0}".format(modes))
        quit()
