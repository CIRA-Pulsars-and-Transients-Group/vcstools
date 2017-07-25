#!/usr/bin/env python
import subprocess
import os
import sys
import glob
import time
import tempfile
import atexit
import hashlib
import datetime
import time
import distutils.spawn
import sqlite3 as lite
from astropy.io import fits as pyfits
from reorder_chans import *


TMPL = """#!/bin/bash

#SBATCH --export=NONE
#SBATCH --output={outfile}
#SBATCH --account=mwaops
#SBATCH --gid=mwaops
#SBATCH --clusters={cluster}
#
{header}

{script}
"""


def tmp(suffix=".sh"):
    t = tempfile.mktemp(suffix=suffix)
    atexit.register(os.unlink, t)
    return t

def submit_slurm(name,commands,slurm_kwargs={},tmpl=TMPL,batch_dir="batch/", depend=0, submit=True, outfile=None, cluster="galaxy"):
        """
        Making this function to cleanly submit slurm jobs using a simple template.
        This will use the <name> to setup both the .batch and .out files into <batch_dir>
        <slurm_kwargs> should be a dictionary of keyword, value pairs of anything the template header is missing
        <commands> is the actual batch script you want to run, and is expecting a list with one entry per line
        """
        
        header = []
        if batch_dir[-1] is not "/":
                batch_dir += "/"

        if not outfile:
                outfile=batch_dir+name+".out"
        for k, v in slurm_kwargs.iteritems():
                if len(k) > 1:
                        k = "--" + k + "="
                else:
                        k = "-" + k + " "
                header.append("#SBATCH {0}{1}".format(k, v))
        header = "\n".join(header)
        
        commands = "\n".join(commands)
        
        tmpl = tmpl.format(script=commands,outfile=outfile, header=header, cluster=cluster)
        

        fh = open(batch_dir+name+".batch","w")
        fh.write(tmpl)
        fh.close()
        
        if depend:
                batch_submit_line = "sbatch --dependency=afterok:{0} {1}".format(depend,batch_dir+name+".batch") # should this just be in the header?
        else:
                batch_submit_line = "sbatch {0}".format(batch_dir+name+".batch")
        
        if submit:
                submit_cmd = subprocess.Popen(batch_submit_line,shell=True,stdout=subprocess.PIPE)
        # submit_cmd = subprocess.Popen(batch_submit_line,shell=True,stdout=subprocess.PIPE)
        # 
        # jobid=""
        # for line in submit_cmd.stdout:
        #         if "Submitted" in line:
        #                 (word1,word2,word3,jobid) = line.split()
        # return jobid

        


def getmeta(service='obs', params=None):
    """
    Function to call a JSON web service and return a dictionary:
    Given a JSON web service ('obs', find, or 'con') and a set of parameters as
    a Python dictionary, return a Python dictionary xcontaining the result.
    Taken verbatim from http://mwa-lfd.haystack.mit.edu/twiki/bin/view/Main/MetaDataWeb
    """
    import urllib
    import urllib2
    import json

    # Append the service name to this base URL, eg 'con', 'obs', etc.
    BASEURL = 'http://mwa-metadata01.pawsey.org.au/metadata/'


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

def mdir(path,description, gid=30832):
    # the default groupID is mwaops which is 30832 in numerical
    # we try and make sure all directories created by process_vcs
    # end up belonging to the user and the group mwaops
    # with rwx permissions and the sticky bit set for both user and group
    try:
        os.makedirs(path)
        # we leave the uid unchanged but change gid to mwaops
        os.chown(path,-1,gid)
        os.chmod(path,0771)
        os.system("chmod -R g+s {0}".format(path))
    except:
        if (os.path.exists(path)):
            print "{0} Directory Already Exists\n".format(description)
        else:
            sys.exit()

def get_user_email():
    command="echo `ldapsearch -x \"uid=$USER\" mail |grep \"^mail\"|cut -f2 -d' '`"
    email = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True).communicate()[0]
    return email.strip()

def ensure_metafits(data_dir, obs_id, metafits_file):
    # TODO: To get the actual ppds file should do this with obsdownload -o <obsID> -m
    if not os.path.exists(metafits_file):
        print "{0} does not exists".format(metafits_file)
        print "Will download it from the archive. This can take a while so please do not ctrl-C."
        obsdownload = distutils.spawn.find_executable("obsdownload.py")
        get_metafits = "{0} -o {1} -d {2} -m".format(obsdownload, obs_id, data_dir.replace(str(obs_id), ''))
        try:
            subprocess.call(get_metafits,shell=True)
        except:
            print "Couldn't download {0}. Aborting.".format(os.basename(metafits_file))
            quit()
        # clean up
        os.remove('obscrt.crt')
        os.remove('obskey.key')
    # make a copy of the file in the product_dir if that directory exists
    # if it doesn't we might have downloaded the metafits file of a calibrator (obs_id only exists on /astro)
    # in case --work_dir was specified in process_vcs call product_dir and data_dir
    # are the same and thus we will not perform the copy
    product_dir = data_dir.replace('/astro/mwaops/vcs/', '/group/mwaops/vcs/') # being pedantic
    if os.path.exists(product_dir) and not os.path.exists(metafits_file):
        print "Copying {0} to {1}".format(metafits_file, product_dir)
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
    target_dir = data_dir + '/' + target_dir
    # check if link exists and whether it already points to where we'd like it to
    if os.path.exists(link):
        if os.path.islink(link):
            if os.readlink(link) == target_dir:
                return
            else:
                print "The link {0} already exists but points at {1} while you asked it to point at {2}. Aborting...".format(link, os.readlink(link), target_dir)
                quit()
        else:
            print "{0} is an existing directory and cannot be turned into a link. Aborting...".format(link)
            quit()
    else:
        #needs to point at vis on scratch for gpubox files
        print "Trying to link {0} against {1}".format(link, target_dir)
        os.symlink(target_dir, link)
    

def obs_max_min(obs_id):
    """
    Small function to query the database and returns the times of the first and last file
    :param obs_id:
    :return:
    """
    from file_maxmin import getmeta

    obsinfo = getmeta(service='obs', params={'obs_id':str(obs_id)})
    times=[file[11:21] for file in obsinfo['files'] if is_number(file[11:21])] #Make a list of gps times excluding non-numbers from list
    obs_start = int(min(times))
    obs_end = int(max(times))
    return obs_start, obs_end

def get_frequencies(metafits,resort=False):
        # TODO: for robustness, this should force the entries to be 3-digit numbers
    hdulist    = pyfits.open(metafits)
    freq_str   = hdulist[0].header['CHANNELS']
    freq_array = [int(f) for f in freq_str.split(',')]
    if resort:
        return sfreq(freq_array)
    else:
        return freq_array

def vcs_download(obsid, start_time, stop_time, increment, head, data_dir, product_dir, parallel, args, ics=False, n_untar=2, keep=""):
        vcs_database_id = database_command(args, obsid)
        print "Downloading files from archive"
        voltdownload = distutils.spawn.find_executable("voltdownload.py")
        obsinfo = getmeta(service='obs', params={'obs_id':str(obsid)})
        data_format = obsinfo['dataquality']
        if data_format == 1:
                target_dir = link = '/raw'
                if ics:
                        print "Data have not been recombined in the archive yet. Exiting"
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
                print "Unable to determine data format from archive. Exiting"
                quit()
        mdir(dl_dir, dir_description)
        create_link(data_dir, target_dir, product_dir, link)
        batch_dir = product_dir+"/batch/"

        for time_to_get in range(start_time,stop_time,increment):
                get_data = "{0} --obs={1} --type={2} --from={3} --duration={4} --parallel={5} --dir={6}".format(voltdownload,obsid, data_type, time_to_get,(increment-1),parallel, dl_dir) #need to subtract 1 from increment since voltdownload wants how many seconds PAST the first one
                if head:
                        log_name="{0}/voltdownload_{1}.log".format(product_dir,time_to_get)
                        with open(log_name, 'w') as log:
                                subprocess.call(get_data, shell=True, stdout=log, stderr=log)
                else:
                        voltdownload_batch = "volt_{0}".format(time_to_get)
                        check_batch = "check_volt_{0}".format(time_to_get)
                        volt_secs_to_run = datetime.timedelta(seconds=300*increment)
                        check_secs_to_run = "15:00"
                        check_nsecs = increment if (time_to_get + increment <= stop_time) else (stop_time - time_to_get + 1)
                        if data_type == 16:
                            tar_batch = "untar_{0}".format(time_to_get)
                            tar_secs_to_run = "10:00:00"
                            body = []
                            untar = distutils.spawn.find_executable('untar.sh')
                            body.append("aprun -n 1 {0} -w {1} -o {2} -b {3} -e {4} -j {5} {6}".format(untar, dl_dir, obsid, time_to_get, time_to_get+check_nsecs-1, n_untar, keep))
                            submit_slurm(tar_batch,body,batch_dir=batch_dir, slurm_kwargs={"time":str(tar_secs_to_run), "partition":"workq"}, \
                                             submit=False, outfile=batch_dir+tar_batch+".out", cluster="galaxy")
                        checks = distutils.spawn.find_executable("checks.py")
                        # Write out the checks batch file but don't submit it
                        commands = []
                        commands.append("newcount=0")
                        commands.append("let oldcount=$newcount-1")
                        commands.append("sed -i -e \"s/oldcount=${{oldcount}}/oldcount=${{newcount}}/\" {0}".format(batch_dir+voltdownload_batch+".batch"))
                        commands.append("oldcount=$newcount; let newcount=$newcount+1")
                        commands.append("sed -i -e \"s/_${{oldcount}}.out/_${{newcount}}.out/\" {0}".format(batch_dir+voltdownload_batch+".batch"))
                        commands.append("{0} -m download -o {1} -w {2} -b {3} -i {4} --data_type {5}".format(checks, obsid, dl_dir, time_to_get, check_nsecs, str(data_type)))
                        commands.append("if [ $? -eq 1 ];then")
                        commands.append("sbatch {0}".format(batch_dir+voltdownload_batch+".batch"))
                        # if we have tarballs we send the untar jobs to the workq
                        if data_type == 16:
                            commands.append("else")
                            commands.append("sbatch {0}.batch".format(batch_dir+tar_batch))
                        commands.append("fi")
                        submit_slurm(check_batch,commands,batch_dir=batch_dir, slurm_kwargs={"time" : check_secs_to_run, "partition" : "copyq", "clusters":"zeus"}, submit=False, outfile=batch_dir+check_batch+"_0.out", cluster="zeus")
                        
                        # Write out the tar batch file if in mode 15
                        #if format == 16:
                        #        body = []
                        #        for t in range(time_to_get, time_to_get+increment):
                        #                body.append("aprun tar -xf {0}/1149620392_{1}_combined.tar".format(dl_dir,t))
                        #        submit_slurm(tar_batch,body,batch_dir=working_dir+"/batch/", slurm_kwargs={"time":"1:00:00", "partition":"gpuq" })
                                
                        
                        
                        body = []
                        body.append("oldcount=0")
                        body.append("let newcount=$oldcount+1")
                        body.append("if [ ${newcount} -gt 10 ]; then")
                        body.append("echo \"Tried ten times, this is silly. Aborting here.\";exit")
                        body.append("fi")
                        body.append("sed -i -e \"s/newcount=${{oldcount}}/newcount=${{newcount}}/\" {0}\n".format(batch_dir+check_batch+".batch"))
                        body.append("sed -i -e \"s/_${{oldcount}}.out/_${{newcount}}.out/\" {0}".format(batch_dir+check_batch+".batch"))
                        body.append("sbatch -d afterany:${{SLURM_JOB_ID}} {0}".format(batch_dir+check_batch+".batch"))
                        body.append(get_data)
                        submit_slurm(voltdownload_batch, body, batch_dir=batch_dir, slurm_kwargs={"time" : str(volt_secs_to_run), "partition" : "copyq", "clusters":"zeus"}, outfile=batch_dir+voltdownload_batch+"_1.out", cluster="zeus")
                        
                        # submit_cmd = subprocess.Popen(volt_submit_line,shell=True,stdout=subprocess.PIPE)
                        continue
                # TODO: what is the below doing here???
                try:
                        os.chdir(product_dir)
                except:
                        print "cannot open working dir:{0}".format(product_dir)
                        sys.exit()


def download_cal(obs_id, cal_obs_id, data_dir, product_dir, args, head=False):
    vcs_database_id = database_command(args, obs_id)
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

    obsdownload = distutils.spawn.find_executable("obsdownload.py")
    get_data = "{0} -o {1} -d {2}".format(obsdownload,cal_obs_id, data_dir)
    if head:
        print "Will download the data from the archive. This can take a while so please do not ctrl-C."
        log_name="{0}/caldownload_{1}.log".format(batch_dir,cal_obs_id)
        with open(log_name, 'w') as log:
            subprocess.call(get_data, shell=True, stdout=log, stderr=log)
        create_link(data_dir, target_dir, product_dir, link)
        #clean up
        try:
            os.remove('obscrt.crt')
            os.remove('obskey.key')
        except:
            pass
    else:
        # we create the link using bash and not our create_link function 
        # as we'd like to do this only once the data have arrived,
        # i.e. the download worked.
        make_link = "ln -s {0}/{1} {2}/{3}".format(data_dir, target_dir, product_dir, link)
        obsdownload_batch = "caldownload_{0}".format(cal_obs_id)
        secs_to_run = "02:00:00" # sometimes the staging can take a while... 
        commands = []
        commands.append('source /group/mwaops/PULSAR/psrBash.profile')
        commands.append('module load setuptools')
        commands.append('cd {0}'.format(data_dir))
        commands.append(get_data)
        commands.append(make_link)
        commands.append('rm obscrt.crt obskey.key')
        submit_slurm(obsdownload_batch,commands,batch_dir=batch_dir, slurm_kwargs={"time" : secs_to_run, "partition" : "copyq"}, cluster="zeus")


def vcs_recombine(obsid, start_time, stop_time, increment, data_dir, product_dir, args):
        vcs_database_id = database_command(args, obsid)
        print "Running recombine on files"
        jobs_per_node = 8
        target_dir = link = 'combined'
        mdir(data_dir + '/' + target_dir, 'Combined')
        create_link(data_dir, target_dir, product_dir, link)
        batch_dir = product_dir+"/batch/"
        recombine = distutils.spawn.find_executable("recombine.py")
        #recombine = "/group/mwaops/stremblay/galaxy-scripts/scripts/recombine.py"
        checks = distutils.spawn.find_executable("checks.py")
        recombine_binary = distutils.spawn.find_executable("recombine")
        for time_to_get in range(start_time,stop_time,increment):
        
                process_nsecs = increment if (time_to_get + increment <= stop_time) else (stop_time - time_to_get + 1)
                if (jobs_per_node > process_nsecs):
                        jobs_per_node = process_nsecs
                nodes = (increment+(-increment%jobs_per_node))//jobs_per_node + 1 # Integer division with ceiling result plus 1 for master node
                recombine_batch = "recombine_{0}".format(time_to_get)
                check_batch = "check_recombine_{0}".format(time_to_get)
                commands = []
                commands.append("newcount=0")
                commands.append("let oldcount=$newcount-1")
                commands.append("sed -i -e \"s/oldcount=${{oldcount}}/oldcount=${{newcount}}/\" {0}".format(batch_dir+recombine_batch+".batch"))
                commands.append("oldcount=$newcount; let newcount=$newcount+1")
                commands.append("sed -i -e \"s/_${{oldcount}}.out/_${{newcount}}.out/\" {0}".format(batch_dir+recombine_batch+".batch"))
                commands.append("{0} -m recombine -o {1} -w {2}/combined/ -b {3} -i {4}".format(checks, obsid, data_dir, time_to_get, process_nsecs))
                commands.append("if [ $? -eq 1 ];then")
                commands.append("sbatch {0}".format(batch_dir+recombine_batch+".batch"))  
                commands.append("fi")
                submit_slurm(check_batch,commands,batch_dir=batch_dir, slurm_kwargs={"time" : "15:00", "partition" : "gpuq"}, submit=False, outfile=batch_dir+check_batch+"_0.out")
                
                commands = []
                commands.append("module switch PrgEnv-cray PrgEnv-gnu")
                commands.append("module load mpi4py")
                commands.append("module load cfitsio")
                commands.append("oldcount=0")
                commands.append("let newcount=$oldcount+1")
                commands.append("if [ ${newcount} -gt 10 ]; then")
                commands.append("echo \"Tried ten times, this is silly. Aborting here.\";exit")
                commands.append("fi")
                commands.append("sed -i -e \"s/newcount=${{oldcount}}/newcount=${{newcount}}/\" {0}".format(batch_dir+check_batch+".batch"))
                commands.append("sed -i -e \"s/_${{oldcount}}.out/_${{newcount}}.out/\" {0}".format(batch_dir+check_batch+".batch"))
                commands.append("sbatch -d afterany:${{SLURM_JOB_ID}} {0}".format(batch_dir+check_batch+".batch")) #TODO: Add iterations?
                commands.append("aprun -n {0} -N {1} python {2} -o {3} -s {4} -w {5} -e {6}".format(process_nsecs,jobs_per_node,recombine,obsid,time_to_get,data_dir,recombine_binary))
                
                submit_slurm(recombine_batch,commands,batch_dir=batch_dir, slurm_kwargs={"time" : "06:00:00", "nodes" : str(nodes), "partition" : "gpuq"}, outfile=batch_dir+recombine_batch+"_1.out")



def vcs_correlate(obsid,start,stop,increment, data_dir, product_dir, ft_res, args, metafits):
        vcs_database_id = database_command(args, obsid)
        print "Correlating files at {0} kHz and {1} milliseconds".format(ft_res[0], ft_res[1])
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
        
        print "Input chan list is" , chan_list
        
        for time_to_get in range(start,stop,increment):
                inc_start = time_to_get
                inc_stop = time_to_get+increment
                for index,channel in enumerate(chan_list):
                        gpubox_label = (index+1)
                        f=[]
                        for time_to_corr in range(inc_start,inc_stop,1):
                                file_to_process = "{0}/combined/{1}_{2}_ch{3:0>2}.dat".format(data_dir,obsid,time_to_corr,channel)
                                #check the file exists
                                if (os.path.isfile(file_to_process) == True):
                                        f.append(file_to_process)
        
                        #now have a full list of files
                        #for this increment 
                        #and this channel
                        if (len(f) > 0):
                                corr_batch = "correlator_{0}_gpubox{1:0>2}".format(inc_start,gpubox_label)
                                body = []
                                body.append("source /group/mwaops/PULSAR/psrBash.profile")
                                body.append("module swap craype-ivybridge craype-sandybridge")
        
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
        
                                        body.append(" aprun -n 1 -N 1 {0} -o {1}/{2} -s {3} -r {4} -i {5} -f 128 -n {6} -c {7:0>2} -d {8}".format("mwac_offline",corr_dir,obsid,unix_time,num_frames,integrations,int(ft_res[0]/10),gpubox_label,file))
                                        to_corr += 1
                                        # with open(corr_batch, 'a') as batch_file:
                                        #     batch_file.write(corr_line)
                                        #     to_corr = to_corr+1
        
                                secs_to_run = str(datetime.timedelta(seconds=2*12*num_frames*to_corr)) # added factor two on 10 April 2017 as galaxy seemed really slow...
                                submit_slurm(corr_batch,body,slurm_kwargs={"time" : secs_to_run, "partition" : "gpuq"}, batch_dir=batch_dir)
                                # batch_submit_line = "sbatch --workdir={0} --time={1} --partition=gpuq --gid=mwaops {2} \n".format(corr_dir,secs_to_run,corr_batch)
                                # submit_cmd = subprocess.Popen(batch_submit_line,shell=True,stdout=subprocess.PIPE)
                                # jobid=""
                                # for line in submit_cmd.stdout:
                                #     if "Submitted" in line:
                                #         (word1,word2,word3,jobid) = line.split()
                        else:
                                print "Couldn't find any recombine files. Aborting here."


def write_rts_in_files(chan_groups,basepath,rts_in_file,chan_type,cal_obs_id,count=0):
    import re
    import numpy as np
    chan_file_dict = {} # used to keep track of which file needs how many nodes
    
    cc = 0
    for c in chan_groups:
        if len(c)==1:
            #single channel group, write its own rts_in file
            subid = str(count+1)

            # NOTE: we've already re-ordered the channels respectively, so no need to do any weird counting
            #if chan_type == "low":
            #    subid = str(count+1)
            #elif chan_type == "high":
            #        subid = str(count+1)
            #else:
            #        print "No channel group given, assuming low"
            #        subid = str(count+1)
         
            if chan_type == "low":
                offset = cc
            elif chan_type == "high":
                offset = 24-int(subid)
            else:
                print "Invalid channel group type: must be \"low\" or \"high\". Aborting!"
                sys.exit(1)

            basefreq = 1.28*(c[0]-offset)-0.625

            # use re.compile to make search expressions
            with open(rts_in_file,'rb') as f:
                string = f.read()
                # find the pattern and select the entire line for replacement
                string = re.sub("ObservationFrequencyBase=.*\n","ObservationFrequencyBase={0}\n".format(basefreq),string)
                # include the SubBandIDs tag 
                string = re.sub("StartProcessingAt=0\n","StartProcessingAt=0\nSubBandIDs={0}\n\n".format(subid),string)

            fname = "{0}/rts_{1}_chan{2}.in".format(basepath,cal_obs_id,c[0])
            chan_file_dict[fname] = 1 # this particular rts_in file has only 1 channel
            with open(fname,'wb') as f:
                f.write(string)
            
            print "Single channel:: (subband id, abs. chan, abs. freq) = ({0}, {1}, {2}) {3}".format(subid,c[0],basefreq,cc)

            count += 1
            cc += 1
        elif len(c)>1:
            # multiple consecutive channels
            subids = [str(count+i+1) for i in range(len(c))]
            
            #if chan_type == "low":
            #    subids = [str(count+i+1) for i in range(len(c))]
            #elif chan_type == "high":
            #        subids = [str(count+i+1) for i in range(len(c))]
            #else:
            #        print "No channel group given, assuming low"
            #    subid = [str(count+i+1) for i in range(len(c))]
        
            if chan_type == "low":
                offset = cc
            elif chan_type == "high":
                offset = np.array([24-int(x) for x in subids])
            else:
                print "Invalid channel group type: must be \"low\" or \"high\". Aborting!"
                sys.exit(1)
  
            freqs = 1.28*(np.array(c)-offset)-0.625
            basefreq = min(freqs)

            # use re.compile to make search expressions
            with open(rts_in_file,'rb') as f:
                string = f.read()
                # find the pattern and select the entire line for replacement
                string = re.sub("ObservationFrequencyBase=.*\n","ObservationFrequencyBase={0}\n".format(basefreq),string)
                # include the SubBandIDs tag
                string = re.sub("StartProcessingAt=0\n","StartProcessingAt=0\nSubBandIDs={0}\n\n".format(",".join(subids)),string)
            
            print "Multiple channels:: (subband ids, abs. chans, abs. freqs) = {0}".format(", ".join("({0}, {1}, {2})".format(i,j,k) for i,j,k in zip(subids,c,freqs)))

            if chan_type == "low":
                fname = "{0}/rts_{1}_chan{2}-{3}.in".format(basepath,cal_obs_id,min(c),max(c))
            elif chan_type == "high":
                fname = "{0}/rts_{1}_chan{2}-{3}.in".format(basepath,cal_obs_id,max(c),min(c))
            else:
                print "No channel group given, assuming low"
                fname = "{0}/rts_{1}_chan{2}-{3}.in".format(basepath,cal_obs_id,min(c),max(c))

            chan_file_dict[fname] = len(c)
            with open(fname,'wb') as f:
                f.write(string)

            count += len(c)
            cc += len(c)
        else:
            print "Reached a channel group with no entries!? Aborting."
            sys.exit(1)


    return chan_file_dict,count


def run_rts(obs_id, cal_obs_id, product_dir, rts_in_file, args, rts_output_dir=None):
    rts_run_file = distutils.spawn.find_executable('run_rts.sh')
    vcs_database_id = database_command(args, obs_id)
    #[BWM] Re-written to incorporate picket-fence mode of calibration (21/02/2017)
    
    batch_dir = product_dir+"/batch"
    if rts_output_dir:
        product_dir = os.path.abspath(rts_output_dir)
    else:
        # as with the other functions product_dir should come as /group/mwaops/obs_id
        product_dir = "{0}/{1}/{2}/{3}".format(product_dir,'cal', cal_obs_id,'rts')
    mdir(product_dir,'RTS output')
  
    # check if the meta-file has already been created, in which case don't query the database again
    metafile = "{0}/{1}.meta".format(product_dir,cal_obs_id)
    metafile_exists = False
    if os.path.isfile(metafile):
        print "Found observation metafile: {0}".format(metafile)
        channels = None
        with open(metafile,'rb') as m:
            for line in m.readlines():
                if line.startswith("channels"):
                    channels = line.strip().split(",")[1:]
        if channels == None :
            print "Channels keyword not found in metafile. Re-querying the database."
        else:
            metafile_exists = True
            channels = [int(c) for c in channels]
    
    if metafile_exists == False:
        print "Querying the database for calibrator obs ID {0}...".format(cal_obs_id)
        obs_info = getmeta(service='obs', params={'obs_id':str(cal_obs_id)})
        channels = obs_info[u'rfstreams'][u"0"][u'frequencies']
        with open(metafile,"wb") as m:
            m.write("#Metadata for obs ID {0} required to determine if: normal or picket-fence\n".format(cal_obs_id))
            m.write("channels,{0}".format(",".join([str(c) for c in channels])))


    #define some of the header/body text to go into the slurm submission script
    body = []
    body.append("module load cudatoolkit")
    body.append("module load cfitsio")

    body.append("cd {0}".format(product_dir))

    # from the channels, first figure out if they are all consecutive
    if channels[-1]-channels[0] == len(channels)-1: #TODO: this assumes also ascending order: is that always true??
        # the channels are consecutive and this is a normal observation
        # submit the jobs as normal
        print "The channels are consecutive: this is a normal observation"
        nnodes = 25 
        rts_batch = "RTS_{0}".format(cal_obs_id)
        slurm_kwargs = {"partition":"gpuq", "workdir":"{0}".format(product_dir), "time":"00:20:00", "nodes":"{0}".format(nnodes)}
        commands = list(body) # make a copy of body to then extend
        commands.append("aprun -n {0} -N 1  rts_gpu {1}".format(nnodes,rts_in_file))
        submit_slurm(rts_batch, commands, slurm_kwargs=slurm_kwargs, batch_dir=batch_dir,submit=True)
    else:
        # it is a picket fence observation, we need to do some magic
        # will use the given rts_in file as a base format
        # ** THIS ASSUMES THAT THE ORIGINAL RTS_IN FILE HAS BEEN EDITTED TO APPROPRIATELY REPRESENT THE CORRELATOR SETUP **
        #    e.g. the channel width, number of channels, correlator dump times, etc.
        print "The channels are NOT consecutive: this is a picket-fence observation"
        # figure out whether channels are "high" or "low"
        hichans = [c for c in channels if c>128]
        lochans = [c for c in channels if c<=128]
        print "High channels:",hichans
        print "Low channels:",lochans

        # get the consecutive chunks within each channel subselection
        print "Grouping individual channels into consecutive chunks (if possible)"
        from itertools import groupby
        from operator import itemgetter
        # create a list of lists with consecutive channels grouped wihtin a list
        hichan_groups = [map(itemgetter(1),g)[::-1] for k,g in groupby(enumerate(hichans), lambda (i, x): i-x)][::-1] # reversed order (both of internal lists and globally)
        lochan_groups = [map(itemgetter(1),g) for k,g in groupby(enumerate(lochans), lambda (i, x): i-x)]
        print "High channels (grouped):",hichan_groups
        print "Low channels (grouped):",lochan_groups

        # for each group of channels, we need to write 1 rts_in file
        print "Mapping GPU box numbers to coarse channels..."
        basepath = os.path.dirname(rts_in_file)

        # write out the RTS in files and keep track of the number of nodes required for each
        count = 0
        lodict,count = write_rts_in_files(lochan_groups,basepath,rts_in_file,"low",cal_obs_id,count)
        hidict,count = write_rts_in_files(hichan_groups,basepath,rts_in_file,"high",cal_obs_id,count)
        chan_file_dict = lodict.copy()
        chan_file_dict.update(hidict)

        #print "Editting run_rts.sh script for group sizes... (creating temporary copies)"
        # the new, channel-based rts_in files should now be in the same location as the original rts_in file
        # we need to now adjust the run_rts.sh script to work for each of the groups
        lolengths = set([len(l) for l in lochan_groups]) # figure out the unique lengths 
        hilengths = set([len(h) for h in hichan_groups])
        lengths = lolengths.union(hilengths) # combine the sets
    
        # Now submit the RTS jobs
        print "Writing and submitting RTS jobs"
        
        for k,v in chan_file_dict.iteritems():
            nnodes = v + 1
            channels = k.split('_')[-1].split(".")[0]
            rts_batch = "RTS_{0}_{1}".format(cal_obs_id,channels)
            slurm_kwargs = {"partition":"gpuq", "workdir":"{0}".format(product_dir), "time":"00:45:00", "nodes":"{0}".format(nnodes)}
            commands = list(body) # make a copy of body to then extend
            commands.append("aprun -n {0} -N 1  rts_gpu {1}".format(nnodes,k))
            submit_slurm(rts_batch, commands, slurm_kwargs=slurm_kwargs, batch_dir=batch_dir,submit=True)

                        

def coherent_beam(obs_id, start, stop, execpath, data_dir, product_dir, metafile, nfine_chan, pointing,
                 args, rts_flag_file=None, bf_format=' -f psrfits_header.txt', DI_dir=None, calibration_type='rts'):
    vcs_database_id = database_command(args, obs_id)
    # Print relevant version numbers to screen
    mwacutils_version_cmd = "{0}/make_beam -V".format(execpath)
    mwacutils_version = subprocess.Popen(mwacutils_version_cmd, stdout=subprocess.PIPE, shell=True).communicate()[0]
    tested_version  = "0.9.0"
    print "Current version of MWACUtils = {0}".format(mwacutils_version.strip())
    print "Tested  version of MWACUtils = {0}".format(tested_version.strip())

    # Need to run get_delays and then the beamformer on each desired coarse channel
    if not DI_dir:
        print "You need to specify the path to the calibrator files, either where the DIJs are or where the offringe calibration_solution.bin file is. Aborting here."
        quit()
    DI_dir = os.path.abspath(DI_dir)
    RA = pointing[0]
    Dec = pointing[1]

    # get_delays requires the start time in UTC, get it from the start GPS time
    # as is done in timeconvert.py
    t=ephem_utils.MWATime(gpstime=float(start))
    utctime = t.strftime('%Y-%m-%dT%H:%M:%S %Z')[:-4]

    print "Running get_delays"
    P_dir = product_dir+"/pointings"
    mdir(P_dir, "Pointings")
    pointing_dir = "{0}/{1}_{2}".format(P_dir, RA, Dec)
    mdir(pointing_dir, "Pointing {0} {1}".format(RA, Dec))

    startjobs = True
    chan_index = 0
    get_delays_batch = "{0}/batch/gd_{1}_{2}.batch".format(product_dir,start, stop)
    bf_adjust_flags = distutils.spawn.find_executable("bf_adjust_flags.py")
    with open(get_delays_batch,'w') as batch_file:
        batch_line = "#!/bin/bash -l\n#SBATCH --export=NONE\n#SBATCH --account=mwaops\n#SBATCH --output={0}/batch/gd_{1}_{2}.out\n#SBATCH --mail-type=ALL\n".format(product_dir,start,stop)
        batch_file.write(batch_line)
        batch_file.write('source /group/mwaops/PULSAR/psrBash.profile\n')
        batch_file.write('module swap craype-ivybridge craype-sandybridge\n')
        for gpubox in ["{0:0>2}".format(i) for i in range(1,25)]:
            #DI_file = "{0}/{1}".format(DI_dir, ?) # Need to finish file path
            pointing_chan_dir = "{0}/{1}".format(pointing_dir,gpubox)
            mdir(pointing_chan_dir, "Pointing {0} {1} gpubox {2}".format(RA, Dec, gpubox))

            if calibration_type == 'rts':
                chan_list = get_frequencies(metafits_file, resort=True)
                DI_file = "{0}/DI_JonesMatrices_node0{1}.dat".format(DI_dir, gpubox)
                jones_option = "-R {0}".format(DI_file)
            elif calibration_type == 'offringa':
                chan_list = get_frequencies(metafits_file, resort=False)
                DI_file = "{0}/calibration_solution.bin".format(DI_dir)
                jones_option = "-O {0} -C {1}".format(DI_file, int(gpubox)-1)
            channel_file = "{0}/channel".format(pointing_chan_dir)
            with open(channel_file,"w") as ch_file:
                ch_line = "{0}".format(chan_list[chan_index]);
                ch_file.write(ch_line)

            #ASSUMES 10kHz channels <beware>

            basefreq = int(chan_list[chan_index]) * 1.28e6 - 5e3 - 640e3  + 5e3
        
            if (os.path.isfile(DI_file)):
                ch_dir_line = "cd {0}\n".format(pointing_chan_dir)
                batch_file.write(ch_dir_line)
                delays_line = "{0}/get_delays -a {1} -b {2} {3} -m {4} -c -i -p -z {5} -o {6} -f {7} -n {8} -w 10000 -r {9} -d {10}\n".format(execpath, pointing_chan_dir, stop-start+1, jones_option, metafile, utctime, obs_id, basefreq, nfine_chan, RA, Dec) 
                batch_file.write(delays_line)
                if rts_flag_file:
                    flags_file = "{0}/flags.txt".format(pointing_chan_dir)
                    flag_line="{0} {1} {2}\n".format(bf_adjust_flags, rts_flag_file, flags_file)
                    batch_file.write(flag_line)
            else:
                print "WARNING: No Calibration Found for Channel {0}! Could not find file {1}".format(gpubox, DI_file)
                startjobs = False

            chan_index = chan_index+1
    submit_line = "sbatch --account=mwaops --time={0} --workdir={1} --partition=gpuq --gid=mwaops --mail-user={2} {3}\n".format("00:45:00", pointing_dir, e_mail, get_delays_batch)
    print submit_line;
    if startjobs:
        output = subprocess.Popen(submit_line, stdout=subprocess.PIPE, shell=True).communicate()[0]
        # output is something like this: Submitted batch job 245521 on cluster zeus
        # as the beamformer can only run once delays are computed need to put in dependency on get_delay jobs.
        dependsOn = output.split(" ")[3].strip()
        #submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        print "Forming coherent beam. \n"
        print "starts once Job {0} is done.\n".format(dependsOn)
    else:
        print "Not submitted. \n"
 
 
    # Run make_beam
    seconds_to_run = 60*(stop-start+1)
    if seconds_to_run > 86399.:
        secs_to_run = datetime.timedelta(seconds=86399)
    else:
        secs_to_run = datetime.timedelta(seconds=60*(stop-start+1))

    # VDIF needs gpu for inverting the pfb, otherwise cpu nodes are fine
    if bf_format == " -v psrfits_header.txt" or bf_format == " -v psrfits_header.txt -f psrfits_header.txt":
        partition = "gpuq"
        n_omp_threads = 8
    else:
        partition = "workq"
        n_omp_threads = 20
    openmp_line = "export OMP_NUM_THREADS={0}\n".format(n_omp_threads)

    # Run one coarse channel per node
    for coarse_chan in range(24):
        make_beam_batch = "{0}/batch/mb_{1}_{2}_ch{3}.batch".format(product_dir, RA, Dec, coarse_chan)
        make_beam_batch_out = make_beam_batch.replace('.batch','.out')
        with open(make_beam_batch, 'w') as batch_file:
            batch_file.write("#!/bin/bash -l\n")
            nodes_line = "#SBATCH --nodes=1\n#SBATCH --export=NONE\n#SBATCH --account=mwaops\n" 
            batch_file.write(nodes_line)
            output_line = "#SBATCH --output={0}\n".format(make_beam_batch_out)
            batch_file.write(output_line)
            time_line = "#SBATCH --time=%s\n" % (str(secs_to_run))
            batch_file.write(time_line)
            batch_file.write('source /group/mwaops/PULSAR/psrBash.profile\n')
            batch_file.write('module swap craype-ivybridge craype-sandybridge\n')
            # The beamformer runs on all files within time range specified with
            # the -b and -e flags
            batch_file.write(openmp_line)
            aprun_line = "aprun -n 1 -d {0} {1}/make_beam -o {2} -b {3} -e {4} -a 128 -n 128 -N {5} -t 1 {6} -c phases.txt -w flags.txt -d {7}/combined -D {8}/ {9} \n".format(n_omp_threads, execpath, obs_id, start, stop, coarse_chan, jones, data_dir, pointing_dir, bf_format)
            batch_file.write(aprun_line)
        
        submit_line = "sbatch --workdir={0} --partition={1} -d afterok:{2} --gid=mwaops --mail-user={3} {4} \n".format(pointing_dir, partition, dependsOn, e_mail, make_beam_batch)
        print submit_line
        if startjobs:
            output = subprocess.Popen(submit_line, stdout=subprocess.PIPE, shell=True).communicate()[0]
            jobID = output.split(" ")[3].strip()
                #submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
            print "Submitted as job {0}".format(jobID)
        else:
            print "Not submitted. \n"


def coherent_beam_new(obs_id, start, stop, execpath, data_dir, product_dir, batch_dir, metafits_file, nfine_chan, pointing,
                      args, rts_flag_file=None, bf_format=' -f psrfits_header.txt', DI_dir=None,
                      calibration_type='rts'):
    """
    This function runs the new version of the beamformer. It is modelled after the old function above and will likely
    be able to be streamlined after working implementation (SET)
    """
    vcs_database_id = database_command(args, obs_id)  # why is this being calculated here? (SET)
    # Print relevant version numbers to screen
    make_beam_version_cmd = "{0}/make_beam_small -V".format(execpath)
    make_beam_version = subprocess.Popen(make_beam_version_cmd, stdout=subprocess.PIPE, shell=True).communicate()[0]
    tested_version = "?.?.?"
    print "Current version of make_beam_small = {0}".format(make_beam_version.strip())
    print "Tested version of make_beam_small = {0}".format(tested_version.strip())

    metafile = "{0}/{1}.meta".format(product_dir, obs_id)
    metafile_exists = False
    if os.path.isfile(metafile):
        print "Found observation metafile: {0}".format(metafile)
        channels = None
        with open(metafile, 'rb') as m:
            for line in m.readlines():
                if line.startswith("channels"):
                    channels = line.strip().split(",")[1:]
        if channels == None:
            print "Channels keyword not found in metafile. Re-querying the database."
        else:
            metafile_exists = True
            channels = [int(c) for c in channels]
    # Grabbing this from the calibrate section for now. This should be streamlined to call an external function (SET)
    if metafile_exists == False:
        print "Querying the database for calibrator obs ID {0}...".format(obs_id)
        obs_info = getmeta(service='obs', params={'obs_id': str(obs_id)})
        channels = obs_info[u'rfstreams'][u"0"][u'frequencies']
        with open(metafile, "wb") as m:
            m.write("#Metadata for obs ID {0} required to determine if: normal or picket-fence\n".format(obs_id))
            m.write("channels,{0}".format(",".join([str(c) for c in channels])))
    hichans = [c for c in channels if c>128]
    lochans = [c for c in channels if c<=128]
    lochans.extend(list(reversed(hichans)))
    ordered_channels = lochans

    # Run for each coarse channel. Calculates delays and makes beam

    if not DI_dir:
        print "You need to specify the path to the calibrator files, either where the DIJs are or where the\
        Offringa calibration_solutions.bin file is. Aborting here"
        quit()
    DI_dir = os.path.abspath(DI_dir)
    RA = pointing[0]
    Dec = pointing[1]

    # make_beam_small requires the start time in UTC, get it from the start
    # GPS time as is done in timeconvert.py
    t = ephem_utils.MWATime(gpstime=float(start))
    utctime = t.strftime('%Y-%m-%dT%H:%M:%S %Z')[:-4]

    print "Running make_beam_small"
    P_dir = product_dir+"/pointings"
    mdir(P_dir, "Pointings")
    pointing_dir = "{0}/{1}_{2}".format(P_dir, RA, Dec)
    mdir(pointing_dir, "Pointing {0} {1}".format(RA, Dec))
    # startjobs = True


    # set up SLURM requirements
    seconds_to_run = 60 * (stop - start + 1)  # This should be able to be reduced after functional testing
    if seconds_to_run > 86399.:
        secs_to_run = datetime.timedelta(seconds=86399)
    else:
        secs_to_run = datetime.timedelta(seconds=seconds_to_run)

    # VDIF will need gpuq if inverting pfb with '-m' option, otherwise cpuq is fine
    # In general this needs to be cleaned up, prefferably to be able to intelligently select a
    # queue and a maximum wall time (SET)
    partition = "gpuq"
    n_omp_threads = 8
    openmp_line = "export OMP_NUM_THREADS={0}".format(n_omp_threads)

    # Run one coarse channel per node
    #for coarse_chan in range(24):
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
            print "Please an accepted calibratin type. Aborting here."
            quit()

        make_beam_small_batch = "mbs_{0}_{1}_ch{2}".format(RA, Dec, coarse_chan)
        commands = []
        commands.append("source /group/mwaops/PULSAR/psrBash.profile")
        commands.append("module swap craype-ivybridge craype-sandybridge")
        commands.append(openmp_line)
        commands.append("cd {0}".format(pointing_dir))
        commands.append("aprun -n 1 -d {0} {1}/make_beam_small -o {2} -b {3} -e {4} -a 128 -n 128 -f {5} {6} -d "
                        "{7}/combined -R {8} -D {9} -r 10000 -m {10} -z {11}".format(n_omp_threads, execpath, obs_id, start,
                        stop, coarse_chan, jones_option, data_dir, RA, Dec, metafits_file, utctime))  # check these
        submit_slurm(make_beam_small_batch, commands, batch_dir=batch_dir,
                     slurm_kwargs={"time": secs_to_run, "partition": partition}, submit=True)

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

    modes=['download', 'download_ics', 'download_cal', 'recombine','correlate','calibrate', 'beamform']
    bf_out_modes=['psrfits', 'vdif', 'both']
    jobs_per_node = 8
    chan_list_full=["ch01","ch02","ch03","ch04","ch05","ch06","ch07","ch08","ch09","ch10","ch11","ch12","ch13","ch14","ch15","ch16","ch17","ch18","ch19","ch20","ch21","ch22","ch23","ch24"]
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

    group_calibrate = OptionGroup(parser, 'Calibration Options')
    group_calibrate.add_option('--rts_in_file', type='string', help="Either relative or absolute path (including file name) to setup file for the RTS.", default=None)
    group_calibrate.add_option('--rts_output_dir', help="Working directory for RTS -- all RTS output files will end up here. "+\
                                   "Default is /group/mwaops/[obsID]/cal/cal_obsID/.", default=None)

    group_beamform = OptionGroup(parser, 'Beamforming Options')
    group_beamform.add_option("-p", "--pointing", nargs=2, help="required, R.A. and Dec. of pointing, e.g. \"19:23:48.53\" \"-20:31:52.95\"")
    group_beamform.add_option("--DI_dir", default=None, help="Directory containing either Direction Independent Jones Matrices (as created by the RTS) " +\
                                  "or calibration_solution.bin as created by Andre Offringa's tools.[no default]")
    group_beamform.add_option("--bf_new", action="store_true", default=False, help="Run the new beamformer (make_beam_small). [default=%defalut]")
    group_beamform.add_option("--bf_out_format", type="choice", choices=['psrfits','vdif','both'], help="Beam former output format. Choices are {0}. [default=%default]".format(bf_out_modes), default='psrfits')
    group_beamform.add_option("--flagged_tiles", type="string", default=None, help="Path (including file name) to file containing the flagged tiles as used in the RTS, will be used to adjust flags.txt as output by get_delays. [default=%default]")
    group_beamform.add_option('--cal_type', type='string', help="Use either RTS (\"rts\") solutions or Andre-Offringa-style (\"offringa\") solutions. Default is \"rts\". If using Offringa's tools, the filename of calibration solution must be \"calibration_solution.bin\".", default="rts")
    group_beamform.add_option("-E", "--execpath", type="string", default='/group/mwaops/PULSAR/bin/', help=SUPPRESS_HELP)

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
    parser.add_option_group(group_download)
    parser.add_option_group(group_correlate)
    parser.add_option_group(group_calibrate)
    parser.add_option_group(group_beamform)

    (opts, args) = parser.parse_args()
    
    if opts.all and (opts.begin or opts.end):
        print "Please specify EITHER (-b,-e) OR -a"
        quit()
    elif opts.all:
        opts.begin, opts.end = obs_max_min(opts.cal_obs if opts.mode == 'download_cal' else opts.obs)
    # make sure we can process increments smaller than 64 seconds when not in calibration related mode
    if opts.mode not in ['download_cal','calibrate']:
        if opts.end - opts.begin +1 < opts.increment:
            opts.increment = opts.end - opts.begin + 1
    e_mail = ""
    if opts.mail:
        e_mail = get_user_email()
        print "Sending info to {0}".format(e_mail)
    if not opts.mode:
      print "Mode required {0}. Please specify with -m or --mode.".format(modes)
      quit()
    if not opts.obs:
        print "Observation ID required, please put in with -o or --obs"
        quit()
    if opts.begin > opts.end:
        print "Starting time is after end time"
        quit()
    if opts.mode == "beamform":
        if not opts.pointing:
            print "Pointing (-p) required in beamformer mode"
            quit()
        if (opts.bf_out_format == 'psrfits'):
            bf_format = " -f psrfits_header.txt"
        elif  (opts.bf_out_format == 'vdif'):
            bf_format = " -v psrfits_header.txt"
        elif (opts.bf_out_format == 'both'):
            bf_format = " -v psrfits_header.txt -f psrfits_header.txt"
            print "Writing out both psrfits and vdif."

        if opts.execpath:
            execpath = opts.execpath
    if opts.work_dir:
        print "YOU ARE MESSING WITH THE DEFAULT DIRECTORY STRUCTURE FOR PROCESSING -- BE SURE YOU KNOW WHAT YOU ARE DOING!"
        time.sleep(5)
        data_dir = product_dir = "{0}/{1}".format(opts.work_dir, opts.obs)
    else:
        data_dir = '/astro/mwaops/vcs/{0}'.format(opts.obs)
        product_dir = '/group/mwaops/vcs/{0}'.format(opts.obs)
    batch_dir = "{0}/batch".format(product_dir)
    mdir(data_dir, "Data")
    mdir(product_dir, "Products")
    mdir(batch_dir, "Batch")
    metafits_file = "{0}/{1}_metafits_ppds.fits".format(data_dir, opts.obs)
    # TODO: modify metafits downloader to not just do a trivial wget 

    print "Processing Obs ID {0} from GPS times {1} till {2}".format(opts.obs, opts.begin, opts.end)

    if opts.mode == 'download_ics':
        print opts.mode
        vcs_download(opts.obs, opts.begin, opts.end, opts.increment, opts.head, data_dir, product_dir, opts.parallel_dl, sys.argv, ics=True)
    elif opts.mode == 'download':
        print opts.mode
        vcs_download(opts.obs, opts.begin, opts.end, opts.increment, opts.head, data_dir, product_dir, opts.parallel_dl, sys.argv, n_untar=opts.untar_jobs, keep='-k' if opts.keep_tarball else "")
    elif opts.mode == 'recombine':
        print opts.mode
        ensure_metafits(data_dir, opts.obs, metafits_file)
        vcs_recombine(opts.obs, opts.begin, opts.end, opts.increment, data_dir, product_dir, sys.argv)
    elif opts.mode == 'correlate':
        print opts.mode 
        ensure_metafits(data_dir, opts.obs, metafits_file)
        vcs_correlate(opts.obs, opts.begin, opts.end, opts.increment, data_dir, product_dir, opts.ft_res, sys.argv, metafits_file)
    elif opts.mode == 'download_cal':
        print opts.mode
        if not opts.cal_obs:
            print "You need to also pass the calibrator observation ID. Aborting here."
            quit()
        if opts.cal_obs == opts.obs:
            print "The calibrator obsID cannot be the same as the target obsID -- there are not gpubox files for VCS data on the archive." 
            quit()
        data_dir = data_dir.replace(str(opts.obs), str(opts.cal_obs))
        mdir(data_dir, "Calibrator Data")
        download_cal(opts.obs, opts.cal_obs, data_dir, product_dir, sys.argv, opts.head)
    elif opts.mode == 'calibrate':
        print opts.mode
        if not opts.rts_in_file:
            print "You have to provide the full path to the setup file for the RTS. Aborting here."
            quit()
        if not os.path.isfile(opts.rts_in_file):
            print "Your are not pointing at a file with your input to --rts_in_file. Aboring here as the RTS will not run..."
            quit()
        if not opts.cal_obs:
            print "You need to also pass the calibrator observation ID (GPS seconds), otherwise we can't query the database. Aborting here."
            quit()
        # turn whatever path we got into an absolute path 
        rts_in_file = os.path.abspath(opts.rts_in_file)
        if opts.rts_output_dir:
            rts_output_dir = os.path.abspath(opts.rts_output_dir)
        run_rts(opts.obs, opts.cal_obs, product_dir, rts_in_file, sys.argv, opts.rts_output_dir)
    elif opts.mode == 'beamform':
        print opts.mode
        if not opts.DI_dir:
            print "You need to specify the path to either where the DIJs are or where the offringe calibration_solution.bin file is. Aborting here."
            quit()
        if opts.flagged_tiles:
            flagged_tiles_file = os.path.abspath(opts.flagged_tiles)
            if not os.path.isfile(opts.flagged_tiles):
                print "Your are not pointing at a file with your input to --flagged_tiles. Aboring here as the beamformer will not run..."
                quit()
        else:
            flagged_tiles_file = None
        ensure_metafits(data_dir, opts.obs, metafits_file)
        from mwapy import ephem_utils
        if opts.bf_new:
            coherent_beam_new(opts.obs, opts.begin, opts.end, opts.execpath, data_dir, product_dir, batch_dir, metafits_file,
                      opts.nfine_chan, opts.pointing, sys.argv, flagged_tiles_file, bf_format, opts.DI_dir, opts.cal_type)
        else:
            coherent_beam(opts.obs, opts.begin, opts.end, opts.execpath, data_dir, product_dir, metafits_file,
                      opts.nfine_chan, opts.pointing, sys.argv, flagged_tiles_file, bf_format, opts.DI_dir, opts.cal_type)
    else:
        print "Somehow your non-standard mode snuck through. Try again with one of {0}".format(modes)
        quit()


