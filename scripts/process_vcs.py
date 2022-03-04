#!/usr/bin/env python3

import subprocess
import os
import sys
import datetime
from time import sleep
from astropy.io import fits as pyfits
from astropy.time import Time
import numpy as np
import logging

#vcstools functions
from vcstools.job_submit import submit_slurm
import vcstools.metadb_utils as meta
from vcstools.general_utils import mdir, gps_to_utc, create_link
from vcstools.pointing_utils import format_ra_dec
from vcstools.config import load_config_file
from vcstools.general_utils import sfreq, setup_logger

logger = logging.getLogger(__name__)


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def get_user_email():
    command="echo `ldapsearch -x \"uid=$USER\" mail |grep \"^mail\"|cut -f2 -d' '`"
    email = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True).communicate()[0]
    return email.strip()


def gps_time_lists(start, stop, chunk):
    time_chunks = []
    while (start + chunk) < stop:
        time_chunks.append((start, start + chunk - 1))
        start += chunk
    time_chunks.append((start, stop))
    return time_chunks


def get_frequencies(metafits,resort=False):
    # TODO: for robustness, this should force the entries to be 3-digit numbers
    hdulist    = pyfits.open(metafits)
    freq_str   = hdulist[0].header['CHANNELS']
    freq_array = [int(f) for f in freq_str.split(',')]
    if resort:
        return sfreq(freq_array)
    else:
        return freq_array


def vcs_download(obsid, start_time, stop_time, increment, data_dir,
                 product_dir, parallel,
                 ics=False, n_untar=2, keep="", vcstools_version="master",
                 nice=0):

    #Load computer dependant config file
    comp_config = load_config_file()

    logger.info("Downloading files from archive")
    voltdownload = "voltdownload.py"
    obsinfo = meta.getmeta(service='obs', params={'obs_id':str(obsid)})
    comb_del_check = meta.combined_deleted_check(obsid, begin=start_time, end=stop_time)
    data_format = obsinfo['dataquality']
    if data_format == 1 or (comb_del_check and data_format == 6):
        # either only the raw data is available (data_format == 1) 
        # or there was combined files but they were deleted (comb_del_check and data_format == 6)
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
    mdir(dl_dir, dir_description, gid=comp_config['gid'])
    create_link(data_dir, target_dir, product_dir, link)
    batch_dir = product_dir+"/batch/"

    for time_to_get in range(start_time,stop_time,increment):
        if time_to_get + increment > stop_time:
            increment = stop_time - time_to_get + 1
        #need to subtract 1 from increment since voltdownload wants how many
        #seconds PAST the first one

        voltdownload_batch = "volt_{0}".format(time_to_get)
        check_batch = "check_volt_{0}".format(time_to_get)
        volt_secs_to_run = datetime.timedelta(seconds=500*increment)
        check_secs_to_run = "15:00"
        if data_type == 16:
            check_secs_to_run = "10:15:00"

        checks = "checks.py"
        # Write out the checks batch file but don't submit it
        commands = []
        commands.append("newcount=0")
        commands.append("let oldcount=$newcount-1")
        commands.append("sed -i -e \"s/oldcount=${{oldcount}}/oldcount=${{newcount}}/\" {0}".\
                        format(batch_dir+voltdownload_batch+".batch"))
        commands.append("oldcount=$newcount; let newcount=$newcount+1")
        commands.append("sed -i -e \"s/_${{oldcount}}.out/_${{newcount}}.out/\" {0}".\
                        format(batch_dir+voltdownload_batch+".batch"))
        checks_command = "-m download -o {0} -w {1} -b {2} -i {3} --data_type {4}".format(obsid,
                            dl_dir, time_to_get, increment, data_type)
        commands.append('{0} {1}'.format(checks, checks_command))
        commands.append("if [ $? -eq 1 ];then")
        commands.append("sbatch {0}".format(batch_dir+voltdownload_batch+".batch"))
        # if we have tarballs we send the untar jobs to the workq
        if data_type == 16:
            commands.append("else")
            untar = 'untar.sh'
            untar_command = "-w {0} -o {1} -b {2} -e {3} -j {4} {5}".format(dl_dir,
                                obsid, time_to_get, time_to_get+increment-1, n_untar,
                                keep)
            commands.append('{0} {1}'.format(untar, untar_command))

            #commands.append("sbatch {0}.batch".format(batch_dir+tar_batch))
        commands.append("fi")

        # Download and checks should be done on Zeus's cpuq. This will only work
        # on Galaxy as the Ozstar workflow is different
        submit_slurm(check_batch, commands, batch_dir=batch_dir,
                        slurm_kwargs={"time": check_secs_to_run,
                                    "nice": nice},
                        vcstools_version=vcstools_version, submit=False,
                        outfile=batch_dir+check_batch+"_0.out",
                        queue="zcpuq", export="NONE", mem=10240,
                        # Manually handing it the module dir as it should only run
                        module_dir='/group/mwa/software/modulefiles')

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
        body.append("{0} {1}".format(voltdownload, voltdownload_command))
        submit_slurm(voltdownload_batch, body, batch_dir=batch_dir,
                        module_list=module_list,
                        slurm_kwargs={"time" : str(volt_secs_to_run),
                                      "nice" : nice},
                        vcstools_version=vcstools_version,
                        outfile=batch_dir+voltdownload_batch+"_1.out",
                        queue="copyq", export="NONE", mem=5120,
                        # Manually handing it the module dir as it should only run
                        module_dir='/group/mwa/software/modulefiles')


def download_cal(obs_id, cal_obs_id, data_dir, product_dir,
                 vcstools_version="master", nice=0):

    #Load computer dependant config file
    comp_config = load_config_file()

    batch_dir = os.path.join(product_dir, 'batch')
    product_dir = os.path.join(product_dir, 'cal', str(cal_obs_id))
    vis_dir = os.path.join(data_dir, 'vis')
    mdir(vis_dir, 'Calibrator vis', gid=comp_config['gid'])
    mdir(product_dir, 'Calibrator product', gid=comp_config['gid'])
    mdir(batch_dir, 'Batch', gid=comp_config['gid'])
    # Downloads the visablities to  /astro/mwavcs/vcs/[cal_obs_id]/vis
    # but creates a link to it here /astro/mwavcs/vcs/[obs_id]/cal/[cal_obs_id]
    csvfile = os.path.join(batch_dir, "{0}_dl.csv".format(cal_obs_id))
    create_link(data_dir, 'vis', product_dir, 'vis')
    obsdownload_batch = "caldownload_{0}".format(cal_obs_id)
    secs_to_run = "03:00:00" # sometimes the staging can take a while...
    module_list = ["manta-ray-client/python3"]
    commands = []
    commands.append("csvfile={0}".format(csvfile))
    commands.append('cd {0}'.format(vis_dir))
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
    commands.append('mwa_client --csv={0} --dir={1}'.format(csvfile,vis_dir))
    #commands.append("ln -sfn {0} {1}/{2}".format(data_dir, product_dir, 'vis'))
    commands.append('unzip *.zip')
    submit_slurm(obsdownload_batch, commands, batch_dir=batch_dir,
                    module_list=module_list,
                    slurm_kwargs={"time": secs_to_run, "nice": nice},
                    vcstools_version=vcstools_version, queue="copyq",
                    export="NONE", mem=4096,
                    # Manually handing it the module dir as it should only run
                    module_dir='/group/mwa/software/modulefiles')


def vcs_recombine(obsid, start_time, stop_time, increment, data_dir, product_dir,
                  vcstools_version="master", nice=0):

    #Load computer dependant config file
    comp_config = load_config_file()

    logger.info("Running recombine on files")
    jobs_per_node = 8
    target_dir = link = 'combined'
    mdir(data_dir + '/' + target_dir, 'Combined', gid=comp_config['gid'])
    create_link(data_dir, target_dir, product_dir, link)
    batch_dir = product_dir+"/batch/"
    recombine = "recombine.py"
    checks = "checks.py"
    recombine_binary = "recombine"
    for time_to_get in range(start_time,stop_time,increment):
        process_nsecs = increment if (time_to_get + increment <= stop_time) \
                                  else (stop_time - time_to_get + 1)
        if (jobs_per_node > process_nsecs):
                jobs_per_node = process_nsecs
        nodes = (increment+(-increment%jobs_per_node))//jobs_per_node + 1 # Integer division with ceiling result plus 1 for master node
        recombine_batch = "recombine_{0}".format(time_to_get)
        check_batch = "check_recombine_{0}".format(time_to_get)
        #module_list = ["module switch PrgEnv-cray PrgEnv-gnu", "python/3.6.3", "numpy/1.13.3", "mwa-voltage/master"]
        module_list = ["mwa-voltage/master"]
        commands = []
        commands.append("newcount=0")
        commands.append("let oldcount=$newcount-1")
        commands.append("sed -i -e \"s/oldcount=${{oldcount}}/oldcount=${{newcount}}/\" {0}".\
                        format(batch_dir+recombine_batch+".batch"))
        commands.append("oldcount=$newcount; let newcount=$newcount+1")
        commands.append("sed -i -e \"s/_${{oldcount}}.out/_${{newcount}}.out/\" {0}".\
                        format(batch_dir+recombine_batch+".batch"))
        checks_command = "-m recombine -o {0} -w {1}/combined/ -b {2} -i {3}".format(obsid,
                          data_dir, time_to_get, process_nsecs)
        commands.append("{0} {1}".format(checks, checks_command))
        commands.append("if [ $? -eq 1 ];then")
        commands.append("sbatch {0}".format(batch_dir+recombine_batch+".batch"))
        commands.append("fi")
        submit_slurm(check_batch, commands, batch_dir=batch_dir,
                     module_list=module_list,
                     slurm_kwargs={"time": "15:00", "nice": nice},
                     vcstools_version=vcstools_version, submit=False,
                     outfile=batch_dir+check_batch+"_0.out",
                     queue='gpuq', export="NONE")

        #module_list = ["module switch PrgEnv-cray PrgEnv-gnu", "python/3.6.3",
        #               "numpy/1.13.3", "mwa-voltage/master", "mpi4py", "cfitsio"]
        module_list = ["mwa-voltage/master", "mpi4py"]
        commands = []
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
        commands.append("srun --export=all {0} {1}".format(recombine, recombine_command))

        submit_slurm(recombine_batch, commands, batch_dir=batch_dir,
                     module_list=module_list,
                     slurm_kwargs={"time": "06:00:00", "nodes": str(nodes),
                                   "ntasks-per-node": jobs_per_node,
                                   "nice": nice},
                     vcstools_version=vcstools_version,
                     outfile=batch_dir+recombine_batch+"_1.out",
                     queue='gpuq', export="NONE")




def vcs_correlate(obsid,start,stop,increment, data_dir, product_dir, ft_res,
                  metafits, vcstools_version="master", nice=0):

    #Load computer dependant config file
    comp_config = load_config_file()

    logger.info("Correlating files at {0} kHz and {1} milliseconds".\
                format(ft_res[0], ft_res[1]))

    batch_dir = product_dir+"/batch/"
    target_dir = link = 'vis'

    if data_dir == product_dir:
        corr_dir = "{0}/cal/{1}/{2}".format(product_dir, obsid, target_dir)
    else:
        corr_dir = "{0}/{1}".format(data_dir, target_dir)
        product_dir = "{0}/cal/{1}/".format(product_dir, obsid)
        mdir(product_dir, "Correlator", gid=comp_config['gid'])
    mdir(corr_dir, "Correlator Product", gid=comp_config['gid'])
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
                to_corr = 0
                for file in f:
                    (current_time,_) = os.path.splitext(os.path.basename(file))
                    (obsid,gpstime,_) = current_time.split('_')
                    t = Time(int(gpstime), format='gps', scale='utc')
                    unix_time = int(t.unix)

                    offline_correlator_command = "-o {0}/{1} -s {2} -r {3} -i {4} -n {5} "\
                                                 "-c {6:0>2} -d {7}".format(corr_dir, obsid,
                                                 unix_time, num_frames, integrations,
                                                 int(ft_res[0]/10), gpubox_label, file)
                    body.append("{0} {1}".format("offline_correlator", offline_correlator_command))
                    to_corr += 1

                #module_list = ["module switch PrgEnv-cray PrgEnv-gnu"]
                module_list = ["offline_correlator/v1.0.0"]
                secs_to_run = str(datetime.timedelta(seconds=2*12*num_frames*to_corr))
                # added factor two on 10 April 2017 as galaxy seemed really slow...
                submit_slurm(corr_batch, body, module_list=module_list,
                             slurm_kwargs={"time": secs_to_run, "nice": nice},
                             queue='gpuq', vcstools_version=vcstools_version,
                             batch_dir=batch_dir, export="NONE")
            else:
                logger.error("Couldn't find any recombine files. Aborting here.")


def coherent_beam(obs_id, start, stop, data_dir, product_dir, batch_dir,
                  metafits_file, nfine_chan, pointing_list,
                  rts_flag_file=None, bf_formats=None, DI_dir=None,
                  execpath=None, calibration_type='rts', ipfb_filter="LSQ12",
                  vcstools_version="master", nice=0, channels_to_beamform=None,
                  beam_version="FEE2016"):
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
    logger.info("Current version of make_beam = {0}".format(make_beam_version.strip()))

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

    P_dir = os.path.join(product_dir, "pointings")
    mdir(P_dir, "Pointings", gid=comp_config['gid'])
    mdir(os.path.join(product_dir, "incoh"), "Incoh", gid=comp_config['gid'])
    # startjobs = True

    # Set up supercomputer dependant parameters
    import socket
    hostname = socket.gethostname()
    if hostname.startswith('john') or hostname.startswith('farnarkle'):
        max_pointing = 120
    else:
        max_pointing = 15
    if comp_config['ssd_dir'] is None:
        temp_mem = None
    else:
        #Work out required SSD size
        obs_length = stop - start + 1.
        temp_mem = int(0.0012 * obs_length * max_pointing + 1)
        temp_mem_single = int(0.0024 * obs_length + 2)
        if "-s" not in bf_format:
            temp_mem = temp_mem * 4
            temp_mem_single = temp_mem_single *4

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
                mdir("{0}/{1}".format(P_dir, pointing), "Pointing {0}".format(pointing), gid=comp_config['gid'])

            n_omp_threads = 1
            if "v" in bf_formats:
                for pointing in pointing_list:
                    make_beam_small_batch = "mb_{0}_ch{1}".format(pointing, coarse_chan)
                    module_list = [comp_config['container_module']]
                    commands = []
                    commands.append("cd {0}/{1}".format(P_dir,pointing))
                    runline = "srun --export=all -n 1"
                    runline += " -c {}".format(n_omp_threads)
                    if comp_config['container_command'] !='':
                        runline += " {} '".format(comp_config['container_command'])
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
                    if beam_version == "ANALYTIC":
                        runline += " -H"
                    if comp_config['container_command'] !='':
                        runline += "'"
                    commands.append(runline)

                    job_id = submit_slurm(make_beam_small_batch, commands,
                                batch_dir=batch_dir, module_list=module_list,
                                slurm_kwargs={"time":secs_to_run, "nice":nice},
                                queue='gpuq', vcstools_version=vcstools_version,#forces olf version with vdif
                                submit=True, export="NONE", gpu_res=1,
                                cpu_threads=n_omp_threads,
                                mem=comp_config['gpu_beamform_mem'], temp_mem=temp_mem_single)
                    job_id_list.append(job_id)

            else:
                make_beam_small_batch = "mb_{0}_{1}_ch{2}".format(pl, time_now, coarse_chan)
                module_list = [comp_config['container_module']]
                commands = []
                if comp_config['ssd_dir'] is None:
                    # Write outputs to SSDs if on Ozstar
                    commands.append("cd {0}".format(P_dir))
                else:
                    commands.append("cd {0}".format(comp_config['ssd_dir']))

                runline = "srun --export=all -n 1"
                runline += " -c {}".format(n_omp_threads)
                if comp_config['container_command'] !='':
                    runline += " {} '".format(comp_config['container_command'])
                runline += " {}".format(make_beam_cmd)
                runline += " -o {}".format(obs_id)
                runline += " -b {}".format(start)
                runline += " -e {}".format(stop)
                runline += " -a 128"
                runline += " -n 128"
                runline += " -f {}".format(coarse_chan)
                runline += " {}".format(jones_option)
                runline += " -d {}/combined".format(data_dir)
                runline += " -P {}".format(pointing_str)
                runline += " -r 10000"
                runline += " -m {}".format(metafits_file)
                runline += " -z {}".format(utctime)
                runline += " {}".format(bf_formats)
                runline += " -F {}".format(rts_flag_file)
                if beam_version == "ANALYTIC":
                    runline += " -H"
                if comp_config['container_command'] !='':
                    runline += "'"
                commands.append(runline)

                commands.append("")
                if comp_config['ssd_dir'] is not None:
                    for pointing in pointing_list:
                        commands.append("cp {0}/{1}/{2}_{3}_{1}_ch{4}_00*.fits "
                                        "{5}/{1}/".format(comp_config['ssd_dir'], pointing,
                                        project_id, obs_id, coarse_chan, P_dir))
                    if 'i' in bf_formats:
                        commands.append("cp {0}/{1}/{2}_{3}_{1}_ch{4}_00*.fits "
                                        "{5}/{1}/".format(comp_config['ssd_dir'], "incoh",
                                        project_id, obs_id, coarse_chan, product_dir))
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

if __name__ == '__main__':

    # Dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING)

    modes=['download', 'download_ics', 'download_cal', 'recombine', 'correlate', 'beamform']
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
    group_beamform.add_argument("--beam_version", type=str, choices=["ANALYTIC", "FEE2016"], help="The version of the beamformer to use. If FEE2016 is selected but unavailable, analytic will be used", default="FEE2016")


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
    logger = setup_logger(logger, log_level=loglevels[args.loglvl])

    logger.info("Using vcstools/{0}".format(args.vcstools_version))
    if args.version:
        try:
            import version
            logger.info(version.__version__)
            sys.exit(0)
        except ImportError as IE:
            logger.error("Couldn't import version.py - have you installed vcstools?")
            logger.error("ImportError: {0}".format(IE))
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
        product_dir = '{0}{1}'.format(comp_config['base_data_dir'],args.obs)
    batch_dir = "{0}/batch".format(product_dir)
    mdir(data_dir, "Data", gid=comp_config['gid'])
    mdir(product_dir, "Products", gid=comp_config['gid'])
    mdir(batch_dir, "Batch", gid=comp_config['gid'])
    metafits_file = "{0}/{1}_metafits_ppds.fits".format(data_dir, args.obs)
    # TODO: modify metafits downloader to not just do a trivial wget

    logger.info("Processing Obs ID {0} from GPS times {1} till {2}".\
                format(args.obs, args.begin, args.end))

    if args.mode == 'download_ics':
        logger.info("Mode: {0}".format(args.mode))
        vcs_download(args.obs, args.begin, args.end, args.increment,
                     data_dir, product_dir, args.parallel_dl,
                     ics=True, vcstools_version=args.vcstools_version,
                     nice=args.nice)
    elif args.mode == 'download':
        logger.info("Mode: {0}".format(args.mode))
        vcs_download(args.obs, args.begin, args.end, args.increment,
                     data_dir, product_dir, args.parallel_dl,
                     n_untar=args.untar_jobs,
                     keep='-k' if args.keep_tarball else "",
                     vcstools_version=args.vcstools_version, nice=args.nice)
    elif args.mode == 'recombine':
        logger.info("Mode: {0}".format(args.mode))
        meta.ensure_metafits(data_dir, args.obs, metafits_file)
        vcs_recombine(args.obs, args.begin, args.end, args.increment, data_dir,
                      product_dir,
                      vcstools_version=args.vcstools_version, nice=args.nice)
    elif args.mode == 'correlate':
        logger.info("Mode: {0}".format(args.mode))
        meta.ensure_metafits(data_dir, args.obs, metafits_file)
        vcs_correlate(args.obs, args.begin, args.end, args.increment, data_dir,
                      product_dir, args.ft_res, metafits_file,
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
        mdir(data_dir, "Calibrator Data", gid=comp_config['gid'])
        download_cal(args.obs, args.cal_obs, data_dir, product_dir,
                     vcstools_version=args.vcstools_version,
                     nice=args.nice)
    elif args.mode == ('beamform' or 'incoh'):
        logger.info("Mode: {0}".format(args.mode))
        if not args.DI_dir:
            args.DI_dir =  os.path.join(data_dir, 'cal', str(args.cal_obs), 'rts')
            logger.info("No DI_dir input. Assuming DI_dir is {0}".format(args.DI_dir))
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
        meta.ensure_metafits(data_dir, args.obs, metafits_file)
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
        logger.debug("Input pointing list: {}".format(pointing_list))

        # Format input pointing list
        ra_dec_list = [i.split("_") for i in pointing_list]
        formatted_ra_dec_list = format_ra_dec(ra_dec_list)
        pointing_list = ["_".join(i) for i in formatted_ra_dec_list]

        logger.debug("Formatted pointing list: {}".format(pointing_list))

        coherent_beam(args.obs, args.begin, args.end,
                      data_dir, product_dir, batch_dir,
                      metafits_file, args.nfine_chan, pointing_list,
                      rts_flag_file=flagged_tiles_file,
                      bf_formats=bf_format, DI_dir=args.DI_dir,
                      calibration_type=args.cal_type,
                      vcstools_version=args.vcstools_version, nice=args.nice,
                      execpath=args.execpath, ipfb_filter=args.ipfb_filter,
                      beam_version=args.beam_version)
    else:
        logger.error("Somehow your non-standard mode snuck through. "
                     "Try again with one of {0}".format(modes))
        sys.exit(0)
