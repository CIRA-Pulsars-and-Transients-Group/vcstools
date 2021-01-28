import subprocess
from vcstools.config import load_config_file
import logging
import socket

logger = logging.getLogger(__name__)

SLURM_TMPL = """{shebag}

#SBATCH --export={export}
#SBATCH --output={outfile}
{account}
#SBATCH --clusters={cluster}
#SBATCH --partition={partition}
#
#SBATCH --cpus-per-task={threads}
#SBATCH --mem-per-cpu={mem}MB
#SBATCH --nice={nice}
{header}

ncpus={threads}
export OMP_NUM_THREADS={threads}

{switches}
module use {module_dir}
{modules}

{script}
"""
# NOTE: --gid option removed after discussion in helpdesk ticket GS-9370


def submit_slurm(name, commands, tmpl=SLURM_TMPL, slurm_kwargs=None,
                 module_list=[], vcstools_version="master",
                 batch_dir="batch/", depend=None, depend_type='afterok',
                 submit=True, outfile=None, queue="cpuq", export="NONE",
                 gpu_res=None, mem=1024, cpu_threads=1, temp_mem=None,
                 nice=0, shebag='#!/bin/bash -l',
                 module_dir=None, load_vcstools=True):
    """
    Making this function to cleanly submit SLURM jobs using a simple template.

    Parameters
    ----------
    name : str
        The base name that is used to create the "`name`.batch" and "`name`.out" files.

    commands : list of strs
        The actual bash script commands you wnat to run.
        Expects a list where each element is a single line of the bash script.

    tmpl : str
        A template header string with format place holders: export, outfile,
        cluster, header and script.
        This is used to create the final string to be written to the job script.
        For this function, it is required to be SLURM compliant.
        Default: `SLURM_TMPL`

    slurm_kwargs : dict [optional]
        A dictionary of SLURM keyword, value pairs to fill in whatever is not
        in the template supplied to `tmpl`.
        Default: `{}` (empty dictionary, i.e. no additional header parameters)

    module_list : list of str [optional]
        A list of module names (including versions if applicable) that will
        be included in the header for the batch
        scripts. e.g. ["vcstools/master", "mwa-voltage/master", "presto/master"] would append
            module load vcstools/master
            module load mwa-voltage/master
            module load presto/master
        to the header of the batch script. This can also invoke "module use ..." commands.
        NOTE: /group/mwa/software/modulefiles is used and vcstools/master is loaded by default.

    vcstools_version :  str
        The version of vcstools to load. Default: master.

    batch_dir : str [optional]
        The LOCAL directory where you want to write the batch scripts
        (i.e. it will write to `$PWD/batch_dir`).
        Default: "batch/"

    depend : list or None [optional]
        A list of the SLURM job IDs that your would like this job to depend on.
        If `None` then it is assumed there is no dependency on any other job.
        Default: `None`

    depend_type : str [optional]
        The type of slurm dependancy required. For example if you wanted the
        job to run after the jobs have been terminated use 'afterany'.
        Default: "afterok"

    submit : boolean [optional]
        Whether to write and submit the job scripts (`True`) or only write the scripts (`False`).
        Default: `True`

    outfile : str [optional]
        The output file name if "`name`.out" is not desirable.
        Default: `None` (i.e. "`batch_dir`/`name`.out")

    queue : str [optional]
        The type of queue you require (cpuq, gpuq or copyq) then the script will
        choose the correct partitions and clusters for the job to run on
        Default: "cpuq"

    export : str [optional]
        Switch that lets SLURM use your login environment on the compute
        nodes ("ALL") or not ("NONE").
        Default: "None"

    gpu_res : int [optional]
        Number of GPUs that the SLURM job will reserve.
        Default: "None"

    mem : int [optional]
        The MB of ram required for your slurm job.
        Default: 8192

    cpu_threads : int [optional]
        The number of cpu threads required for your slurm job.
        Default: 1


    Returns
    -------
    jobid : int
        The unique SLURM job ID associated with the submitted job.
    """
    if slurm_kwargs is None:
        slurm_kwargs={}

    #Work out which partition and cluster to use based on the supercomputer
    #(in config file) and queue required
    comp_config = load_config_file()
    if queue == 'cpuq':
        cluster   = comp_config['cpuq_cluster']
        partition = comp_config['cpuq_partition']
    elif queue == 'gpuq':
        cluster   = comp_config['gpuq_cluster']
        partition = comp_config['gpuq_partition']
        if gpu_res is None:
            # No gpus reserved so change it to a default of 1
            gpu_res = 1
    elif queue == 'copyq':
        cluster   = comp_config['copyq_cluster']
        partition = comp_config['copyq_partition']
    elif queue == 'zcpuq':
        # Download and checks should be done on Zeus's cpuq. This will only work
        # on Galaxy as the Ozstar workflow is different
        cluster   = comp_config['zcpuq_cluster']
        partition = comp_config['zcpuq_partition']
    else:
        logger.error("No queue found, please use cpuq, gpuq or copyq")



    header = []

    if batch_dir.endswith("/") is False:
        batch_dir += "/"

    # define file names (both the batch job file and the output file)
    jobfile = batch_dir + name + ".batch"
    if not outfile:
        outfile = batch_dir + name + ".out"

    # create the header from supplied arguments
    for k,v in slurm_kwargs.items():
        if len(k) > 1:
            k = "--" + k + "="
        else:
            k = "-" + k + " "

        header.append("#SBATCH {0}{1}".format(k, v))

    # check if there are dependencies, and if so include that in the header
    if depend is not None:
        #assumes append is a list but if not will make an educated guess of how to reformat it
        if isinstance(depend, int):
            #assume it's ben given a single job id
            header.append("#SBATCH --dependency={0}:{1}".format(depend_type, depend))
        if isinstance(depend, str):
            if ":" in depend:
                #assume it has been given an already formated string
                if depend.startswith(":"):
                    depend = depend[1:]
            #or a single jobid
            header.append("#SBATCH --dependency={0}:{1}".format(depend_type, depend))
        if isinstance(depend, list):
            depend_str = ""
            for job_id in depend:
                 depend_str += ":" + str(job_id)
            header.append("#SBATCH --dependency={0}{1}".format(depend_type, depend_str))

    # add a gpu res to header
    if gpu_res is not None:
        header.append('#SBATCH --gres=gpu:{0}'.format(gpu_res))

    # add temp SSD memory to combat I/O issues. Only availble on Ozstar
    hostname = socket.gethostname()
    if temp_mem is not None:
        header.append("#SBATCH --tmp={0}GB".format(temp_mem))

    if module_dir is None:
        module_dir = comp_config['module_dir']

    # now join the header into one string
    header = "\n".join(header)

    # construct the module loads
    if load_vcstools:
        modules = ["module load vcstools/{0}\n".format(vcstools_version)]
    else:
        modules = []
    switches = []
    for m in module_list:
        if m == "vcstools":
            # don't do anything as vcstools is loaded automatically
            continue
        if "module switch" in m:
            # if a module switch command is included rather than just a module name, then add it to a separate list
            switches.append(m)
        elif "module" in m:
            modules.append("{0}\n".format(m))
        else:
            modules.append("module load {0}\n".format(m))

    # join the module loads and switches into a single string
    switches = "\n".join(switches)
    modules = "\n".join(modules)

    # join the commands into a single string
    commands = "\n".join(commands)

    # some little hacks to make jobs work on the shanghai server
    if hostname.startswith('x86') or hostname.startswith('arm'):
        if vcstools_version == 'master':
            vcstools_version = 'cpu-master'
        if export == "NONE":
            export = "ALL"
        if shebag == "#!/bin/bash -l":
            shebag = "#!/bin/bash"

    # format the template script
    tmpl = tmpl.format(shebag=shebag, script=commands, outfile=outfile, header=header,
                       switches=switches, modules=modules,
                       cluster=cluster, partition=partition,
                       export=export, account=comp_config['group_account'][queue],
                       module_dir=module_dir,
                       threads=cpu_threads, mem=mem, nice=nice)

    # write the formatted template to the job file for submission
    with open(jobfile, "w") as fh:
        fh.write(tmpl)

    # submit the jobs
    batch_submit_line = "sbatch {0}".format(jobfile)
    jobid = None
    if submit:
        submit_cmd = subprocess.Popen(batch_submit_line, shell=True, stdout=subprocess.PIPE)
        for line in submit_cmd.stdout:
            if b"Submitted" in line:
                jobid = str(line.split(b" ")[3].decode())
        if jobid is None:
            logger.debug(batch_submit_line)
            logger.debug(submit_cmd.stdout)
            return
        else:
            return jobid
    else:
        return
