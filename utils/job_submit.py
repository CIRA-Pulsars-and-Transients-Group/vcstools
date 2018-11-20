#!/usr/bin/env python

import subprocess


SLURM_TMPL = """#!/bin/bash -l

#SBATCH --export={export}
#SBATCH --output={outfile}
#SBATCH --account=mwaops
#SBATCH --clusters={cluster}
#
{header}

{switches}
module use /group/mwa/software/modulefiles
module load vcstools/{version}
{modules}

{script}
"""
# NOTE: --gid option removed after discussion in helpdesk ticket GS-9370


def submit_slurm(name, commands, tmpl=SLURM_TMPL, slurm_kwargs={}, module_list=[], vcstools_version="master",
                    batch_dir="batch/", depend=None, depend_type='afterok', submit=True, 
                    outfile=None, cluster="galaxy", export="NONE"):
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
        A template header string with format place holders: export, outfile, cluster, header and script.
        This is used to create the final string to be written to the job script.
        For this function, it is required to be SLURM compliant. 
        Default: `SLURM_TMPL`

    slurm_kwargs : dict [optional]
        A dictionary of SLURM keyword, value pairs to fill in whatever is not in the template supplied to `tmpl`.
        Default: `{}` (empty dictionary, i.e. no additional header parameters)

    module_list : list of str [optional]
        A list of module names (including versions if applicable) that will be included in the header for the batch
        scripts. e.g. ["vcstools/master", "mwa-voltage/master", "presto/master"] would append
            module load vcstools/master
            module load mwa-voltage/master
            module load presto/master
        to the header of the batch script. This can also invoke "module use ..." commands.
        NOTE: /group/mwa/software/modulefiles is used and vcstools/master is loaded by default.

    vcstools_version :  str
        The version of vcstools to load. Default: master.

    batch_dir : str [optional]
        The LOCAL directory where you want to write the batch scripts (i.e. it will write to `$PWD/batch_dir`).
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

    cluster : str [optional]
        The compute cluster that the job is to run on.
        Default: "galaxy"

    export : str [optional]
        Switch that lets SLURM use your login environment on the compute nodes ("ALL") or not ("NONE").
        Default: "NONE"

    Returns
    -------
    jobid : int
        The unique SLURM job ID associated with the submitted job.
    """

    header = []

    if batch_dir.endswith("/") is False:
        batch_dir += "/"

    # define file names (both the batch job file and the output file)
    jobfile = batch_dir + name + ".batch"
    if not outfile:
        outfile = batch_dir + name + ".out"
    
    # create the header from supplied arguments
    for k,v in slurm_kwargs.iteritems():
        if len(k) > 1:
            k = "--" + k + "="
        else:
            k = "-" + k + " "
        
        header.append("#SBATCH {0}{1}".format(k, v))

    # check if there are dependencies, and if so include that in the header
    if depend is not None:
        #assumes append is a list but if not will make an educated guess of how to reformat it
        if type(depend) is int:
            #assume it's ben given a single job id
            header.append("#SBATCH --dependency={0}:{1}".format(depend_type, depend))
        if type(depend) is str:
            if ":" in depend:
                #assume it has been given an already formated string
                if depend.startswith(":"):
                    depend = depend[1:]
            #or a single jobid 
            header.append("#SBATCH --dependency={0}:{1}".format(depend_type, depend))
        if type(depend) is list:
            depend_str = ""
            for job_id in depend:
                 depend_str += ":" + str(job_id)
            header.append("#SBATCH --dependency={0}{1}".format(depend_type, depend_str))

    # now join the header into one string
    header = "\n".join(header)

    # construct the module loads
    modules = []
    switches = []
    for m in module_list:
        if m == "vcstools":
            # don't do anything as vcstools is loaded automatically
            continue
        if "module switch" in m:
            # if a module switch command is included rather than just a module name, then add it to a separate list
            switches.append(m)
        else:
            modules.append("module load {0}\n".format(m))

    # join the module loads and switches into a single string
    switches = "\n".join(switches)
    modules = "\n".join(modules)

    # join the commands into a single string
    commands = "\n".join(commands)

    # format the template script
    tmpl = tmpl.format(script=commands, outfile=outfile, header=header, switches=switches, modules=modules,
                       version=vcstools_version, cluster=cluster, export=export)

    # write the formatted template to the job file for submission
    with open(jobfile, "w") as fh:
        fh.write(tmpl)

    # submit the jobs
    batch_submit_line = "sbatch {0}".format(jobfile)
    jobid = None
    if submit:
        submit_cmd = subprocess.Popen(batch_submit_line, shell=True, stdout=subprocess.PIPE)
        for line in submit_cmd.stdout:
            if "Submitted" in line:
                jobid = str(line.split(" ")[3])
        if jobid is None:
            print batch_submit_line
            print submit_cmd.stdout
            return
        else:
            return jobid
    else:
        return