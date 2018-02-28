#!/usr/bin/env python

SLURM_TMPL = """#!/bin/bash -l

#SBATCH --export={export}
#SBATCH --output={outfile}
#SBATCH --account=mwaops
#SBATCH --gid=mwaops
#SBATCH --clusters={cluster}
#
{header}

{script}
"""

def submit_slurm(name, commands, tmpl=SLURM_TMPL, slurm_kwargs={}, 
                    batch_dir="batch/", depend=None, submit=True, 
                    outfile=None, cluster="galaxy", export="NONE"):
    """
    Making this function to cleanly submit slurm jobs using a simple template.

    Parameters
    ----------
    name : str
        The base name that is used to create the "`name.batch`" and "`name.out`" files.

    commands : list of strs
        The actual bash script commands you wnat to run. 
        Expects a list where each element is a single line of the bash script.

    tmpl : str
        A template header string with format place holders: export, outfile, cluster, header and script.
        This is used to create the final string to be written to the job script `outfile`.
        For this function, it is required to be SLURM compliant. 
        Default: `SLURM_TMPL`

    slurm_kwargs : dict [optional]
        A dictionary of SLURM keyword, value pairs to fill in whatever is not in the template supplied to `tmpl`.
        Default: `{}` (empty dictionary, i.e. no additional header parameters)

    batch_dir : str [optional]
        The LOCAL directory where you want to write the batch scripts (i.e. it will write to `$PWD/batch_dir`).
        Default: "`batch/`"

    depend : int or None [optional]
        The job ID that the desired job depends on to finish with the "afterok" qualifier from SLURM.
        If `None` then it is assumed there is no dependency on any other job.
        Default: `None`

    submit : boolean [optional]
        Whether to write and submit the job scripts (`True`) or only write the scripts (`False`).
        Default: `True`

    outfile : str [optional]
        The output file name if "`name.out`" is not desirable.
        Default: `None` (i.e. "`batch_dir/name.out`")

    cluster : str [optional]
        The compute cluster that the job is to run on.
        Default: "`galaxy`"

    export : str [optional]
        Switch that lets SLURM use your login environment on the compute nodes ("`ALL`") or not ("`NONE`").
        Default: "`NONE`"

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
        header.append("#SBATCH --dependency=afterok:{0}".format(depend))
    
    # now join everything into one string
    header = "\n".join(header)
    commands = "\n".join(commands)
    tmpl = tmpl.format(script=commands, outfile=outfile, header=header, cluster=cluster, export=export)

    # write the formatted template to the job file for submission
    with open(jobfile, "w") as fh:
        fh.write(tmpl)

    # submit the jobs
    batch_submit_line = "sbatch {0}".format(jobfile)
    if submit:
        submit_cmd = subprocess.Popen(batch_submit_line, shell=True, stdout=subprocess.PIPE)
   
        # the standard output to stdout after submitting a job is:
        #   Submitted batch job [jobid] on cluster [cluster]
        # thus the 4th "word" is the job ID
        # TODO: confirm and then uncomment
#        for line in submit_cmd.stdout:
#            if "Submitted" in line:
#                jobid = str(line.split(" ")[3])

#    return jobid
