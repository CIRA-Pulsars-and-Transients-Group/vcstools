#!/usr/bin/env python

TMPL = """#!/bin/bash -l

#SBATCH --export={export}
#SBATCH --output={outfile}
#SBATCH --account=mwaops
#SBATCH --gid=mwaops
#SBATCH --clusters={cluster}
#
{header}

{script}
"""

def submit_slurm(name, commands, slurm_kwargs={}, tmpl=TMPL, 
                    batch_dir="batch/", depend=None, submit=True, 
                    outfile=None, cluster="galaxy", export="NONE"):
    """
    Making this function to cleanly submit slurm jobs using a simple template.
    This will use the <name> to setup both the .batch and .out files into <batch_dir>
    <slurm_kwargs> should be a dictionary of keyword, value pairs of anything the template header is missing
    <commands> is the actual batch script you want to run, and is expecting a list with one entry per line
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
    
    # now join everything into one string
    header = "\n".join(header)
    commands = "\n".join(commands)
    tmpl = tmpl.format(script=commands, outfile=outfile, header=header, cluster=cluster, export=export)

    # write the formatted template to the job file for submission
    with open(jobfile,"w") as fh:
        fh.write(tmpl)

    # if there are dependencies, ensure they are taken into account on the command line
    if depend is not None:
        batch_submit_line = "sbatch --dependency=afterok:{0} {1}".format(depend, jobfile) # should this just be in the header?
    else:
        batch_submit_line = "sbatch {0}".format(jobfile)

    # submit the jobs
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
