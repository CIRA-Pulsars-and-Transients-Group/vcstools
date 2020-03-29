#! /usr/bin/env python3

import subprocess
import argparse
import logging
import os
import numpy as np

import prof_utils

logger = logging.getLogger(__name__)

def subprocess_pdv(archive, outfile="archive.txt", pdvops="-FTt"):

    myoutput = open(outfile,'w+')
    commands=["pdv"]
    commands.append(pdvops)
    commands.append(achive)
    subprocess.run(commands, stdout=myoutput)
    myoutput.close()

def rmfit_quad(archive, phase_min, phase_max):

    commands=["rmfit"]
    commands.append("-m")
    commands.append("-10,10,20")
    commands.append("-w")
    commands.append("{0},{1}".format(phase_min, phase_max))
    commands.append("-Y")
    commands.append(archive)
    subprocess.run(commands)

def find_max_L_ratio(I, L):

    #gaussian fit to find on_pulse alpha
    prof_dict = prof_utils.auto_gfit(I)
    _, only_noise = prof_utils.sigmaClip(I, alpha=prof_dict["alpha"])
    for i, val in enumerate(only_noise):
        if np.isnan(val):
            only_noise[i] = 0
    only_noise = prof_utils.fill_clipped_prof(only_noise, nan_type=0)
    only_pulse = []
    for i, val in enumerate(only_noise):
        if val == 0:
            only_pulse.append(I[i])
        else:
            only_pulse.append(0)

    _, component_idx = prof_utils.find_components(only_pulse)
    max_L_ratio = 0.
    max_I = max(I)
    for comp_no in component_idx.keys():
        for i in component_idx[comp_no]:
            if I[i] >= 0.1*max_I and L[i] > max_L_ratio:
                max_L_ratio = L[i]
                max_L_idx = i
    return max_L_ratio, max_L_idx

def rm_synth_pipe(archive, work_dir="./", pulsar=None, obsid=None):

    #Write the .ar file to text archive
    ascii_archive = os.path.join(work_dir, "archive.txt")
    subprocess_pdv(archive, ascii_archive)
    I, Q, U, V, _ = prof_utils.get_stokes_from_ascii(ascii_archive)
    os.remove(ascii_archive) #remove the archive file, we don't need it anymore

    #find the location of maximum linear pol ratio:
    L = list(np.sqrt(np.array(Q)**2 + np.array(U)**2))
    max_L_ratio, max_L_idx = find_max_L_ratio(I, L)
    max_L_phase = max_L_idx/len(I)
    logger.info("Maximum linear polarisation ratio: {0} occurs at phase {1}".format(max_L_ratio, max_L_phase))
    phase_min = max_L_phase - 2/len(I)
    phase_max = max_L_phase + 2/len(I)

    #run rmfit with -w option
    rmfit_quad(archive, phase_min, phase_max)
    os.remove(os.path.join(work_dir, "rmfit_results.out"))
    os.remove(os.path.join(work_dir, "all_bestRMs.out"))

    #format the QUVflux.out file



if __name__ == '__main__':
    #dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR = logging.ERROR)

    #Arguments
    parser = argparse.ArgumentParser(description="A script that performs the RM synthesis technique",\
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group("Required Inputs:")
    required.add_argument("-f", "--archive", type=str, help="The name of the archive (.ar) file to work with")

    optional = parser.add_argument_group("Optional Inputs:")
    optional.add_argument("-p", "--pulsar", type=str, help="The name of the pulsar to be fit")
    optional.add_argument("-o", "--obsid", type=int, help="The name of the MWA obsid")
    optional.add_argument("-d", "--work_dir", type=str, default="./", help="The directory to work in")
    optional.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level. Default: INFO", choices=loglevels.keys())

    args = parser.parse_args()
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.propagate = False

    rm_synth_pipe(args.archive, work_dir=args.work_dir)