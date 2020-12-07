#! /usr/bin/env python3

import subprocess
import argparse
import logging
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random

from vcstools import rm_synth
from vcstools import prof_utils

logger = logging.getLogger(__name__)

def rmfit_quad(archive, phase_min, phase_max):
    """
    Runs the PRSCHIVE rmfit command as a python subprocess using the -w option for a quadratic fit

    Parameters:
    -----------
    archive: string
        The name of the archive file to take as an input
    phase_min: float
        The minimum phase to begin the fit, should be a float between 0 and 1
    phase_max: float
        The maximum phase to use to fit, should be a float between 0 and 1
    """
    commands=["rmfit"]
    commands.append("-m")
    commands.append("-10,10,20")
    commands.append("-w")
    commands.append(f"{phase_min},{phase_max}")
    commands.append("-Y")
    commands.append(archive)
    subprocess.run(commands)

def find_on_pulse_ranges(I, **kwargs):
    """
    Find ranges of pulse components from a pulse profile by fitting a gaussian distribution

    Parameters:
    -----------
    I:list
        The pulse profile
    **kwargs:
        Keyword arguments for prof_utils.auto_gfit

    Returns:
    --------
    phases: list
        A list of phases (from 0 to 1) corresponding to the on-pulse components
    """
    prof_dict = prof_utils.auto_gfit(I, **kwargs)
    phases = []
    for comp_no in prof_dict["comp_idx"].keys():
        phases.append(min(prof_dict["comp_idx"][comp_no])/len(I))
        phases.append(max(prof_dict["comp_idx"][comp_no])/len(I))

    return phases

def read_rmfit_QUVflux(QUVflux):
    """
    Reads the freq, I, Q and U values with their errors from a QUVflux.out file generated from rmfit

    Parameters:
    -----------
    QUVflux: string
        The QUVflux.out file

    Returns:
    --------
    List containing:
        freq_hz: list
            Frequency values in Hz
        I: list
            Stokes I values
        I_e: list
            Stokes I uncertainties
        Q: list
            Stokes Q values
        Q_e: list
            Stokes Q uncertainties
        U: list
            Stokes U values
        U_e: list
            Stokes U uncertainties
    """
    f = np.genfromtxt(QUVflux)
    freq_hz = np.array([i[1] for i in f])*1e6
    I       = np.array([i[2] for i in f])
    I_e     = np.array([i[3] for i in f])
    Q       = np.array([i[4] for i in f])
    Q_e     = np.array([i[5] for i in f])
    U       = np.array([i[6] for i in f])
    U_e     = np.array([i[7] for i in f])

    return [freq_hz, I, I_e, Q, Q_e, U, U_e]

def write_rm_to_file(filename, rm_dict):
    """
    Writes the rotation measure and its error to a file

    Parameters:
    -----------
    filename: str
        The pathname of the file to write to
    rm_dict: dictionary
        A dictionary with the RM information formatted like the output from rm_synth_pipe()
    """
    f = open(filename, "w+")
    for i in rm_dict.keys():
        rm = rm_dict[i]["rm"]
        rm_e = rm_dict[i]["rm_e"]
        phase_range = rm_dict[i]["phase_range"]
        f.write("Mearurement {} #######################\n".format(i))
        f.write("Phase:                 {0} -   {1}\n".format(*phase_range))
        f.write("Rotation Measure:    {0} +/- {1}\n".format(rm, rm_e))
        f.write("######################################\n")
    f.close()

def find_best_range(I, Q, U, phase_ranges):
    """
    Finds the phase range with the largest ratio of lin pol to stokes I

    Parameters:
    -----------
    I: list
        Stokes I
    Q: list
        Stokes Q
    U: list
        Stokes U
    phase_ranges: list
        A list where every 2 entries is a phase range from 0 to 1. ie. [0.1 0.3 0.4 0.8]

    Returns:
    ---------
    best_phase_range: list
        The most suitable phase range
    best_max: float
        The maximum linear polarisation value
    """
    lin_pol = np.sqrt(np.array(Q)**2 + np.array(U)**2)
    I_pol_ratio = list(lin_pol / np.array(I))

    best_max = 0
    for _, phase_min, phase_max in zip(range(len(phase_ranges)//2), phase_ranges[0::2], phase_ranges[1::2]):
        int_min = int(len(I)*phase_min)
        int_max = int(len(I)*phase_max)
        print(max(I_pol_ratio[int_min:int_max]))
        if max(I_pol_ratio[int_min:int_max]) > best_max:
            best_phase_range = [phase_min, phase_max]
            best_max = max(I_pol_ratio[int_min:int_max])

    return best_phase_range, best_max

def IQU_rm_synth(freq_hz, I, Q, U, I_e, Q_e, U_e, phase_range=None, force_single=False, title=None, plotname=None, phi_range=(-300, 300), phi_steps=10000):
    """
    Performs RM synthesis on input data

    Parameters:
    -----------
    freq_hz: list
        Frequency values in Hz
    I: list
        Stokes I values
    I_e: list
        Stokes I uncertainties
    Q: list
        Stokes Q values
    Q_e: list
        Stokes Q uncertainties
    U: list
        Stokes U values
    U_e: list
        Stokes U uncertainties
    phi_range: tuple
        OPTIONAL - The range of RMs to search over. Default: (-300, 300)
    phi_steps: int
        OPTINAL - The bumber of RM steps to search over. Default: 10000
    phase_range: tuple
        OPTIONAL - The phase range of the profile used in fitting. Will be displayed on plot. Default: None
    plotname: string
        OPTIONAL - The name of the plot. If None, will not plot: Default: None

    Returns:
    --------
    rm: float
        The rotation measure
    rm_e: float
        The uncertainty in the rotation measure
    """
    p = rm_synth.PolObservation(freq_hz, (I, Q, U), IQUerr=(I_e, Q_e, U_e))
    phi_axis = np.linspace(*phi_range, phi_steps)
    p.rmsynthesis(phi_axis)
    p.rmclean(cutoff=3.)
    p.get_fdf_peak()
    p.print_rmstats()
    rm = p.cln_fdf_peak_rm
    rm_e = p.cln_fdf_peak_rm_err
    norm_factor = max(abs(p.fdf))

    if plotname:
        plt.figure(figsize=(12, 10))
        rm_string='RM: {0:7.3f}+/-{1:6.3f}'.format(p.cln_fdf_peak_rm, p.cln_fdf_peak_rm_err)
        phi_min = min(phi_range)
        plt.text(min(phi_range)+0.05*abs(phi_min), 0.95, rm_string, fontsize=12)
        if phase_range:
            phase_string = "Phase range: {0:7.3f} - {1:7.3f}".format(float(phase_range[0]), float(phase_range[1]))
            plt.text(min(phi_range)+0.05*abs(phi_min), 0.9, phase_string, fontsize=12)

        plt.plot(p.rmsf_phi,abs(p.rmsf),'k-', linewidth=0.5)
        plt.plot(p.phi,abs(p.fdf/norm_factor),'r-', linewidth=0.5)
        plt.plot(p.phi,abs(p.rm_cleaned/norm_factor),'b-', linewidth=0.5)
        plt.plot(p.phi,abs(p.rm_comps/norm_factor),'g-', linewidth=0.5)
        plt.legend(('RMSF','Dirty FDF','Clean FDF','FDF Model Components'), loc='upper right')
        plt.xlabel('RM (rad/m2)', fontsize=12)
        plt.ylabel('Amplitude (Arbitrary Units)', fontsize=12)
        plt.xlim(*phi_range)
        plt.ylim(0, 1)
        if title:
            plt.title(title, fontsize=16)
        plt.savefig(plotname, bbox_inches='tight')

    return rm, rm_e

def rm_synth_pipe(kwargs):
    """
    Performs all the nexessary operations on an archive file to attain a rotation measure through the RM synthesis technique

    Parameters:
    -----------
    archive: string
        The name of the archive file to use
    work_dir: string
        OPTIONAL - The directory to work in. Default: './'
    pulsar: string
        OPTIONAL - The name of the puslar (used for naming purposes). Default: None
    obsid: int
        OPTIONAL - The observation ID (used for naming purposes). Default: None
    plot: boolean
        OPTIONAL - If True, will ouput a plot when finished. Default: False
    write: boolean
        OPTIONAL - If True, will write the result to a file. Default: False
    label: string
        OPTIONAL - A label used to identify the output files. If None, will generate a 64 bit string
    keep_QUV: boolean
        OPTIONAL - If True, will keep the QUVflux.out file generated from rmfit
    force_single: boolean
        OPTIONAL - If True, will find the phase range with the greates lin_pol ratio and fit with only this (only if kwargs_rms['phase_range'] is None). Default: False
    kwagrs_rms: dict
        keyword arguments for IQU_rm_synth()
    kwargs_gfit: dict
        keyword arguments for prof_utils.auto_gfit()

    Returns:
    --------
    rm_dict: dictionary
        Contains the following keys:
        i: dictionary
            There are i entries where i is the number of different phase ranges
            Contains the following keys:
            rm: float
                The rotation measure of this run
            rm_e: float
                The uncertainty in rm
            phase_range: tuple
                The range of phases used for this run
            plotname: str
                The name of the output plot. None if no plot
            label: str
                The label used
    filename: str
        The path of the file that was written to. None if not written
    """
    randomgen = random.getrandbits(64)
    if not kwargs["label"]:
        kwargs["label"] = randomgen
    logger.info(f"Applying label: {kwargs['label']}")

    #Move to directory and make temporary working dir
    os.chdir(kwargs["work_dir"])
    os.mkdir("{}".format(randomgen))
    os.chdir("{}".format(randomgen))
    archive = os.path.join("..", kwargs["archive"])

    #Write the .ar file to text archive
    ascii_archive = f"{kwargs['label']}_archive.txt"
    prof_utils.subprocess_pdv(archive, ascii_archive)
    I, Q, U, _, _ = prof_utils.get_stokes_from_ascii(ascii_archive)
    os.remove(ascii_archive) #remove the archive file, we don't need it anymore

    #find the phase range(s) to fit:
    if not kwargs["phase_ranges"]:
        kwargs["phase_ranges"] = find_on_pulse_ranges(I, cliptype=kwargs["cliptype"])
        if kwargs["force_single"]:
            logger.info("Forcing use of a single phase range")
            kwargs["phase_ranges"], _ = find_best_range(I, Q, U, kwargs["phase_ranges"])
    logger.info(f"Using phase ranges: {kwargs['phase_ranges']}")

    #run rmfit with -w option
    rm_dict = {}
    for i, phase_min, phase_max in zip(range(len(kwargs["phase_ranges"])//2), kwargs["phase_ranges"][0::2], kwargs["phase_ranges"][1::2]):
        logger.info(f"Performing RM synthesis on phase range: {phase_min} - {phase_max}")
        rmfit_quad(archive, phase_min, phase_max)
        #read the QUVflux.out file
        QUVflux = "QUVflux.out"
        rmfreq, rmI, rmI_e, rmQ, rmQ_e, rmU, rmU_e, = read_rmfit_QUVflux(QUVflux)

        #make plot name if needed
        if kwargs["plot"]:
            plotname = f"{kwargs['label']}_"
            if len(kwargs["phase_ranges"])>2:
                plotname += f"{i}_"
            plotname += "RMsynthesis.png"
            kwargs["plotname"] = plotname
        else:
            plotname = None
        kwargs["phase_range"] = (phase_min, phase_max)
        kwargs["title"] = ""
        if kwargs["label"]:
            kwargs["title"] = f"{kwargs['label']} "
        kwargs["title"] += "RM Synthesis"

        #perform RM synthesis
        rm, rm_e = IQU_rm_synth(rmfreq, rmI, rmQ, rmU, rmI_e, rmQ_e, rmU_e,
                   phase_range=kwargs["phase_range"], force_single=kwargs["force_single"], title=kwargs["title"],
                   plotname=kwargs["plotname"], phi_range=kwargs["phi_range"], phi_steps=kwargs["phi_steps"])
        rm_dict[str(i)] = {}
        rm_dict[str(i)]["rm"]           = rm
        rm_dict[str(i)]["rm_e"]         = rm_e
        rm_dict[str(i)]["phase_range"]  = kwargs["phase_range"]
        rm_dict[str(i)]["plotname"]     = plotname
        rm_dict[str(i)]["label"]        = kwargs["label"]

        #move plot out of working dir
        if plotname:
            os.rename(plotname, f"../{plotname}")

        #keep quvflux if needed
        if kwargs["keep_QUV"]:
            quvflux_name = kwargs["label"]
            if len(kwargs["phase_ranges"])>2:
                quvflux_name += f"_{i}"
            quvflux_name += "_QUVflux.out"
            os.rename("QUVflux.out", f"../{quvflux_name}")

    if kwargs["write"]:
        filename = f"{kwargs['label']}_"
        filename += "RMsynthesis.txt"
        write_rm_to_file(filename, rm_dict)
        os.rename(filename, f"../{filename}")
    else:
        filename = None

    os.chdir("../")
    allfiles = glob.glob(f"{randomgen}/*")
    for afile in allfiles:
        os.remove(afile)
    os.rmdir(f"{randomgen}")

    return rm_dict, filename


def rm_synth_main(kwargs):
    rm_dict, _ = rm_synth_pipe(kwargs)
    for i in rm_dict.keys():
        logger.info("For phase range: {0} - {1}".format(*rm_dict[i]["phase_range"]))
        logger.info("RM: {0:7.3f} +/- {1:6.3f}".format(rm_dict[i]["rm"], rm_dict[i]["rm_e"]))


if __name__ == '__main__':

    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR = logging.ERROR)

    parser = argparse.ArgumentParser(description="A script that performs the RM synthesis technique",\
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group("Required Inputs:")
    required.add_argument("-f", "--archive", required=True, type=str, help="The name of the archive (.ar) file to work with")

    fitting = parser.add_argument_group("Fitting Options:")
    fitting.add_argument("--phase_ranges", type=float, nargs="+", help="The phase range(s) to fit the RM to. If unsupplied, will find the on-pulse and fit that range.\
                         Supports multiple ranges. eg. 0.1 0.15 0.55 0.62 will fit from 0.1 to 0.15 and from 0.55 or 0.62.")
    fitting.add_argument("--phi_steps", type=int, default=10000, help="The number of rm steps to use for synthesis.")
    fitting.add_argument("--phi_range", type=float, default=(-300, 300), nargs="+", help="The range of RMs so synthsize. Giving a smaller window will speed up operations.")
    fitting.add_argument("--force_single", action="store_true", help="use this tag to force using only a single phase range (if phase_ranges is unsupplied)")

    output = parser.add_argument_group("Output Options:")
    output.add_argument("--label", type=str, help="A label for the output.")
    output.add_argument("--write", action="store_true", help="Use this tag to write the results to a labelled file")
    output.add_argument("--plot", action="store_true", help="Use this tag to plot the result.")
    output.add_argument("--keep_QUV", action="store_true", help="Use this tag to keep the QUVflux.out file from rmfit.")

    gfit = parser.add_argument_group("Gaussian Fit Options")
    gfit.add_argument("--cliptype", type=str, default="regular", help="Verbosity of clipping for gaussian fitting. Options - regular, noisy, verbose")

    optional = parser.add_argument_group("additional Inputs:")
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
    kwargs = vars(args)
    rm_synth_main(kwargs)