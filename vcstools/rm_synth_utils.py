import subprocess
import logging
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random
import yaml

from vcstools import rm_synth
from vcstools import prof_utils
from vcstools.gfit import gfit
from vcstools.config import load_config_file

logger = logging.getLogger(__name__)


def remove_allfiles(directory):
    """Delete all files in a directory then the directory itself.

    Parameters
    ----------
    directory : `str`
        The directory to delete.
    """
    allfiles = glob.glob(f"{directory}/*")
    for afile in allfiles:
        os.remove(afile)
    os.rmdir(f"{directory}")


def rmfit_quad(archive, phase_min, phase_max):
    """Runs the PRSCHIVE rmfit command as a python subprocess using the -w option for a quadratic fit.

    Parameters
    ----------
    archive : `str`
        The name of the archive file to take as an input,
    phase_min : `float`
        The minimum phase to begin the fit, should be a float between 0 and 1,
    phase_max : `float`
        The maximum phase to use to fit, should be a float between 0 and 1,
    """
    comp_config = load_config_file()
    commands = [comp_config["prschive_container"]]
    commands.append("rmfit")
    commands.append("-m")
    commands.append("-10,10,20")
    commands.append("-w")
    commands.append(f"{phase_min},{phase_max}")
    commands.append("-Y")
    commands.append(f"{archive}")
    subprocess.run(commands)


def find_on_pulse_ranges(I, clip_type="regular", plot_name=None):
    """Find ranges of pulse components from a pulse profile by fitting a gaussian distribution.

    Parameters
    ----------
    I : `list`
        The pulse profile.
    clip_type : `str`
        The clipping verbosity for the Gaussian fitter. Choose between regular, noisy and verbose.
    plot_name : `str`
        The name of the ouput plot. If none, will not produce one.

    Returns
    -------
    phases : `list`
        A list of phases (from 0 to 1) corresponding to the on-pulse components.
    """
    g_fitter = gfit(I, clip_type=clip_type, plot_name=plot_name)
    g_fitter.auto_gfit()
    prof_dict = g_fitter.fit_dict
    if plot_name:
        g_fitter.plot_fit()
    phases = []
    for comp_no in prof_dict["comp_idx"].keys():
        phases.append(min(prof_dict["comp_idx"][comp_no])/len(I))
        phases.append(max(prof_dict["comp_idx"][comp_no])/len(I))

    return phases


def read_rmfit_QUVflux(QUVflux):
    """Reads the freq, I, Q and U values with their errors from a QUVflux.out file generated from rmfit.

    Parameters
    ----------
    QUVflux : `str`
        The QUVflux.out file

    Returns
    -------
    freq_hz : `list`
        Frequency values in Hz.
    I : `list`
        Stokes I values.
    I_e : `list`
        Stokes I uncertainties.
    Q : `list`
        Stokes Q values.
    Q_e : `list`
        Stokes Q uncertainties.
    U : `list`
        Stokes U values.
    U_e : `list`
        Stokes U uncertainties.
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


def write_rm_to_yaml(filename, rm_dict):
    """Writes the rotation measure dictionary to a yaml file.

    Parameters
    ----------
    filename : `str`
        The pathname of the file to write to.
    rm_dict : `dict`
        A dictionary with the RM information formatted like the output from :py:meth:`vcstools.rm_synth_utils.rm_synth_pipe`.
    """
    with open(filename, "w+") as f:
        yaml.dump(rm_dict, f, default_flow_style=False)
    f.close()


def find_best_range(I, Q, U, phase_ranges):
    """Finds the phase range with the largest ratio of lin pol to stokes I.

    Parameters
    ----------
    I : `list`
        Stokes I.
    Q : `list`
        Stokes Q.
    U : `list`
        Stokes U.
    phase_ranges : `list`
        A list where every 2 entries is a phase range from 0 to 1. ie. [0.1 0.3 0.4 0.8].

    Returns
    ---------
    best_phase_range : `list`
        The most suitable phase range.
    best_max : `float`
        The maximum linear polarisation value.
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
    """Performs RM synthesis on input data.

    Parameters
    ----------
    freq_hz : `list`
        Frequency values in Hz.
    I : `list`
        Stokes I values.
    I_e : `list`
        Stokes I uncertainties.
    Q : `list`
        Stokes Q values.
    Q_e : `list`
        Stokes Q uncertainties.
    U : `list`
        Stokes U values.
    U_e : `list`
        Stokes U uncertainties.
    phi_range : `tuple`, optional
        The range of RMs to search over. |br| Default: (-300, 300).
    phi_steps : `int`, optional
        The bumber of RM steps to search over. |br| Default: 10000.
    phase_range : tuple, optional
        The phase range of the profile used in fitting. Will be displayed on plot. |br| Default: None.
    plotname : `str`, optional
        The name of the plot. If None, will not plot: |br| Default: None.

    Returns
    -------
    rm : `float`
        The rotation measure.
    rm_e : `float`
        The uncertainty in the rotation measure.
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


def read_rmsynth_out(filename):
    """Reads the ouput file from :py:meth:`vcstools.rm_synth_utils.rm_synth_pipe`

    Parameters
    ----------
    filename : `str`
        The name of the file ouput from :py:meth:`vcstools.rm_synth_utils.rm_synth_pipe`.

    Returns
    -------
    rm_dict : `dict`
        Contains the following keys:

        ``"i"``
            There are i entries where i is the number of different phase ranges (`dict`). Contains the following keys:

            ``"rm"``
                The rotation measure of this run (`float`).
            ``"rm_e"``
                The uncertainty in rm (`float`).
            ``"phase_ranges"``
                The range of phases used for this run (`tuple`).
    """
    rm_dict={}
    f = open(filename, "r")
    lines = f.readlines()
    j=0
    for line in lines:
        if line.split()[0] == "Phase:":
            rm_dict[str(j)] = {}
            rm_dict[str(j)]["phase_range"] = (float(line.split()[1]), float(line.split()[-1]))
        if line.split()[0] == "Rotation":
            rm_dict[str(j)]["rm"] = float(line.split()[2])
            rm_dict[str(j)]["rm_e"] = float(line.split()[-1])
            j += 1
    return rm_dict


def rm_synth_pipe(kwargs):
    """Performs all the nexessary operations on an archive file to attain a rotation measure through the RM synthesis technique.

    Parameters
    ----------
    archive : `str`
        The name of the archive file to use.
    work_dir : `str`, optional
        The directory to work in. |br| Default: './'.
    pulsar : `str`, optional
        The name of the puslar (used for naming purposes). |br| Default: None.
    obsid : `int`, optional
        The observation ID (used for naming purposes). |br| Default: None.
    plot : `boolean`, optional
        If True, will ouput a plot when finished. |br| Default: False.
    write : `boolean`
        If True, will write the result to a file. |br| Default: False.
    label : `str`, optional
        A label used to identify the output files. If None, will generate a 64 bit string.
    keep_QUV : `boolean`, optional
        If True, will keep the QUVflux.out file generated from rmfit.
    force_single : `boolean`, optional
        If True, will find the phase range with the greates lin_pol ratio and fit with only this (only if kwargs_rms['phase_range'] is None). |br| Default: False.
    kwagrs_rms : `dict`
        keyword arguments for :py:meth:`vcstools.rm_synth_utils.IQU_rm_synth`.
    kwargs_gfit : `dict`
        keyword arguments for :py:meth:`vcstools.prof_utils.auto_gfit`.

    Returns
    -------
    rm_dict: dictionary
        Contains the following keys:

        ``"i"``
            There are i entries where i is the number of different phase ranges (`dict`). Contains the following keys:

            ``"rm"``
                The rotation measure of this run (`float`).
            ``"rm_e"``
                The uncertainty in rm (`float`).
            ``"phase_ranges"``
                The range of phases used for this run (`tuple`).
            ``"plotname"``
                The name of the output plot. None if no plot (`str`).
            ``"label"``
                The label used (`str`).
    filename : `str`
        The path of the file that was written to. None if not written.
    """
    # Make a randomly generated label for temp directory
    randomgen = random.getrandbits(64)
    if not kwargs["label"]:
        kwargs["label"] = randomgen
    logger.info(f"Applying label: {kwargs['label']}")

    # Move to directory and make temporary working dir
    os.chdir(kwargs["work_dir"])
    os.mkdir("{}".format(randomgen))

    ## WORKING FROM A DIFFERENT DIRECTORY ##
    os.chdir("{}".format(randomgen))
    archive_path = os.path.join("..", kwargs["archive"])

    # Read the ascii file
    ascii_archive = f"{kwargs['label']}_archive.txt"
    prof_utils.subprocess_pdv(archive_path, ascii_archive)
    I, Q, U, _, _ = prof_utils.get_stokes_from_ascii(ascii_archive)
    os.remove(ascii_archive) # Remove the ascii file, we don't need it anymore

    # Find the phase range(s) to fit:
    if not kwargs["phase_ranges"]:
        gaussian_name = None
        if kwargs["plot_gfit"]:
            gaussian_name = f"{kwargs['label']}_RMsynth_gaussian_fit.png"
        try:
            kwargs["phase_ranges"] = find_on_pulse_ranges(I, clip_type=kwargs["cliptype"], plot_name=gaussian_name)
        except prof_utils.NoFitError as e:
            os.chdir("../")
            remove_allfiles(randomgen)
            raise(prof_utils.NoFitError("No profiles could be fit. Please supply on-pulse range"))
        if kwargs["plot_gfit"]: # Move this out of the working dir
            os.rename(gaussian_name, f"../{gaussian_name}")
        if kwargs["force_single"]:
            logger.info("Forcing use of a single phase range")
            kwargs["phase_ranges"], _ = find_best_range(I, Q, U, kwargs["phase_ranges"])
    logger.info(f"Using phase ranges: {kwargs['phase_ranges']}")

    # Run rmfit with -w option
    rm_dict = {}
    for i, phase_min, phase_max in zip(range(len(kwargs["phase_ranges"])//2), kwargs["phase_ranges"][0::2], kwargs["phase_ranges"][1::2]):
        logger.info(f"Performing RM synthesis on phase range: {phase_min} - {phase_max}")
        rmfit_quad(archive_path, phase_min, phase_max)
        # Read the QUVflux.out file
        QUVflux = "QUVflux.out"
        rmfreq, rmI, rmI_e, rmQ, rmQ_e, rmU, rmU_e, = read_rmfit_QUVflux(QUVflux)

        # Make plot name if needed
        if kwargs["plot_rm"]:
            plotname = f"{kwargs['label']}_"
            if len(kwargs["phase_ranges"])>2:
                plotname += f"{i}_"
            plotname += "RMsynthesis.png"
            kwargs["plotname"] = plotname
        else:
            kwargs["plotname"] = None
        kwargs["phase_range"] = (phase_min, phase_max)
        kwargs["title"] = f"{kwargs['label']} RM Synthesis"

        # Perform RM synthesis
        rm, rm_e = IQU_rm_synth(rmfreq, rmI, rmQ, rmU, rmI_e, rmQ_e, rmU_e,
                   phase_range=kwargs["phase_range"], force_single=kwargs["force_single"], title=kwargs["title"],
                   plotname=kwargs["plotname"], phi_range=kwargs["phi_range"], phi_steps=kwargs["phi_steps"])
        rm_dict[str(i)] = {}
        rm_dict[str(i)]["rm"]           = float(rm)
        rm_dict[str(i)]["rm_e"]         = float(rm_e)
        rm_dict[str(i)]["phase_range"]  = list(kwargs["phase_range"])
        rm_dict[str(i)]["plotname"]     = kwargs["plotname"]
        rm_dict[str(i)]["label"]        = kwargs["label"]

        # Move plot out of working dir
        if kwargs["plotname"]:
            os.rename(kwargs["plotname"], f"../{kwargs['plotname']}")

        # Keep quvflux if needed
        if kwargs["keep_QUV"]:
            quvflux_name = kwargs["label"]
            if len(kwargs["phase_ranges"])>2:
                quvflux_name += f"_{i}"
            quvflux_name += "_QUVflux.out"
            os.rename("QUVflux.out", f"../{quvflux_name}")

    # Write results to file
    if kwargs["write"]:
        filename = f"{kwargs['label']}_RMsynthesis.yaml"
        write_rm_to_yaml(filename, rm_dict)
        os.rename(filename, f"../{filename}")
    else:
        filename = None

    os.chdir("../")
    ## BACK TO ORIGINAL DIRECTORY ##
    # Clean directory
    remove_allfiles(randomgen)
    return rm_dict, filename