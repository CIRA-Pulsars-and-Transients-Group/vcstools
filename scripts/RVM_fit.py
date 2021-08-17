#!/usr/bin/env python3

from scipy.optimize import curve_fit
import numpy
import argparse
import logging
import os
import numpy as np
import json

from vcstools import prof_utils
from vcstools.general_utils import setup_logger

logger = logging.getLogger(__name__)

def get_phase_ranges(archive: str):
    p, _ = prof_utils.get_from_ascii(archive)
    try:
        pdict = prof_utils.auto_gfit(p)
        phase_ranges=[]
        for component in pdict["comp_idx"].keys():
            phase_ranges.append(pdict["comp_idx"][component][0]/len(pdict["profile"]))
            phase_ranges.append(pdict["comp_idx"][component][-1]/len(pdict["profile"]))
    except prof_utils.NoFitError as e:
        logger.warn("No suitable fit found for profile. Will use entire phase range.")
        phase_ranges = [0, 1]

    return phase_ranges

def check_input_kwargs(kwargs):
    if not os.path.isfile(kwargs["archive"]):
        raise AssertionError(f"File {kwargs['archive']} does not exist or is not a file")
    if kwargs["phase_ranges"]:
        if len(kwargs["phase_ranges"])%2 != 0:
            raise AssertionError(f"Unever number of phases given. Cannot interpret phase range: {kwargs['phase_ranges']}")
        if any(i > 1 for i in kwargs["phase_ranges"]) or any(i < 0 for i in kwargs['phase_ranges']):
            raise AssertionError(f"Phase ranges should be between 0 and 1: {kwargs['phase_ranges']}")
    if kwargs["outfile"]:
        dirname = os.path.dirname(kwargs["outfile"])
        if not os.path.exists(dirname) and not dirname == "":
            raise AssertionError(f"Directory to store outfile does not exist: {kwargs['outfile']}")

    return kwargs

def process_args(kwargs):
    #get PA, error and X positions, removing phase values
    prof_utils.subprocess_pdv(kwargs["archive"], outfile="RVM_fit_archive.txt", pdvops="-FTtlZ")
    ar_read = np.genfromtxt("RVM_fit_archive.txt", skip_header=1)
    if not kwargs["phase_ranges"]:
        kwargs["phase_ranges"] = get_phase_ranges("RVM_fit_archive.txt")
    logger.info(f"Using ranges: {kwargs['phase_ranges']}")
    PA_ = np.array([i[8] for i in ar_read])
    PA_e = np.array([i[9] for i in ar_read])
    X_PA = np.linspace(0, 1, len(PA_))
    #remove values lying out of the phases
    ranges = []
    for i, j in zip(kwargs["phase_ranges"][0::2], kwargs["phase_ranges"][1::2]):
        ranges.append((i, j))
    for i, val in enumerate(X_PA):
        if not any(lower <= val <= upper for (lower, upper) in ranges):
            X_PA[i] = 0 #if the value isn't within range, set to 0
    #a[b != 0] clips all values in a for indices where b == 0
    # clip zero values from PA
    PA = PA_[PA_ != 0]
    PA_e = PA_e[PA_ != 0]
    X_PA = X_PA[PA_ != 0]
    #clip zero values from phase ranges
    PA = PA[X_PA != 0]
    PA_e = PA_e[X_PA != 0]
    X_PA = X_PA[X_PA != 0]

    fit_dict = {"PA":PA, "PA_e":PA_e, "X_PA":X_PA, "outfile":kwargs["outfile"]}
    return fit_dict

def check_processed_kwargs(kwargs):
    print(kwargs["X_PA"])
    print(kwargs["PA"])
    if  len(kwargs["X_PA"])==0 or len(kwargs["PA"])==0:
        raise AssertionError("No nonzero values in phase range")

def _fit(kwargs):
    PA = kwargs["PA"]
    PA_e = kwargs["PA_e"]
    X_PA = kwargs["X_PA"]
    outfile = kwargs["outfile"]
    def analytic_pa(phi, alpha, beta, psi_0, phi_0):
        #Inputs should be in radians
        numerator = np.sin(alpha) * np.sin(phi - phi_0)
        denominator = np.sin(beta + alpha) * np.cos(alpha) - np.cos(beta + alpha) * np.sin(alpha) * np.cos(phi - phi_0)
        return np.arctan2(numerator,denominator) + psi_0

    min_rad_45 = np.deg2rad(-45) #alpha/beta
    max_rad_45 = np.deg2rad(45)
    min_rad_90 = np.deg2rad(-90) #psi_0
    max_rad_90 = np.deg2rad(90)
    min_rad_180 = np.deg2rad(-180) #phi_0
    max_rad_180 = np.deg2rad(180)

    initial = (0,0,0,0)
    lower_bounds = (min_rad_45, min_rad_45, min_rad_90, min_rad_180)
    upper_bounds = (max_rad_45, max_rad_45, max_rad_90, max_rad_180)
    bounds = [lower_bounds, upper_bounds]
    popt, pcov = curve_fit(analytic_pa, X_PA, PA, p0=initial, bounds=bounds, sigma=PA_e)
    logger.info(f"Alpha:    {np.rad2deg(popt[0])}")
    logger.info(f"Beta:     {np.rad2deg(popt[1])}")
    logger.info(f"Psi_0:    {np.rad2deg(popt[2])}")
    logger.info(f"Phi_0:    {np.rad2deg(popt[3])}")
    if outfile:
        outdict = {"alpha":popt[0], "beta":popt[1], "psi_0":popt[2], "phi_0":popt[3], "covariance_matrix":np.array(pcov).reshape(4,4).tolist()}
        with open(outfile, "w") as f:
            json.dump(outdict, f)
    return popt, pcov

def fit(inkwargs):
    inkwargs = check_input_kwargs(inkwargs)
    prokwargs = process_args(inkwargs)
    check_processed_kwargs(prokwargs)
    popt, pcov = _fit(prokwargs)
    return popt, pcov

if __name__ == '__main__':

    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR = logging.ERROR)

    parser = argparse.ArgumentParser(description="A script that performs a least squares fit to the Rotating Vector Model",\
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--archive", required=True, type=str, help="The pathname of the archive (.ar) file to attempt to fit")
    parser.add_argument("--outfile", type=str, help="The pathname of the file to write the ouput to. If None, will not write. Note: file outputs are in radians")
    parser.add_argument("--phase_ranges", type=float, nargs="+", help="The phase range(s) to fit the RM to. If unsupplied, will find the on-pulse and fit that range.\
                         Supports multiple ranges. eg. 0.1 0.15 0.55 0.62 will fit values from 0.1 to 0.15 and from 0.55 to 0.62.")
    parser.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbostity level")
    args = parser.parse_args()

    # set up the logger for stand-alone execution
    logger = setup_logger(logger, log_level=loglevels[args.loglvl])

    fit(vars(args))