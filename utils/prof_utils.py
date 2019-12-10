#!/usr/bin/env python3

import numpy as np
import re
import sys
import logging
import argparse
from scipy.interpolate import UnivariateSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.time import Time

from stickel import Stickel

logger = logging.getLogger(__name__)

def get_from_bestprof(file_loc):
    """
    Get info from a bestprof file

    Parameters:
    -----------
    file_loc: string
        The path to the bestprof file

    Returns:
    --------
    [obsid, pulsar, dm, period, period_uncer, obsstart, obslength, profile, bin_num]: list
        obsid: int
            The observation ID
        pulsar: string
            The name of the pulsar
        dm: float
            The dispersion measure of the pulsar
        period: float
            The period of the pulsar
        period_uncer: float
            The uncertainty in the period measurement
        obsstart: int
            The beginning time of the observation
        obslength: float
            The length of the observation in seconds
        profile: list
            A list of floats containing the profile data
        bin_num: int
            The number of bins in the profile
    """

    with open(file_loc,"r") as bestprof:
        lines = bestprof.readlines()
        # Find the obsid by finding a 10 digit int in the file name
        obsid = re.findall(r'(\d{10})', lines[0])[0]
        try:
            obsid = int(obsid)
        except ValueError:
            obsid = None

        pulsar = str(lines[1].split("_")[-1][:-1])
        if not (pulsar.startswith('J') or pulsar.startswith('B')):
            pulsar = 'J{0}'.format(pulsar)

        dm = lines[14][22:-1]

        period = lines[15][22:-1]
        period, period_uncer = period.split('  +/- ')

        mjdstart = Time(float(lines[3][22:-1]), format='mjd', scale='utc')
        # Convert to gps time
        obsstart = int(mjdstart.gps)

        # Get obs length in seconds by multipling samples by time per sample
        obslength = float(lines[6][22:-1])*float(lines[5][22:-1])

        # Get the pulse profile
        orig_profile = []
        for l in lines[27:]:
            orig_profile.append(float(l.split()[-1]))
        bin_num = len(orig_profile)
        profile = np.zeros(bin_num)

        # Remove min
        min_prof = min(orig_profile)
        for p, _ in enumerate(orig_profile):
            profile[p] = orig_profile[p] - min_prof

    return [obsid, pulsar, dm, period, period_uncer, obsstart, obslength, profile, bin_num]

def get_from_ascii(file_loc):
    """
    Retrieves the profile from an ascii file

    Parameters:
    -----------
    file_loc: string
        The location of the ascii file

    Returns:
    --------
    [profile, bin_num]: list
        profile: list
            A list of floats containing the profile data
        len(profile): int
            The number of bins in the profile
    """

    f = open(file_loc)
    lines = iter(f.readlines())
    next(lines) #skip first line
    profile=[]
    for line in lines:
        thisline=line.split()
        profile.append(float(thisline[3]))

    return [profile, len(profile)]

def sigmaClip(data, alpha=3, tol=0.1, ntrials=10):
    """
    Sigma clipping operation:
    Compute the data's median, m, and its standard deviation, sigma.
    Keep only the data that falls in the range (m-alpha*sigma,m+alpha*sigma) for some value of alpha, and discard everything else.

    This operation is repeated ntrials number of times or until the tolerance level is hit.

    Parameters:
    -----------
    data: list
        A list of floats - the data to clip
    alpha: float
        OPTIONAL - Determines the number of sigmas to use to determine the upper nad lower limits. Default=3
    tol: float
        OPTIONAL - The fractional change in the standard deviation that determines when the tolerance is hit. Default=0.1
    ntrils: int
        OPTIONAL - The maximum number of times to apply the operation. Default=10

    Returns:
    --------
    oldstd: float
        The std of the clipped data
    x: list
        The data list that contains only noise, with nans in place of 'real' data
    """
    x = np.copy(data)
    oldstd = np.nanstd(x)
    #When the x[x<lolim] and x[x>hilim] commands encounter a nan it produces a
    #warning. This is expected because it is ignoring flagged data from a
    #previous trial so the warning is supressed.
    old_settings = np.seterr(all='ignore')
    for trial in range(ntrials):
        median = np.nanmedian(x)

        lolim = median - alpha * oldstd
        hilim = median + alpha * oldstd
        x[x<lolim] = np.nan
        x[x>hilim] = np.nan

        newstd = np.nanstd(x)
        tollvl = (oldstd - newstd) / newstd

        if tollvl <= tol:
            logger.debug("Took {0} trials to reach tolerance".format(trial+1))
            np.seterr(**old_settings)
            return oldstd, x

        if trial + 1 == ntrials:
            logger.warn("Reached number of trials without reaching tolerance level")
            np.seterr(**old_settings)
            return oldstd, x

        oldstd = newstd

def fill_clipped_prof(clipped_prof, profile, search_scope=None):
    """
    Intended for use on noisy profiles. Fills nan values that are surrounded by non-nans to avoid discontinuities in the profile

    Parameters:
    -----------
    clipped_prof: list
        The on-pulse profile
    profile: list
        The original profile
    search_scope: int
        The number of bins to search for non-nan values. If None, will search 5% of the total bin number. Default:None.

    Returns:
    --------
    clipped_prof: list
        The clipped profile with nan gaps filled in
    """
    #sanity check
    if len(clipped_prof) != len(profile):
        logger.error("profiles not the same length. Exiting")
        sys.exit(1)

    if search_scope is None:
        #Search 5% ahead for non-nans
        length = len(profile)
        search_scope = round(length*0.05)
    search_scope = np.linspace(1, search_scope, search_scope, dtype=int)

    #loop over all values in clipped profile
    for i, val in enumerate(clipped_prof):
        if val==0. and not (i+max(search_scope)) >= length:
            #look 'search_scope' indices ahead for non-nans
            for j in sorted(search_scope, reverse=True):
                #fill in nans
                if clipped_prof[i+j]!=0:
                    for j in search_scope:
                        clipped_prof[i+j]=profile[i+j]
                    break

    return clipped_prof

def find_components(profile, min_comp_len=5):
    """
    Given a profile in which the noise is clipped to 0, finds the components that are clumped together.

    Parameters:
    -----------
    profile: list
        A list of floats describing the profile where the noise has been clipped to zero
    min_comp_len: float
        OPTIONAL - Minimum length of a component to be considered real. Measured in bins. Default: 5

    Returns:
    --------
    component_dict: dictionary
        dict["component_x"] contains an array of the component x
    component_idx: dictionary
        dict["component_x"] contains an array of indexes of the original profile corresponding to component x
    """
    component_dict={}
    component_idx={}
    num_components=0
    for i, val in enumerate(profile):
        if val!=0.:
            if profile[i-1]==0 or i==0:
                num_components+=1
                comp_key = "component_{}".format(num_components)
                component_dict[comp_key]=[]
                component_idx[comp_key]=[]
            component_dict[comp_key].append(val)
            component_idx[comp_key].append(i)

    del_comps = []
    for comp_key in component_dict.keys():
        if len(component_dict[comp_key]) < min_comp_len:
            del_comps.append(comp_key)
    for i in del_comps:
        del component_dict[i]
        del component_idx[i]

    return component_dict, component_idx

def find_minima_maxima(profile, ignore_threshold=0.02, min_comp_len=5):
    """
    Finds all minima and maxima of the input profile. Assumes that the profile has noise zero-clipped.

    profile: list
        The profile with noise zero-clipped
    ignore_threshold: float
        OPTIONAL -  Default: 0.02
    min_comp_len: float
        OPTIONAL - Minimum length of a component to be considered real. Measured in bins. Default: 5

    Returns:
    --------
    minima: list
        A list of floats corresponding to the bin location of the profile minima
    maxima: list
        A list of floats corresponding to the bin location of the profile maxima
    """
    #If there is more than one component, find each one
    comp_dict, comp_idx = find_components(profile, min_comp_len)

    maxima=[]
    minima=[]
    #loop over each profile component
    for key in comp_dict.keys():
        x = np.linspace(0, len(comp_dict[key])-1, len(comp_dict[key]), dtype=int)
        spline = UnivariateSpline(x, comp_dict[key], s=0.0, k=4)
        comp_roots = spline.derivative().roots()
        # These are the roots, we want to split maxima and minima ^^
        comp_maxima=[]
        comp_minima=[]
        for i, root in enumerate(comp_roots):
            idx = int(root)
            left = comp_dict[key][idx-1]
            if left>comp_dict[key][idx]:
                comp_minima.append(root)
            else:
                comp_maxima.append(root)
        #Turn the root locations into locations on profile, not on component
        for root in comp_minima:
            abs_root = root + comp_idx[key][0]
            minima.append(abs_root)
        for root in comp_maxima:
            abs_root = root + comp_idx[key][0]
            maxima.append(abs_root)

    ignore_idx = []
    for i, mx in enumerate(maxima):
        if max(profile[int(mx-1):int(mx+1)]) < ignore_threshold*max(profile):
            ignore_idx.append(i)
    for i in sorted(ignore_idx, reverse=True):
        del maxima[i]

    return minima, maxima

def find_widths(profile):
    """
    Attempts to find the W_10, W_50 and equivalent width of a profile by using a spline approach

    Parameters:
    -----------
    profile: list
        The profile to find the widths of

    Return:
    W10: float
        The W10 width of the profile measured in number of bins
    W50: float
        The W50 width of the profile measured in number of bins
    Weq: float
        The equivalent width of the profile measured in number of bins
    """
    #perform spline operations
    x = np.array(list(range(len(profile))))
    spline0 = UnivariateSpline(x, profile, s=0)
    spline10 = UnivariateSpline(x, profile - np.full(len(x), 0.1), s=0)
    spline50 = UnivariateSpline(x, profile - np.full(len(x), 0.5), s=0)

    #find Weq
    integral = spline0.integral(0, len(profile)-1)
    Weq = integral/max(profile)

    #find W10 and W50
    W10_roots = spline10.roots()
    W50_roots = spline50.roots()
    W10=0
    W50=0
    for i in range(0, len(W10_roots), 2):
        W10 += W10_roots[i+1] - W10_roots[i]
    for i in range(0, len(W50_roots), 2):
        W50 += W50_roots[i+1] - W50_roots[i]

    return W10, W50, Weq

def regularize_prof(profile, reg_param=5e-8):
    """
    Applies the following operations to a pulse profile:
    1) Sigma Clip
    2) Find the on-pulse bins from Sigma Clip
    3) Fill unnecessary nans from the clipping process
    4) Apply Stickel 2010 regularization
    5) Normalize the profile
    6) Set the off-pulse to zero

    Parameters:
    -----------
    profile: list
        The pulse profile
    reg_param: float
        OPTIONAL - the regularization parameter used for smoothing the profile. Default:5e-8

    Returns:
    --------
    on_pulse_prof: list
        The pulse profile after the above operations
    noise_mean: float
        The mean of the profile noise
    """
    #clip the profile to find the on-pulse
    noise_std, clipped_prof = sigmaClip(profile)
    #Reverse the nan-profile. Then replace nans with zero:
    on_pulse_prof = []
    on_pulse_bins = []
    for i, val in enumerate(clipped_prof):
        if np.isnan(val):
            on_pulse_prof.append(profile[i])
            on_pulse_bins.append(i)
        else:
            on_pulse_prof.append(0)

    #Some noisy profiles have nans in the middle of the actual pulse... fix:
    on_pulse_prof = fill_clipped_prof(on_pulse_prof, profile)

    #regularization with Stickel - smoothing
    x = np.linspace(0,len(profile)-1, len(profile), dtype=int)
    y = np.array(on_pulse_prof)
    data = np.column_stack((x, y))
    stkl = Stickel(data)
    stkl.smooth_y(reg_param)
    on_pulse_prof = stkl.yhat

    #subract noise mean and normalize profile
    noise_mean = np.nanmean(clipped_prof)
    on_pulse_prof = np.array(on_pulse_prof) - noise_mean
    on_pulse_prof = on_pulse_prof/np.max(on_pulse_prof)
    #reset noise to 0
    for i, val in enumerate(on_pulse_prof):
        if val<0:
            on_pulse_prof[i]=0

    return on_pulse_prof, noise_mean, noise_std

def prof_eval(profile, reg_param=5e-8, plot_name=None, ignore_threshold=None, min_comp_len=5):
    """
    Transforms a profile by removing noise and applying a regularization process and subsequently finds W10, W50, Weq and maxima

    Parameters:
    -----------
    profile: list
        The pulse profile to evaluate
    reg_param: float
        OPTIONAL - the regularization parameter used for smoothing the profile. Default:5e-8
    plot_name: string
        OPTIONAL - The name of the plot to output. If None, will not plot anything. Default: None
    min_comp_len: float
        OPTIONAL - Minimum length of a component to be considered real. Measured in bins. Default: 5

    Returns:
    --------
    [W10, W10_e, W50, W50_e, Weq, Weq_e, maxima, minima]: list
        W10: float
            The W10 width of the profile measured in number of bins
        W10_e: float
            The uncertainty in the W10
        W50: float
            The W50 width of the profile measured in number of bins
        W50_e: float
            The uncertainty in the W50
        Weq: float
            The equivalent width of the profile measured in number of bins
        Weq_e: float
            The uncertainty in the equivalent width
        maxima: list
            A list of floats corresponding to the bin location of any profile maxima
        minima: list
            A list of floats correspodning to the bin location of any profile minima
    """

    #regularize the profile
    profile = np.array(profile)/max(profile) #for meaningful noise mean and std prints
    on_pulse_prof, noise_mean, noise_std = regularize_prof(profile, reg_param=reg_param)
    logger.debug("Profile noise mean: {}".format(noise_mean))
    logger.debug("Profile noise STD: {}".format(noise_std))

    #find widths:
    W10, W50, Weq = find_widths(on_pulse_prof)
    #Width uncertainties: Find the average difference between widths the profiles with noise std added and subtracted. Then add 0.5 bins
    less_std_prof, _, _  = regularize_prof(np.array(profile) - noise_std, reg_param=reg_param)
    W10_less, W50_less, Weq_less = find_widths(less_std_prof)
    more_std_prof, _, _ = regularize_prof(np.array(profile) + noise_std, reg_param=reg_param)
    W10_more, W50_more, Weq_more = find_widths(more_std_prof)
    W10_e = (abs(W10-W10_less) + abs(W10-W10_more))/2 + 0.5
    W50_e = (abs(W50-W50_less) + abs(W50-W50_more))/2 + 0.5
    Weq_e = (abs(Weq-Weq_less) + abs(Weq-Weq_more))/2 + 0.5

    #find max, min. Root uncertainties are 0.5 bins
    minima, maxima = find_minima_maxima(on_pulse_prof, ignore_threshold=ignore_threshold, min_comp_len=min_comp_len)

    logger.info("W10:                   {0} +/- {1}".format(W10, W10_e))
    logger.info("W50:                   {0} +/- {1}".format(W50, W50_e))
    logger.info("Weq:                   {0} +/- {1}".format(Weq, Weq_e))
    logger.info("Maxima:                {0}".format(maxima))
    logger.info("Minima:                {0}".format(minima))

    #plotting
    if plot_name:
        std_prof = np.array(profile) - noise_mean
        std_prof = std_prof/np.max(std_prof)
        plt.figure(figsize=(20, 12))
        plt.plot(on_pulse_prof, label="Regularized profile")
        plt.plot(std_prof, label="Original Profile")
        plt.title(plot_name.split("/")[-1], fontsize=22)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlim(0, len(std_prof))
        plt.xlabel("Bins", fontsize=20)
        plt.ylabel("Intensity", fontsize=20)
        for mx in maxima:
            plt.axvline(x=mx, ls=":", lw=2, color="gray")
        plt.legend(loc="upper right", prop={'size': 16})
        plt.savefig(plot_name)

    return [W10, W10_e, W50, W50_e, Weq, Weq_e, maxima, minima]

if __name__ == '__main__':

    loglevels = dict(DEBUG=logging.DEBUG,\
                    INFO=logging.INFO,\
                    WARNING=logging.WARNING,\
                    ERROR=logging.ERROR)

    parser = argparse.ArgumentParser(description="""A utility file for calculating a variety of pulss profile properties""",\
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bestprof", type=str, help="The pathname of the file containing the pulse profile. Use if in .pfd.bestprof format")
    parser.add_argument("--ascii", type=str, help="The pathname of the file containing the pulse profile. Use if in ascii text format")
    parser.add_argument("--reg_param", type=float, default=5e-8,\
                        help="The value of the regularization parameter used for a regularization process dscribed by Stickel 2010")
    parser.add_argument("--plot_name", type=str, help="The name of the output plot file. If empty, will not plot anything")
    parser.add_argument("--ignore_threshold", type=float, default=0.02, help="Maxima with values below this fraction of the profile maximum will be ignored.")
    parser.add_argument("--min_comp_len", type=int, default=5, help="The minimum length in bins that a profile component must be to be considered real")
    parser.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbostity level")
    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if args.bestprof:
        profile = get_from_bestprof(args.bestprof)[-2]
    elif args.ascii:
        profile = get_from_ascii(args.ascii)[0]
    else:
        logger.error("Please supply either an ascii or bestprof profile")
        sys.exit(1)

    a = prof_eval(profile, reg_param=args.reg_param, plot_name="on_pulse.png", ignore_threshold=args.ignore_threshold,\
                                    min_comp_len=args.min_comp_len)