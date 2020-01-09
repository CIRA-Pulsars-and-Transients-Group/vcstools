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
from scipy.optimize import curve_fit

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
    [profile, len(profile)]: list
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

def get_stokes_from_ascii(file_loc):
    """
    Retrieves the all stokes components from an ascii file

    Parameters:
    -----------
    file_loc: string
        The location of the ascii file

    Returns:
    --------
    [I, Q, U, V, len(profile)]: list
        I: list
            Stokes I
        Q: list
            Stokes Q
        U: list
            Stokes U
        V: list
            Stokes V
        len(profile): int
            The number of bins in the profile
    """
    f = open(file_loc)
    lines = iter(f.readlines())
    next(lines) #skip first line
    I=[]
    Q=[]
    U=[]
    V=[]
    for line in lines:
        thisline=line.split()
        I.append(float(thisline[3]))
        Q.append(float(thisline[4]))
        U.append(float(thisline[5]))
        V.append(float(thisline[6]))

    return [I, Q, U, V, len(I)]

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

def fill_clipped_prof(clipped_prof, search_scope=None, nan_type=0.):
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
    length = len(clipped_prof)
    if search_scope is None:
        #Search 5% ahead for non-nans
        search_scope = round(length*0.05)
    search_scope = np.linspace(1, search_scope, search_scope, dtype=int)

    #loop over all values in clipped profile
    for i, val in enumerate(clipped_prof):
        if val==nan_type and not (i+max(search_scope)) >= length:
            #look 'search_scope' indices ahead for non-nans
            for j in sorted(search_scope, reverse=True):
                #fill in nans
                if clipped_prof[i+j]==nan_type:
                    for k in range(j):
                        clipped_prof[i+k]=nan_type
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

def find_minima_maxima(profile, ignore_threshold=0, min_comp_len=0):
    """
    Finds all minima and maxima of the input profile. Assumes that the profile has noise zero-clipped.

    profile: list
        The profile with noise zero-clipped
    ignore_threshold: float
        OPTIONAL -  Maxima with values below this number will be ignored. Default: 0
    min_comp_len: float
        OPTIONAL - Minimum length of a component to be considered real. Measured in bins. Default: 0

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

def find_maxima_err(model, noise_std, maxima):
    """
    Finds the uncertainty in a maxima by using the noise_std to find the distance to the points on the model that intersect at maxima-noise_std

    Parameters:
    -----------
    model: list
        The model to test
    noise_std: float
        The standard deviation of the noise from the profile used to make the model
    maxima: list
        The location of the maximum points in bin numbers

    Returns:
    --------
    maxima_err: list
        A list with length equal to that of the maxima list containing an uncertainty value for each maximum point
    """
    x = np.linspace(0, 1, len(model))
    maxima_err = []
    for mx in maxima:
        spline = UnivariateSpline(x, model - np.full(len(x), model[int(mx+0.5)]-noise_std), s=0)
        roots = list(spline.roots())
        roots.append(mx/len(model))
        roots = sorted(roots)
        index = roots.index(mx/len(model))
        mx_less = roots[index-1]
        mx_more = roots[index+1]
        err = (mx_more - mx_less)/2
        maxima_err.append(err*len(model))

    return maxima_err

def find_widths(profile, std=None):
    """
    Attempts to find the W_10, W_50 and equivalent width of a profile by using a spline approach.
    W10 and W50 errors are estimated by using: sigma_x = sigma_y/(dy/dx)
    Weq errors are estimated by finding the average difference in Weq when you add and subtract the std from the on-pulse profile

    Parameters:
    -----------
    profile: list
        The profile to find the widths of
    std: float
        OPTIONAL - The standard deviation of the noise. If unsupplied, will return Nones for uncertainty values. Default: None

    Returns:
    --------
    [W10, W50, Weq, W10_e, W50_e, Weq_e]: list
        W10: float
            The W10 width of the profile measured in number of bins
        W50: float
            The W50 width of the profile measured in number of bins
        Weq: float
            The equivalent width of the profile measured in number of bins
        W10_e: float
            The uncertainty in W10
        W50_e:
            The uncertainty in W50
        Weq_e:
            The uncertainty in Weq
    """
    #perform spline operations
    profile = np.array(profile)
    x = np.array(list(range(len(profile))))
    amp_y = max(profile) - min(profile)
    spline10 = UnivariateSpline(x, profile - np.full(len(x), 0.1*amp_y), s=0)
    spline50 = UnivariateSpline(x, profile - np.full(len(x), 0.5*amp_y), s=0)

    #find Weq
    _, off_pulse = sigmaClip(profile)
    on_pulse=[]
    for i, data in enumerate(off_pulse):
        if np.isnan(data):
            on_pulse.append(profile[i])
    x = np.array(list(range(len(on_pulse))))
    spline0 = UnivariateSpline(x, on_pulse, s=0)
    integral = spline0.integral(0, len(on_pulse)-1)
    Weq = integral/max(on_pulse)

    #find W10 and W50
    W10_roots = spline10.roots()
    W50_roots = spline50.roots()
    W10 = W10_roots[-1] - W10_roots[0]
    W50 = W50_roots[-1] - W50_roots[0]

    if std:
        #W10 root errors
        W10_roots_e = []
        for root in W10_roots:
            W10_roots_e.append(std/spline10.derivatives(root)[1]) #index 1 is the first order derivative
        W10_e = abs(W10_roots_e[-1]) + abs(W10_roots_e[0]) + 0.5 #add 0.5 bins to error for time resolution

        #W50 root errors
        W50_roots_e = []
        for root in W50_roots:
            W50_roots_e.append(std/spline50.derivatives(root)[1])
        W50_e = abs(W50_roots_e[-1]) + abs(W50_roots_e[0]) + 0.5

        #Weq errors
        on_pulse_less = (on_pulse - std).clip(min=0)
        on_pulse_more = (on_pulse + std).clip(min=0)

        spline0 = UnivariateSpline(x, on_pulse_less, s=0)
        integral = spline0.integral(0, len(profile)-1)
        Weq_less = integral/max(on_pulse - std)

        spline0 = UnivariateSpline(x, on_pulse_more, s=0)
        integral = spline0.integral(0, len(profile)-1)
        Weq_more = integral/max(on_pulse + std)
        Weq_e = (abs(Weq-Weq_less) + abs(Weq-Weq_more))/2

    else:
        W10_e = W50_e = Weq_e = None

    return [W10, W50, Weq, W10_e, W50_e, Weq_e]

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

    #Some noisy profiles have nans in the middle of the actual pulse... fix:
    for i, val in enumerate(clipped_prof):
        if np.isnan(val):
            clipped_prof[i]=0.

    clipped_prof = fill_clipped_prof(clipped_prof)

    #Reverse the clipped profile to attain the on-pulse
    on_pulse_prof = []
    on_pulse_bins = []
    for i, val in enumerate(clipped_prof):
        if val==0.:
            on_pulse_prof.append(profile[i])
            on_pulse_bins.append(i)
        else:
            on_pulse_prof.append(0)

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

def multi_gauss(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return y

def fit_gaussian(profile, max_N=6, chi_threshold=0, min_comp_len=0, plot_name=None):
    """
    Fits multiple gaussian components to a pulse profile and finds the best number to use for a fit.
    Will always fit at least one gaussian per profile component.
    Profile components are defined by find_components().
    Each gaussian is defined by the following: y = amp * np.exp( -((x - ctr)/wid)**2)

    Parameters:
    -----------
    profile: list
        A list containing the profile data
    max_N: int
        OPTIONAL - The maximum number of gaussain components to attempt to fit. Default: 6
    chi_threshold: float
        OPTIONAL - The script will stop trying new fits when the reduced chi-squared is within this amount to unity. Default: 0
    plot_name: string
        OPTIONAL - If not none, will make a plot of the best fit with this name. Default: None

    Returns:
    --------
    [fit, chisq, popt, pcov]: list
        fit: list
            The data containing the multi-component gaussian fit to the input profile
        chisq: float
            The reduced chi-sqaured value of the fit
        popt: list
            A list of floats where each 3 numbers describes a single gaussain and are 'ctr', 'amp' and 'wid' respectively
        pcov: numpy matrix
            The covariance matrix generated by the curve_fit function
    """
    #chi sqaured evaluation
    def chsq(observed_values, expected_values, err):
        test_statistic=0
        for observed, expected in zip(observed_values, expected_values):
            test_statistic+=((float(observed)-float(expected))/float(err))**2
        return test_statistic

    #Take noise mean and normalize the profile
    _, clipped = sigmaClip(profile)
    y = np.array(profile) - np.nanmean(np.array(clipped))
    max_y = max(y)
    len_y = len(y)
    y = np.array(y)/max_y
    noise_std = np.nanstd(np.array(clipped)/max_y)

    #Find profile components
    clipped = fill_clipped_prof(clipped, search_scope=int(len(profile)/100))
    on_pulse=[]
    for i, val in enumerate(clipped):
        if not np.isnan(val):
            on_pulse.append(0)
        else:
            on_pulse.append(y[i])
    comp_dict, comp_idx = find_components(on_pulse, min_comp_len=min_comp_len)

    #Estimate gaussian parameters based on profile components
    comp_centres = []
    comp_max = []
    comp_width = []
    for i in range(max_N//len(comp_idx.keys())+1):
        for key in comp_idx.keys():
            comp_centres.append(np.mean(comp_idx[key])/len_y)
            comp_width.append((max(comp_idx[key])-min(comp_idx[key]))/(4*len_y))
            comp_max.append(max(comp_dict[key])*0.8)
    centre_guess = iter(comp_centres)
    width_guess=iter(comp_width)
    max_guess=iter(comp_max)

    n_comps=len(comp_dict.keys())
    logger.debug("Number of profile components: {0} ({1})".format(n_comps, comp_centres[:n_comps]))

    #Fit up to max_N gaussians to the profile. Evaluate profile fit using reduced chi squared
    x=np.linspace(0, 1, len(y))
    bounds=(0, 1)
    guess = []
    fit_dict = {}
    for i in range(n_comps-1):
        guess+=[next(centre_guess), next(width_guess), next(max_guess)]
    for num in range(n_comps-1, max_N):
        fit_dict[str(num+1)]={"popt":[], "pcov":[], "fit":[], "chisq":[]}
        guess += [next(centre_guess), next(width_guess), next(max_guess)]
        popt, pcov = curve_fit(multi_gauss, x, y, bounds=bounds,  p0=guess, maxfev=100000)
        fit = multi_gauss(x, *popt)
        chisq = chsq(y, fit, noise_std)
        fit_dict[str(num+1)]["popt"] = popt
        fit_dict[str(num+1)]["pcov"] = pcov
        fit_dict[str(num+1)]["fit"] = fit
        fit_dict[str(num+1)]["chisq"] = chisq/(len(y)-1)
        logger.debug("Reduced chi squared for {0} components: {1}".format(num+1, fit_dict[str(num+1)]["chisq"]))
        if abs(1-chisq/(len(y)-1)) < chi_threshold:
            break

    #Find the best fit
    best_chi_diff=np.inf
    best_fit=None
    for key in fit_dict.keys():
        diff = abs(1 - fit_dict[key]["chisq"])
        if diff < best_chi_diff:
            best_chi_diff = diff
            best_fit = key
    logger.info("Fit {0} gaussians for a reduced chi sqaured of {1}".format(best_fit, fit_dict[best_fit]["chisq"]))
    popt = fit_dict[best_fit]["popt"]
    pcov = fit_dict[best_fit]["pcov"]
    fit = fit_dict[best_fit]["fit"]
    chisq = fit_dict[best_fit]["chisq"]

    #plot the best fit
    if plot_name:
        plt.figure(figsize=(20, 12))
        plt.plot(x, y, label="Observed")
        for j in range(0, 3*int(best_fit), 3):
            z = multi_gauss(x, *popt[:j+3])
            plt.plot(x, z, "--", label="Gaussian {}".format(int((j+3)/3)))
        plt.plot(x, fit, label="Model")
        plt.legend()
        plt.savefig(plot_name)

    return [fit, chisq, popt, pcov]

def prof_eval_stickel(profile, reg_param=5e-8, plot_name=None, ignore_threshold=None, min_comp_len=None):
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
        OPTIONAL - Minimum length of a component to be considered real. Measured in bins. If none, will use 1% of total profile lengths. Default: None
    ignore_threshold: float
        OPTIONAL -  Maxima with values below this number will be ignored. Default: 0

    Returns:
    --------
    [W10, W10_e, W50, W50_e, Weq, Weq_e, maxima, maxima_e]: list
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
        maxima_e: list
            A list of floats, each correspinding to the error of the maxima of the same index. Measured in bins
    """
    #initialize minimum component length
    if min_comp_len is None:
        min_comp_len = int(len(profile)/100)

    #regularize the profile
    profile = np.array(profile)/max(profile) #for meaningful noise mean and std prints
    on_pulse_prof, noise_mean, noise_std = regularize_prof(profile, reg_param=reg_param)
    logger.debug("Profile noise mean: {}".format(noise_mean))
    logger.debug("Profile noise STD: {}".format(noise_std))

    #find widths:
    W10, W50, Weq, W10_e, W50_e, Weq_e = find_widths(on_pulse_prof, noise_std)

    #find max, min, error
    _, maxima = find_minima_maxima(on_pulse_prof, ignore_threshold=ignore_threshold, min_comp_len=min_comp_len)
    maxima_e = find_maxima_err(on_pulse_prof, noise_std, maxima)

    logger.info("W10:                   {0} +/- {1}".format(W10, W10_e))
    logger.info("W50:                   {0} +/- {1}".format(W50, W50_e))
    logger.info("Weq:                   {0} +/- {1}".format(Weq, Weq_e))
    logger.info("Maxima:                {0}".format(maxima))
    logger.info("Maxima error:          {0}".format(maxima_e))

    #plotting
    if plot_name:
        std_prof = np.array(profile) - noise_mean
        std_prof = std_prof/np.max(std_prof)
        plt.figure(figsize=(20, 12))
        plt.title(plot_name.split("/")[-1], fontsize=22)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlim(0, len(std_prof))
        plt.xlabel("Bins", fontsize=20)
        plt.ylabel("Intensity", fontsize=20)
        for i, mx in enumerate(maxima):
            plt.axvline(x=(mx + maxima_e[i])/len(profile), ls=":", lw=2, color="gray")
            plt.axvline(x=(mx - maxima_e[i])/len(profile), ls=":", lw=2, color="gray")
        plt.plot(std_prof, label="Original Profile", color="black")
        plt.plot(on_pulse_prof, label="Regularized profile", color="red")
        plt.legend(loc="upper right", prop={'size': 16})
        plt.savefig(plot_name)

    return [W10, W10_e, W50, W50_e, Weq, Weq_e, maxima, maxima_e]

def prof_eval_gfit(profile, max_N=6, chi_threshold=0.05, ignore_threshold=0.02, plot_name=None, min_comp_len=None):
    """
    Fits multiple gaussians to a profile and subsequently finds W10, W50, Weq and maxima

    Parameters:
    -----------
    profile: list
        The pulse profile to evaluate
    chi_threshold: float
        OPTIONAL - The script will stop trying new fits when the reduced chi-squared is within this amount to unity. Default: 0.05
    plot_name: string
        OPTIONAL - If not none, will make a plot of the best fit with this name. Default: None
    ignore_threshold: float
        OPTIONAL -  Maxima with values below this number will be ignored. Default: 0.02
    min_comp_len: float
        OPTIONAL - Minimum length of a component to be considered real. Measured in bins. If none, will use 1% of total profile lengths. Default: None

    Returns:
    --------
    [W10, W10_e, W50, W50_e, Weq, Weq_e, maxima, maxima_e]: list
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
        maxima_e: list
            A list of floats, each correspinding to the error of the maxima of the same index. Measured in bins
    """
    #initialize minimum component length
    if min_comp_len is None:
        min_comp_len = int(len(profile)/100)

    #Normalize, find the std
    y = np.array(profile)/max(profile)
    noise_std, clipped = sigmaClip(y)
    y = y - np.nanmean(clipped)
    y = y/max(y)

    #fit gaussians
    fit, _, popt, _ = fit_gaussian(y, max_N=6, chi_threshold=chi_threshold, min_comp_len=min_comp_len)
    fit = np.array(fit)

    #Find widths + error
    W10, W50, Weq, W10_e, W50_e, Weq_e = find_widths(fit, noise_std)

    #find max, min, error
    for i, val in enumerate(clipped):
        if np.isnan(val):
            clipped[i]=0.
    clipped = fill_clipped_prof(clipped, search_scope=int(len(profile)/100))
    on_pulse=[]
    for i, val in enumerate(clipped):
        if val!=0:
            on_pulse.append(0)
        else:
            on_pulse.append(fit[i])
    _, maxima = find_minima_maxima(on_pulse, ignore_threshold=max(on_pulse)/100, min_comp_len=int(len(profile)/100))
    maxima_e = find_maxima_err(fit, noise_std, maxima)

    logger.info("W10:                   {0} +/- {1}".format(W10, W10_e))
    logger.info("W50:                   {0} +/- {1}".format(W50, W50_e))
    logger.info("Weq:                   {0} +/- {1}".format(Weq, Weq_e))
    logger.info("Maxima:                {0}".format(maxima))
    logger.info("Maxima error:          {0}".format(maxima_e))

    #plotting
    if plot_name:
        x = np.linspace(0, 1, len(y))
        plt.figure(figsize=(30, 18))
        for j in range(0, len(popt), 3):
            z = multi_gauss(x, *popt[j:j+3])
            plt.plot(x, z, "--", label="Gaussian Component {}".format(int((j+3)/3)))
        plt.title(plot_name.split("/")[-1].split(".")[0], fontsize=22)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlim(0, 1)
        plt.xlabel("Bins", fontsize=20)
        plt.ylabel("Intensity", fontsize=20)
        for i, mx in enumerate(maxima):
            plt.axvline(x=(mx + maxima_e[i])/len(y), ls=":", lw=2, color="gray")
            plt.axvline(x=(mx - maxima_e[i])/len(y), ls=":", lw=2, color="gray")
        plt.plot(x, y, label="Original Profile", color="black")
        plt.plot(x, fit, label="Gaussian Model", color="red")
        plt.legend(loc="upper right", prop={'size': 16})
        plt.savefig(plot_name)

    return [W10, W10_e, W50, W50_e, Weq, Weq_e, maxima, maxima_e]

if __name__ == '__main__':

    loglevels = dict(DEBUG=logging.DEBUG,\
                    INFO=logging.INFO,\
                    WARNING=logging.WARNING,\
                    ERROR=logging.ERROR)

    parser = argparse.ArgumentParser(description="""A utility file for calculating a variety of pulse profile properties""",\
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    inputs = parser.add_argument_group("Inputs")
    inputs.add_argument("--bestprof", type=str, help="The pathname of the file containing the pulse profile. Use if in .pfd.bestprof format")
    inputs.add_argument("--ascii", type=str, help="The pathname of the file containing the pulse profile. Use if in ascii text format")
    inputs.add_argument("--stickel", action="store_true", help="Use a regularization process to model the pulse profile")
    inputs.add_argument("--gaussian", action="store_true", help="Model the profile using gaussian components")
    inputs.add_argument("--ignore_threshold", type=float, default=0.02, help="Maxima with values below this fraction of the profile maximum will be ignored.")
    inputs.add_argument("--min_comp_len", type=int, default=None,\
                        help="Minimum length of a component to be considered real. Measured in bins. If none, will use 1 percent of total profile length")

    stickel_inputs = parser.add_argument_group("Stickel Inputs")
    stickel_inputs.add_argument("--reg_param", type=float, default=5e-8,\
                        help="The value of the regularization parameter used for a regularization process dscribed by Stickel 2010")

    g_inputs = parser.add_argument_group("Gaussian Inputs")
    g_inputs.add_argument("--chi_threshold", type=float, default=0.05,\
                        help="The script will stop trying new fits when the reduced chi-squared is within this amount to unity")
    g_inputs.add_argument("--max_N", type=int, default=6, help="The maximum number of gaussian components to attempt to fit")

    other_inputs = parser.add_argument_group("Other Inputs")
    other_inputs.add_argument("--plot_name", type=str, help="The name of the output plot file. If none, will not plot anything")
    other_inputs.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbostity level")
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

    if args.gaussian:
        prof_eval_gfit(profile, max_N=args.max_N, chi_threshold=args.chi_threshold, ignore_threshold=args.ignore_threshold,\
                        plot_name=args.plot_name, min_comp_len=args.min_comp_len)
    elif args.stickel:
        prof_eval_stickel(profile, reg_param=args.reg_param, plot_name=args.plot_name, ignore_threshold=args.ignore_threshold,\
                        min_comp_len=args.min_comp_len)
    else:
        logger.error("Please specify either a gaussian or regularization based fitting method")
        sys.exit(1)