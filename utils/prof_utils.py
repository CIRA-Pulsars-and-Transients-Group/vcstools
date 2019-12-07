#!/usr/bin/env python3

import numpy as np
import re
import sys
from scipy.interpolate import UnivariateSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from stickel import Stickel
from astropy.time import Time

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
            print("Took {0} trials to reach tolerance".format(trial+1))
            np.seterr(**old_settings)
            return oldstd, x

        if trial + 1 == ntrials:
            print("Reached number of trials without reaching tolerance level")
            np.seterr(**old_settings)
            return oldstd, x

        oldstd = newstd

def fill_clipped_prof(clipped_prof, profile):
    """
    Intended for use on noisy profiles. Fills nan values that are surrounded by non-nans to avoid discontinuities in the profile

    Parameters:
    -----------
    clipped_prof: list
        The on-pulse profile
    profile: list
        The original profile

    Returns:
    --------
    clipped_prof: list
        The clipped profile with nan gaps filled in
    """

    if len(clipped_prof) != len(profile):
        print("profiles not the same length. Exiting")
        sys.exit(1)

    #Search 5% ahead for non-nans
    length = len(profile)
    search_scope = round(length*0.05)
    search_scope = list(range(search_scope+1))[:-1]

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
        obsid: the observation ID
        pulsar: the name of the pulsar
        dm: the dispersion measure of the pulsar
        period: the period of the pulsar
        period_uncer: the uncertainty in the period measurement
        obsstart: the beginning time of the observation
        obslength: the length of the observation in seconds
        profile: a list containing the profile data
        bin_num: the number of binsn in the profile
    """

    with open(file_loc,"r") as bestprof:
        lines = bestprof.readlines()
        # Find the obsid by finding a 10 digit int in the file name
        obsid = re.findall(r'(\d{10})', lines[0])[0]
        try:
            obsid = int(obsid)
        except:
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
        for p in range(len(orig_profile)):
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
        profile: a list containing the profile data
        bin_num: the number of binsn in the profile 
    """
    with open(file_loc,"r") as f_ascii:
        lines = f_ascii.readlines()

        orig_profile = []
        for l in lines[1:]:
            orig_profile.append(float(l.split(" ")[-1]))
        bin_num = len(orig_profile)
        profile = np.zeros(bin_num)
        min_prof = min(orig_profile)
        for p in range(len(orig_profile)):
            profile[p] = orig_profile[p] - min_prof
        #maybe centre it around the pulse later
    return [profile, bin_num]

def prof_eval(profile):
    
    #clip the profile to find the on-pulse
    _, clipped_prof = sigmaClip(profile)
    #Reverse the nans and non-nans:
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
    x = np.array(list(range(len(profile))))
    y = np.array(on_pulse_prof)
    data = np.column_stack((x, y))
    stkl = Stickel(data)
    stkl.smooth_y(1e-7)
    on_pulse_prof = stkl.yhat

    #subract noise mean and normalize profile
    noise_mean = np.nanmean(clipped_prof)
    std_prof = np.array(profile) - noise_mean
    std_prof = std_prof/np.max(std_prof)
    on_pulse_prof = np.array(on_pulse_prof) - noise_mean
    on_pulse_prof = on_pulse_prof/np.max(on_pulse_prof)
    
    #set noise=0
    for i, val in enumerate(on_pulse_prof):
        if val<0:
            on_pulse_prof[i]=0
    
    #perform spline operations
    spline1 = UnivariateSpline(x, on_pulse_prof - np.full(len(x), 0.02), s=0)
    spline10 = UnivariateSpline(x, on_pulse_prof - np.full(len(x), 0.1), s=0)
    spline50 = UnivariateSpline(x, on_pulse_prof - np.full(len(x), 0.5), s=0)

    w_eq = spline1.roots()
    w10_roots = spline10.roots()
    w50_roots = spline50.roots()
    #print("Spline 1 roots: {}".format(spline1.roots()))
    print("Spline 10 roots: {}".format(spline10.roots()))
    print("Spline 50 roots: {}".format(spline50.roots()))
    
    #find maxima
    spline_maxima = UnivariateSpline(x, on_pulse_prof, s=0.0, k=4)
    maxima = spline_maxima.derivative()
    #rint("Maxima: {}".format(maxima.roots()))
    print(on_pulse_prof)

    #plt.plot(spline(on_pulse_prof))
    plt.figure(figsize=(20, 12))
    plt.plot(maxima(x))
    plt.plot(on_pulse_prof)
    plt.plot(std_prof)
    plt.savefig("on_pulse.png")

    return spline10


if __name__ == '__main__':
    #ar = get_from_bestprof("/group/mwaops/vcs/1256407632/data_products/04:59:51.94_-02:10:06.62/1256407632_1024_bins_PSR_0459-0210.pfd.bestprof")
    ar = get_from_bestprof("/group/mwaops/vcs/1257617424/data_products/05:14:52.19_-44:08:37.38/1257617424_512_bins_PSR_0514-4408.pfd.bestprof")
    profile = ar[-2]
    spline = prof_eval(profile)