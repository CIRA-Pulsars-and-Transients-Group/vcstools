import numpy as np
import re
import subprocess
import logging
from scipy.interpolate import UnivariateSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.time import Time
from scipy.optimize import curve_fit

from vcstools.config import load_config_file

logger = logging.getLogger(__name__)

#---------------------------------------------------------------
class LittleClipError(Exception):
    """Raise when not enough data is clipped"""
    pass
class LargeClipError(Exception):
    """Raise when too much data is clipped"""
    pass
class NoComponentsError(Exception):
    """Raise when there are no feasible profile components"""
    pass
class ProfileLengthError(Exception):
    """Raise when a profile's legnth is too small to be fit"""
    pass
class NoFitError(Exception):
    """Raise when no gaussian fits have been found"""
    pass

#---------------------------------------------------------------
def subprocess_pdv(archive, outfile="archive.txt", pdvops="-FTt"):
    """
    Runs the pdv commnand from PSRCHIVE as a python subprocess.
    NOTE: Requires singularity loaded in environment (module load singularity)

    Parameters
    ----------
    archive: string
        The name of the archive file to run the command on
    outfile: string
        OPTIONAL - The name of the text file to write the output to. Default: outfile
    pdvops: string
        OPTIONAL - Additional options for the pdv. Default: -FTt
    """
    comp_config = load_config_file()
    with open(outfile,'w+') as f:
        commands = [comp_config["prschive_container"]]
        commands.append("pdv")
        commands.append(pdvops)
        commands.append(archive)
        a = subprocess.check_output(commands)
        f.write(a.decode("utf-8"))

#---------------------------------------------------------------
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

#---------------------------------------------------------------
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
    f.close()
    profile=[]
    for line in lines:
        thisline=line.split()
        profile.append(float(thisline[3]))

    return [profile, len(profile)]

#---------------------------------------------------------------
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
    f.close()
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

#---------------------------------------------------------------
def sigmaClip(data, alpha=3., tol=0.1, ntrials=10):
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
        OPTIONAL - Determines the number of sigmas to use to determine the upper and lower limits. Default=3
    tol: float
        OPTIONAL - The fractional change in the standard deviation that determines when the tolerance is hit. Default=0.1
    ntrials: int
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
            logger.info("Reached number of trials without reaching tolerance level")
            np.seterr(**old_settings)
            return oldstd, x

        oldstd = newstd

#---------------------------------------------------------------
def check_clip(clipped_prof, toomuch=0.8, toolittle_frac=0., toolittle_absolute=4):
    """
    Determines whether a clipped profile from sigmaClip() has been appropriately clipped by checking the number of nans.
    Raises a LittleClipError or a LargeClipError if too little or toomuch of the data has been clipped respectively.

    Parameters:
    -----------
    clipped_prof: list
        The clipped profile from sigmaClip()
    toomuch: float
        OPTIONAL - The fraction of the clipped profile beyond which is considered overclipped. Default: 0.8
    toolittle: float
        OPTIONAL - The fraction of the clipped profile below which is considered underclipped. Default: 0.
    toolittle_absolute: int
        OPTIONAL - If a profile has this many or less on-pulse bins, it is deemed not sufficient. Default: 4
    """
    num_nans = 0
    for i in clipped_prof:
        if np.isnan(i):
            num_nans += 1
    if num_nans <= toolittle_frac*len(clipped_prof) or num_nans <= toolittle_absolute:
        raise LittleClipError("Not enough data has been clipped. Condsier trying a smaller alpha value when clipping.")
    elif num_nans >= toomuch*len(clipped_prof):
        raise LargeClipError("A large portion of the data has been clipped. Condsier trying a larger alpha value when clipping.")

#---------------------------------------------------------------
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

#---------------------------------------------------------------
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
        if len(component_dict[comp_key]) < min_comp_len or max(component_dict[comp_key]) < 0.:
            del_comps.append(comp_key)
    for i in del_comps:
        del component_dict[i]
        del component_idx[i]

    if len(component_dict.keys()) == 0:
        raise NoComponentsError("No profile components have been found")

    return component_dict, component_idx

#---------------------------------------------------------------
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

def find_minima_maxima_gauss(popt, pcov, x_length, min_comp_len=0):
    """
    Finds all roots of a gaussian function

    Parameters:
    -----------
    popt: list
        A list of length 3N where N is the number of gaussians. This list contains the parameters amp, mean, centre respectively
    pcov: np.matrix
        The covariance matric corresponding to the parameters from popt
    x_length: int
        The length of the list used to fit the gaussian

    Returns:
    --------
    minima: list
        A list of the minimum points of the fit
    maxima: list
        A list of the maximum points of the fit
    """
    #Create the derivative list and spline it to find roots
    x = np.linspace(0, x_length-1, x_length)
    dy = multi_gauss_ddx(x, *popt)
    spline_dy = UnivariateSpline(x, dy, s=0)
    roots = spline_dy.roots()

    #Find which are max and min
    maxima = []
    minima = []
    for root in roots:
        idx = int(root + 0.5)
        if dy[idx-1] > dy[idx]:
            maxima.append(root)
        else:
            minima.append(root)

    minima_e = find_x_err(minima, popt, pcov)
    maxima_e = find_x_err(maxima, popt, pcov)

    return minima, maxima, minima_e, maxima_e

#---------------------------------------------------------------
def find_x_err(x, popt, pcov):
    """
    Finds the error in the horizontal position of a gaussian fit at the point x.
    Uses the equation sigma_x = sigma_y/d2ydx2 where:
    sigma_x = error in x
    d2ydx2 = second derivative of the gaussian function at point x
    sigma_y = J*C*J_T
    J = Jacobian evalutated at point x
    C = covariance matrix of gaussian fit
    J_T = transposed jacobian

    Parameters:
    -----------
    x: list
        A list of points to evaluate the error at
    popt: list
        The parameters used to describe the gaussian fit
    pcov: numpy.matrix
        The covariance matrix corresponding to popt

    Returns:
    --------
    x_err: list
        The error evaluated at each point, x
    """
    x_err = []
    for _, point in enumerate(x):
        J = jacobian_slope(point, *popt)
        d2dx2 = multi_gauss_d2dx2(point, *popt)
        JC = np.matmul(J, pcov)
        sigma_y = np.sqrt( np.matmul(JC, np.transpose(J)).item(0) )
        x_err.append(sigma_y / abs(d2dx2))
    return x_err

#---------------------------------------------------------------
def find_widths(profile, popt, pcov, alpha=3):
    """
    Attempts to find the W_10, W_50 and equivalent width of a profile by using a spline approach.
    W10 and W50 errors are estimated by using: sigma_x = sigma_y/(dy/dx)
    Weq errors are estimated by finding the average difference in Weq when you add and subtract the std from the on-pulse profile

    Parameters:
    -----------
    profile: list
        A list of floats that make up the profile
    popt: list
        The parameters that are used to create the multi-gaussian fit
    pcov: np.matrix
        The covariance matrix corresponding to the parameters from popt
    alpha: float
        OPTIONAL - The alpha value to be used in sigmaClip(). Default: 3

    Returns:
    --------
    [W10, W50, Weq, Wscat, W10_e, W50_e, Weq_e, Wscat_e]: list
        W10: float
            The W10 width of the profile measured in number of bins
        W50: float
            The W50 width of the profile measured in number of bins
        Weq: float
            The equivalent width of the profile measured in number of bins
        Wscat: float
            The scattering width of the profile measured in number of bins
        W10_e: float
            The uncertainty in W10
        W50_e: float
            The uncertainty in W50
        Weq_e: float
            The uncertainty in Weq
        Wscar_e: float
            The unceratinty in Wscat
    """
    def error_in_x_pos(pcov, popt, x):
        J = jacobian_slope(x, *popt)
        JC = np.matmul(J, pcov)
        sigma_y = np.sqrt(np.matmul(JC, np.transpose(J)).item(0))
        ddx = multi_gauss_ddx(x, *popt)
        return sigma_y/ddx

    #perform spline operations on the fit
    x = np.linspace(0, len(profile)-1, len(profile))
    fit = multi_gauss(x, *popt)
    amp_fit = max(fit) - min(fit)
    spline10 = UnivariateSpline(x, fit - np.full(len(x), 0.1*amp_fit), s=0)
    spline50 = UnivariateSpline(x, fit - np.full(len(x), 0.5*amp_fit), s=0)
    spline_s = UnivariateSpline(x, fit - np.full(len(x), 1/np.exp(1)*amp_fit), s=0)

    #find Weq using the real profile
    std, off_pulse = sigmaClip(profile, alpha=alpha)
    check_clip(off_pulse)
    on_pulse=[]
    for i, data in enumerate(off_pulse):
        if np.isnan(data):
            on_pulse.append(profile[i])
    x = np.linspace(0, len(on_pulse)-1, len(on_pulse))
    spline0 = UnivariateSpline(x, on_pulse, s=0)
    integral = spline0.integral(0, len(on_pulse)-1)
    Weq = integral/max(on_pulse)

    #find W10, W50 and Wscat
    W10_roots = spline10.roots()
    W50_roots = spline50.roots()
    Wscat_roots = spline_s.roots()
    W10 = W10_roots[-1] - W10_roots[0]
    W50 = W50_roots[-1] - W50_roots[0]
    Wscat = Wscat_roots[-1] - Wscat_roots[0]

    #W10 root errors
    err_10_1 = error_in_x_pos(pcov, popt, W10_roots[0])
    err_10_2 = error_in_x_pos(pcov, popt, W10_roots[-1])
    W10_e = np.sqrt(err_10_1**2 + err_10_2**2)

    #W50 root errors
    err_50_1 = error_in_x_pos(pcov, popt, W50_roots[0])
    err_50_2 = error_in_x_pos(pcov, popt, W50_roots[-1])
    W50_e = np.sqrt(err_50_1**2 + err_50_2**2)

    #Wscat root errors
    err_scat_1 = error_in_x_pos(pcov, popt, Wscat_roots[0])
    err_scat_2 = error_in_x_pos(pcov, popt, Wscat_roots[-1])
    Wscat_e = np.sqrt(err_scat_1**2 + err_scat_2**2)

    #Weq errors - using covariance formula
    on_pulse_less = (on_pulse - std).clip(min=0)
    spline0 = UnivariateSpline(x, on_pulse_less, s=0)
    integral = spline0.integral(0, len(profile)-1)
    dwdint = 1/max(on_pulse)**2
    dwdmax = -integral/max(on_pulse)**2
    int_e = abs(integral/max(on_pulse - std) - integral/max(on_pulse))
    max_e = std
    Weq_e = np.sqrt( dwdint**2 * int_e**2 + dwdmax**2 * max_e**2 + 2*dwdint*dwdmax*int_e*max_e )

    return [W10, W50, Weq, Wscat, W10_e, W50_e, Weq_e, Wscat_e]

#---------------------------------------------------------------
def est_sn_from_prof(prof_data, period, alpha=3.):
    """
    Estimates the signal to noise ratio from a pulse profile
    Based on code oringally writted by Nick Swainston.

    Parameters:
    -----------
    prof_data: string
        A list of floats that contains the pulse profile
    period: float
        The pulsar's period in ms
    alpha: float
        OPTIONAL - The alpha value to be used in sigmaClip(). Default: 3

    Returns:
    --------
    [sn, sn_e, scattered]
    sn: float
        The estimated signal to noise ratio
    u_sn: float
        The uncertainty in sn
    scattered: boolean
        When true, the profile is highly scattered
    """
    # Check profile is normalised
    prof_data = prof_data / max(prof_data)

    #centre the profile around the max
    shift = -int(np.argmax(prof_data))+int(len(prof_data))//2
    prof_data = np.roll(prof_data, shift)

    #find std and check if profile is scattered
    sigma, flags = sigmaClip(prof_data, tol=0.01, ntrials=100, alpha=alpha)
    check_clip(flags)
    bot_prof_min = (max(prof_data) - min(prof_data)) * .1 + min(prof_data)
    scattered=False
    if (np.nanmin(flags) > bot_prof_min) or ( not np.isnan(flags).any() ):
        logger.warning("The profile is highly scattered. S/N estimate cannot be calculated")
        scattered=True
        sn = sn_e = None
    else:
        #prof_e = 500. #this is when it's not normalised
        prof_e = 0.0005 #this is an approximation
        non_pulse_bins = 0
        #work out the above parameters
        for i, _ in enumerate(prof_data):
            if not np.isnan(flags[i]):
                non_pulse_bins += 1
        sigma_e = sigma / np.sqrt(2 * non_pulse_bins - 2)
        #now calc S/N
        sn = max(prof_data)/sigma
        sn_e = sn * np.sqrt(prof_e/max(prof_data)**2 + (sigma_e/sigma)**2)
        logger.debug("max prof: {} +/- {} ".format(max(prof_data), prof_e ))
        logger.debug("sigma   : {} +/- {} ".format(sigma, sigma_e ))
        logger.debug("sn      : {} +/- {} ".format(sn, sn_e ))

    return [sn, sn_e, scattered]

#---------------------------------------------------------------
def analyse_pulse_prof(prof_data, period, alpha=3):
    """
    Estimates the signal to noise ratio and many other properties from a pulse profile.
    Based on code oringally writted by Nick Swainston.
    This is the old version of 'est_sn_from_prof' but is useful when we can't fit gaussians

    Parameters:
    -----------
    prof_data: list
        A list of floats that contains the pulse profile.
    period: float
        The pulsar's period in ms
    alpha: float
        OPTIONAL - The alpha value to use when clipping using sigmaClip(). Default: 3

    Returns:
    --------
    prof_dict: dictionary
        contains keys:
        sn: float
            The estimated signal to noise ratio
        u_sn: float
            The estimated signal to noise ratio's its uncertainty
        flags: list
            A list of flagged data points
        w_equiv_bins: float
            The equivalent width of the profile measured in bins
        u_w_equiv_bins: float
            The uncertaintiy in w_equiv_bins
        w_equiv_ms: float
            The equivalent width of the profile measured in ms
        u_w_equiv_ms: float
            The uncertainty in w_equiv_ms
        scattering: float
            The scattering width in ms
        u_scattering: float
            The uncertainty in the scattering width in ms
        scattered: boolean
            When true, the profile is highly scattered
    """
    prof_dict = {}
    nbins = len(prof_data)
    #centre the profile around the max
    shift = -int(np.argmax(prof_data))+int(nbins)//2
    prof_data = np.roll(prof_data, shift)

    #find sigma and check if profile is scattered
    sigma, flags = sigmaClip(prof_data, alpha=alpha, tol=0.01, ntrials=100)
    check_clip(flags)

    bot_prof_min = (max(prof_data) - min(prof_data)) * .1 + min(prof_data)
    scattered=False
    if (np.nanmin(flags) > bot_prof_min) or ( not np.isnan(flags).any() ):
        logger.info("The profile is highly scattered. S/N estimate cannot be calculated")
        scattered=True
        #making a new profile with the only bin being the lowest point
        prof_min_i = np.argmin(prof_data)
        flags = []
        for fi, _ in enumerate(prof_data):
            if fi == prof_min_i:
                flags.append(prof_data[fi])
            else:
                flags.append(np.nan)

        flags = np.array(flags)
        prof_data -= min(prof_data)
        #Assuming width is equal to pulsar period because of the scattering
        w_equiv_ms = period
        u_w_equiv_ms = period/nbins
        sn = None
        u_sn = None
    else:
        u_prof = 500. #this is an approximation
        pulse_width_bins = 0
        non_pulse_bins = 0
        p_total = 0.
        u_p = 0.
        #work out the above parameters
        for i, data in enumerate(prof_data):
            if np.isnan(flags[i]):
                pulse_width_bins += 1
                p_total += data
                u_p = np.sqrt(u_p**2 + u_prof**2)
            else:
                non_pulse_bins += 1
        u_simga = sigma / np.sqrt(2 * non_pulse_bins - 2)

        #now calc S/N
        sn = max(prof_data)/sigma
        u_sn = sn * np.sqrt(u_prof/max(prof_data)**2 + (u_simga/sigma)**2)

    if not scattered:
        off_pulse_mean = np.nanmean(flags)
        prof_data -= off_pulse_mean
        flags -= off_pulse_mean

        prof_max = max(prof_data)
        w_equiv_bins = p_total / prof_max
        w_equiv_ms = w_equiv_bins / nbins * period # in ms
        u_w_equiv_bins = np.sqrt(p_total /prof_max)**2 +\
                                   (p_total * u_prof / (prof_max)**2)**2
        u_w_equiv_ms = u_w_equiv_bins / nbins * period # in ms

    else:
        w_equiv_ms = period
        u_w_equiv_ms = period/nbins
        w_equiv_bins = w_equiv_ms/period*nbins
        u_w_equiv_bins = (u_w_equiv_ms/w_equiv_ms)*w_equiv_bins

    #calc scattering
    scat_height = max(prof_data) / 2.71828
    scat_bins = 0
    for p in prof_data:
        if p > scat_height:
            scat_bins = scat_bins + 1
    scattering = float(scat_bins + 1) * float(period) /1000. #in s
    u_scattering = 1. * float(period) /1000. # assumes the uncertainty is one bin

    prof_dict["sn"] = sn
    prof_dict["sn_e"] = u_sn
    prof_dict["flags"] = flags
    prof_dict["w_equiv_bins"] = w_equiv_bins
    prof_dict["w_equiv_bins_e"] = u_w_equiv_bins
    prof_dict["w_equiv_ms"] = w_equiv_ms
    prof_dict["w_equiv_ms_e"] = u_w_equiv_ms
    prof_dict["scattering"] = scattering
    prof_dict["scattering_e"] = u_scattering
    prof_dict["scattered"] = scattered

    return prof_dict

def auto_analyse_pulse_prof(prof_data, period):
    """
    Automatically finds the best alpha value to use for analyse_pulse_prof() and returns the best resulting dictionary

    Parameters:
    -----------
    prof_data: list
        A list of floats that contains the pulse profile.
    preiod: float
        The period of the pulsar

    Returns:
    --------
    fit_dict: dictionary
        contains keys:
        sn: float
            The estimated signal to noise ratio
        u_sn: float
            The estimated signal to noise ratio's its uncertainty
        flags: list
            A list of flagged data points
        w_equiv_bins: float
            The equivalent width of the profile measured in bins
        u_w_equiv_bins: float
            The uncertaintiy in w_equiv_bins
        w_equiv_ms: float
            The equivalent width of the profile measured in ms
        u_w_equiv_ms: float
            The uncertainty in w_equiv_ms
        scattering: float
            The scattering width in ms
        u_scattering: float
            The uncertainty in the scattering width in ms
        scattered: boolean
            When true, the profile is highly scattered

    """
    if not isinstance(period, float):
        period = float(period)

    alphas = np.linspace(1, 5, 9)
    attempts_dict = {}

    loglvl = logger.level
    logger.setLevel(logging.WARNING) #squelch logging for the loop

    #loop over the gaussian evaluation fucntion, excepting in-built errors
    for alpha in alphas:
        try:
            prof_dict = analyse_pulse_prof(prof_data, period, alpha=alpha)
            attempts_dict[alpha] = prof_dict
        except(LittleClipError, LargeClipError, NoComponentsError, ProfileLengthError) as e:
            logger.setLevel(loglvl)
            logger.info(e)
            logger.info("Skipping alpha value: {}".format(alpha))
            logger.setLevel(logging.WARNING) #squelch logging for the loop
    logger.setLevel(loglvl)

    #Evaluate the best profile based on the SN error.
    sne = []
    sne_alphas = []
    scattered_trials = []
    if attempts_dict:
        for alpha_key in attempts_dict.keys():
            scattered_trials.append(attempts_dict[alpha_key]["scattered"])
            if not attempts_dict[alpha_key]["scattered"]:
                sne.append(attempts_dict[alpha_key]["sn_e"])
                sne_alphas.append(alpha_key)

    if sne: #there is an SN estimate available. Only look through these analyses as candidates
        best_sne = min(sne)
        best_alpha = sne_alphas[sne.index(best_sne)]
        fit_dict = attempts_dict[best_alpha]
    elif scattered_trials:
        mse = []
        scattered_alphas = []
        for alpha_key in attempts_dict.keys():
            if attempts_dict[alpha_key]["scattered"]:
                #use mean square of width errors as a metric for the best fit
                mse.append(np.sqrt(attempts_dict[alpha_key]["w_equiv_bins_e"]**2 + attempts_dict[alpha_key]["scattering_e"]**2))
                scattered_alphas.append(alpha_key)
        best_mse = min(mse)
        best_alpha = scattered_alphas[mse.index(best_mse)]
        fit_dict = attempts_dict[best_alpha]

    if not attempts_dict: #sometimes things go wrong :/
        logger.error("Profile could not be fit. Returning empty dictionary!")
        return {}

    logger.info("Best profile analysis using an alpha value of {}".format(best_alpha))

    return fit_dict

#---------------------------------------------------------------
def integral_multi_gauss(*params):
    y=0
    for i in range(0, len(params), 3):
        a = params[i]
        c = params[i+2]
        y = y + a*c*np.sqrt(2*np.pi)
    return y

def multi_gauss(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        a = params[i]
        b = params[i+1]
        c = params[i+2]
        y = y +  a * np.exp( -(((x-b)**2) / (2*c**2)) )
    return y

def multi_gauss_ddx(x, *params):
    #derivative of gaussian
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        a = params[i]
        b = params[i+1]
        c = params[i+2]
        y = y - a/c**2 * (x - b) * np.exp( -(((x-b)**2) / (2*c**2)) )
    return y

def multi_gauss_d2dx2(x, *params):
    #double derivative of gaussian
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        a = params[i]
        b = params[i+1]
        c = params[i+2]
        y = y + (multi_gauss(x, a, b, c) / c**2) * (((x - b)**2)/(c**2) - 1)
    return y

def partial_gauss_dda(x, a, b, c):
        return np.exp((-(b - x)**2)/(2*c**2))
def partial_gauss_ddb(x, a, b, c):
        return a*(x - b) * np.exp((-(b - x)**2)/(2*c**2))/c**2
def partial_gauss_ddc(x, a, b, c):
        return a*(x - b)**2 * np.exp((-(b - x)**2)/(2*c**2))/c**3

#---------------------------------------------------------------
def jacobian_slope(x, *params):
    """
    Evaluates the Jacobian matrix of a gaussian slope at a single point, x

    Parameters:
    -----------
    x: float
        The point to evaluate
    *params: list
        A list containing three parameters per gaussian component in the order: Amp, Mean, Width

    Returns:
    --------
    J: numpy.matrix
        The Jacobian matrix
    """
    def dda(a, b, c, x):
        return -multi_gauss(x, a, b, c) * (x - b)/(c**2)/a
    def ddb(a, b, c, x):
        return multi_gauss(x, a, b, c) * (1 - (x - b)**2/(c**2))/c**2
    def ddc(a, b, c, x):
        return multi_gauss(x, a, b, c) * (x - b)/(c**3) * (2 - (x-b)**2/(c**2))
    J = []
    for i in range(0, len(params), 3):
        a = params[i]
        b = params[i+1]
        c = params[i+2]
        mypars = [a, b, c, x]
        J.append(dda(*mypars))
        J.append(ddb(*mypars))
        J.append(ddc(*mypars))
    J = np.asmatrix(J)
    return J

#---------------------------------------------------------------
def plot_fit(plot_name, y, fit, popt, maxima=None, maxima_e=None):
    x = np.linspace(0, len(y)-1, len(y))
    plt.figure(figsize=(30, 18))

    for j in range(0, len(popt), 3):
        z = multi_gauss(x, *popt[j:j+3])
        plt.plot(x, z, "--", label="Gaussian Component {}".format(int((j+3)/3)))
    if maxima:
        for i, mx in enumerate(maxima):
            plt.axvline(x=(mx + maxima_e[i]), ls=":", lw=2, color="gray")
            plt.axvline(x=(mx - maxima_e[i]), ls=":", lw=2, color="gray")

    plt.title(plot_name.split("/")[-1].split(".")[0], fontsize=22)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlim(0, len(y))
    plt.xlabel("Bins", fontsize=20)
    plt.ylabel("Intensity", fontsize=20)

    plt.plot(x, y, label="Original Profile", color="black")
    plt.plot(x, fit, label="Gaussian Model", color="red")
    plt.legend(loc="upper right", prop={'size': 16})
    plt.savefig(plot_name)
    plt.close()

#---------------------------------------------------------------
def fit_gaussian(profile, max_N=6, min_comp_len=0, plot_name=None, alpha=3.):
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
    min_comp_len: float
        OPTIONAL - Minimum length of a component to be considered real. Measured in bins. Default: 0
    plot_name: string
        OPTIONAL - If not none, will make a plot of the best fit with this name. Default: None
    alpha: float
        OPTIONAL - The alpha value to be used in sigmaClip(). Default: 3

    Returns:
    --------
    [fit, redchisq, best_bic, popt, pcov]: list
        fit: list
            The data containing the multi-component gaussian fit to the input profile
        redchisq: float
            The reduced chi-sqaured value of the fit
        best_bic: float
            The bayesian information criterion for the fit
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

    #Take noise mean and normalize the profile and check the clipped profile
    _, clipped = sigmaClip(profile, alpha=alpha)
    check_clip(clipped)

    y = np.array(profile) - np.nanmean(np.array(clipped))
    max_y = max(y)
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
            comp_centres.append(np.mean(comp_idx[key]))
            comp_max.append(max(comp_dict[key])*0.5)
            comp_width.append((max(comp_idx[key])-min(comp_idx[key])))
    centre_guess = iter(comp_centres)
    width_guess=iter(comp_width)
    max_guess=iter(comp_max)

    n_comps=len(comp_dict.keys())
    logger.debug("Number of profile components: {0} ({1})".format(n_comps, comp_centres[:n_comps]))

    #Fit from 1 to max_N gaussians to the profile. Evaluate profile fit using bayesian information criterion
    x=np.linspace(0, len(y)-1, len(y))
    bounds_arr=[[],[]]
    guess = []
    fit_dict = {}

    for num in range(1, max_N):
        guess += [next(max_guess), next(centre_guess), next(width_guess)]
        bounds_arr[0].append(0)
        bounds_arr[0].append(0)
        bounds_arr[0].append(0)
        bounds_arr[1].append(max(y))
        bounds_arr[1].append(len(y))
        bounds_arr[1].append(len(y))
        bounds_tuple=(tuple(bounds_arr[0]), tuple(bounds_arr[1]))
        popt, pcov = curve_fit(multi_gauss, x, y, bounds=bounds_tuple, p0=guess, maxfev=100000)
        fit = multi_gauss(x, *popt)
        chisq = chsq(y, fit, noise_std)
        #Bayesian information criterion for gaussian noise
        k = 3*(num+1)
        bic = chisq + k*np.log(len(y))
        fit_dict[str(num+1)]={"popt":[], "pcov":[], "fit":[], "chisq":[], "bic":[]}
        fit_dict[str(num+1)]["popt"] = popt
        fit_dict[str(num+1)]["pcov"] = pcov
        fit_dict[str(num+1)]["fit"] = fit
        fit_dict[str(num+1)]["redchisq"] = chisq/(len(y)-1)
        fit_dict[str(num+1)]["bic"] = bic
        logger.debug("Reduced chi squared for               {0} components: {1}".format(num+1, fit_dict[str(num+1)]["redchisq"]))
        logger.debug("Bayesian Information Criterion for    {0} components: {1}".format(num+1, fit_dict[str(num+1)]["bic"]))

    #Find the best fit according to the BIC
    best_bic = np.inf
    best_fit = None
    for n_components in fit_dict.keys():
        if fit_dict[n_components]["bic"] < best_bic:
            best_bic = fit_dict[n_components]["bic"]
            best_fit = n_components
    logger.info("Fit {0} gaussians for a reduced chi sqaured of {1}".format(best_fit, fit_dict[best_fit]["redchisq"]))
    popt = fit_dict[best_fit]["popt"]
    pcov = fit_dict[best_fit]["pcov"]
    fit = fit_dict[best_fit]["fit"]
    redchisq = fit_dict[best_fit]["redchisq"]

    return [fit, redchisq, best_bic, popt, pcov, comp_dict, comp_idx]

#---------------------------------------------------------------
def prof_eval_gfit(profile, max_N=6, ignore_threshold=None, min_comp_len=None, plot_name=None, alpha=3., period=None):
    """
    Fits multiple gaussians to a profile and subsequently finds W10, W50, Weq and maxima

    Parameters:
    -----------
    profile: list
        The pulse profile to evaluate
    max_N: int
        OPTIONAL - The maximum number of gaussian components to attempt to fit. Default: 6
    ignore_threshold: float
        OPTIONAL -  Maxima with values below this number will be ignored. If none, will use 3*noise. Default: None
    min_comp_len: float
        OPTIONAL - Minimum length of a component to be considered real. Measured in bins. If None, will use 1% of total profile lengths + 2, max 50. Default: None
    plot_name: string
        OPTIONAL - If not none, will make a plot of the best fit with this name. Default: None
    alpha: float
        OPTIONAL - The alpha value passed to the sigmaClip() function. Default: 3
    period: float
        OPTIONAL - The puslar's period in ms. If not none, will attempt a S/N calculation. Default: None

    Returns:
    --------
    fit_dict: dictionary
        contains the following keys:
        W10: float
            The W10 width of the profile measured in phase
        W10_e: float
            The uncertainty in the W10
        W50: float
            The W50 width of the profile measured in phase
        W50_e: float
            The uncertainty in the W50
        Weq: float
            The equivalent width of the profile measured in phase
        Weq_e: float
            The uncertainty in the equivalent width
        Wscat: float
            The scattering width of the profile measured in phase
        Wscat_e: float
            The uncertainty in the scattering width
        maxima: list
            A lost of floats corresponding to the bin location of each maximum point
        maxima_e: list
            A list of floats, each correspinding to the error of the maxima of the same index. Measured in bins
        redchisq: float
            The reduced chi sqared of the fit
        num_gauss: int
            The number of gaussian components used in the best fit
        bic: float
            The Bayesian Information Criterion for the best fit
        gaussian_params: list
            A list of length 3*N there N is num_gauss. Each set of 3 parameters corresponds to the amp, centre and width of a guassian component
        cov_mat: np.matrix
            The covariance matrix from the fit
        comp_dict: dictionary
            dict["component_x"] contains an array of the component x
        comp_idx: dictionary
            dict["component_x"] contains an array of indexes of the original profile corresponding to component x
        alpha: float
            The alpha value used in sigmaClip()
        profile: list
            The input profile
        fit: list
            The best fit made into a list form
        sn: float
            The estimated signal to noise ratio, obtained from the profile. Will be None is period unsupplied
        sn_e: float
            The uncertainty in sn. Will be None is period unsupplied
        scattered: boolean
            True is the profile is scattered. Will be None is period unsupplied
    """
    #initialize minimum component length and ignore threshold
    if min_comp_len is None:
        min_comp_len = int(len(profile)/100 + 0.5) + 2
        if min_comp_len > 100:
            min_comp_len = 100
    if min_comp_len < 3:
        min_comp_len = 3

    #Normalize, find the std
    y = np.array(profile)/max(profile)
    noise_std, clipped = sigmaClip(y, alpha=alpha)
    check_clip(clipped)

    if ignore_threshold is None:
        ignore_threshold = 3 * noise_std

    y = y - np.nanmean(clipped)
    y = y/max(y)
    proflen = len(y)

    #fit gaussians
    fit, chisq, bic, popt, pcov, comp_dict, comp_idx = fit_gaussian(y, max_N=max_N, min_comp_len=min_comp_len, alpha=alpha)
    fit = np.array(fit)
    n_rows, _ = np.shape(pcov)
    num_gauss = n_rows/3

    #Find widths + error
    W10, W50, Weq, Wscat, W10_e, W50_e, Weq_e, Wscat_e = find_widths(y, popt, pcov, alpha=alpha)
    #Convert from bins to phase
    W10 = W10/proflen
    W50 = W50/proflen
    Weq = Weq/proflen
    Wscat = Wscat/proflen
    W10_e = W10_e/proflen
    W50_e = W50_e/proflen
    Weq_e = Weq_e/proflen
    Wscat_e = Wscat_e/proflen

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

    _, maxima, _, maxima_e = find_minima_maxima_gauss(popt, pcov, len(fit), min_comp_len=min_comp_len)

    #estimate SN
    if period:
        sn, sn_e, scattered = est_sn_from_prof(y, period, alpha=alpha)
    else:
        sn = sn_e = scattered = None

    #plotting
    if plot_name:
        plot_fit(plot_name, y, fit, popt, maxima=maxima, maxima_e=maxima_e)

    fit_dict = {"W10":W10, "W10_e":W10_e, "W50":W50, "W50_e":W50_e, "Wscat":Wscat, "Wscat_e":Wscat_e,\
                "Weq":Weq, "Weq_e":Weq_e, "maxima":maxima, "maxima_e":maxima_e, "redchisq":chisq,\
                "num_gauss":num_gauss, "bic":bic, "gaussian_params":popt, "cov_mat":pcov, "comp_dict":comp_dict,\
                "comp_idx":comp_idx, "alpha":alpha, "profile":y, "fit":fit, "sn":sn, "sn_e":sn_e, "scattered":scattered}

    logger.info("W10:                   {0} +/- {1}".format(W10, W10_e))
    logger.info("W50:                   {0} +/- {1}".format(W50, W50_e))
    logger.info("Wscat:                 {0} +/- {1}".format(Wscat, Wscat_e))
    logger.info("Weq:                   {0} +/- {1}".format(Weq, Weq_e))
    logger.info("Maxima:                {0}".format(maxima))
    logger.info("Maxima error:          {0}".format(maxima_e))
    if fit_dict["sn"]:
        logger.info("S/N estimate:          {0} +/- {1}".format(fit_dict["sn"], fit_dict["sn_e"]))

    return fit_dict

def auto_gfit(profile, max_N=6, plot_name=None, ignore_threshold=None, min_comp_len=None, period=None, cliptype="regular"):
    """
    runs the gaussian fit evaluation for a range of values of alpha. This is necessary as there is no way to know
    a priori which alpha to use beforehand. Alpha is the input for sigmaClip() and can be interpreted as the level
    of verbosity in clipping.

    Parameters:
    -----------
    profile: list
        A list containing the pulse profile to evaluate
    max_N: int
        OPTIONAL - The maximum number of gaussian components to use when fitting
    plot_name: string
        OPTIONAL - If not none, will make a plot of the best fit with this name. Default: None
    ignore_threshold: float
        OPTIONAL -  Maxima with values below this number will be ignored. If none, will use 3*noise. Default: None
    min_comp_len: float
        OPTIONAL - Minimum length of a component to be considered real. Measured in bins. If none, will use 1% of total profile lengths + 2, max 50. Default: None

    Returns:
    --------
    fit_dict: dictionary
        The dictionary of the best fit from prof_eval_gfit
    """
    if len(profile)<100:
        raise ProfileLengthError("Profile must have length > 100")

    if cliptype == "regular":
        alphas = np.linspace(1, 5, 9)
    elif cliptype == "noisy":
        alphas = np.linspace(1, 3, 17)
    elif cliptype == "verbose":
        alphas = np.linspace(1, 5, 33)
    else:
        raise ValueError("cliptype not recognised. Options are: 'regular', 'noisy' or 'verbose'.")

    attempts_dict = {}

    loglvl = logger.level
    logger.setLevel(logging.WARNING) #squelch logging for the loop

    #loop over the gaussian evaluation fucntion, excepting in-built errors
    for alpha in alphas:
        try:
            prof_dict = prof_eval_gfit(profile, max_N=6, ignore_threshold=ignore_threshold, min_comp_len=min_comp_len, alpha=alpha, period=period)
            attempts_dict[alpha] = prof_dict
        except(LittleClipError, LargeClipError, NoComponentsError, ProfileLengthError):
            logger.setLevel(loglvl)
            logger.info("Skipping alpha value: {}".format(alpha))
            logger.setLevel(logging.WARNING) #squelch logging for the loop
    logger.setLevel(loglvl)

    #Evaluate the best profile based on reduced chi-squared.
    chi_diff = []
    alphas = []
    for alpha_key in attempts_dict.keys():
        chi_diff.append(abs(1 - attempts_dict[alpha_key]["redchisq"]))
        alphas.append(alpha_key)

    if not chi_diff: #sometimes things go wrong :/
        raise NoFitError("No suitable profile fit could be found!")

    best_chi = min(chi_diff)
    best_alpha = alphas[chi_diff.index(best_chi)]
    fit_dict = attempts_dict[best_alpha]

    if plot_name:
        plot_fit(plot_name, fit_dict["profile"], fit_dict["fit"], fit_dict["gaussian_params"], maxima=fit_dict["maxima_e"], maxima_e=fit_dict["maxima"])

    logger.info("### Best fit results ###")
    logger.info("Best model found with BIC of {0} and reduced Chi of {1} using an alpha value of {2}"\
                .format(fit_dict["bic"], fit_dict["redchisq"], best_alpha))
    logger.info("W10:                   {0} +/- {1}".format(fit_dict["W10"], fit_dict["W10_e"]))
    logger.info("W50:                   {0} +/- {1}".format(fit_dict["W50"], fit_dict["W50_e"]))
    logger.info("Wscat:                 {0} +/- {1}".format(fit_dict["Wscat"], fit_dict["Wscat_e"]))
    logger.info("Weq:                   {0} +/- {1}".format(fit_dict["Weq"], fit_dict["Weq_e"]))
    logger.info("Maxima:                {0}".format(fit_dict["maxima"]))
    logger.info("Maxima error:          {0}".format(fit_dict["maxima_e"]))
    if fit_dict["sn"]:
        logger.info("S/N estimate:          {0} +/- {1}".format(fit_dict["sn"], fit_dict["sn_e"]))

    return fit_dict
