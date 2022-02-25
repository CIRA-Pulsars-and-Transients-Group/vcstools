import numpy as np
import sys
import re
import subprocess
import re
import logging
from astropy.time import Time
from scipy.interpolate import UnivariateSpline
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import itertools

from vcstools.config import load_config_file
from vcstools.stickel import Stickel

logger = logging.getLogger(__name__)

# Custom Errors


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


class BadFitError(Exception):
    """Raise when a particular fit is not suitable"""


def subprocess_pdv(archive, outfile="archive.txt", pdvops="-FTt"):
    """Runs the pdv commnand from PSRCHIVE as a python subprocess.
    NOTE: Requires singularity loaded in environment (module load singularity)

    Parameters
    ----------
    archive : `str`
        The name of the archive file to run the command on
    outfile : `str`, optional
        The name of the text file to write the output to. |br| Default: archive.txt.
    pdvops : `str`, optional
        Additional options for the pdv. |br| Default: -FTt.
    """
    comp_config = load_config_file()
    with open(outfile, 'w+') as f:
        commands = [comp_config["prschive_container"]]
        commands.append("pdv")
        commands.append(pdvops)
        commands.append(archive)
        a = subprocess.check_output(commands)
        f.write(a.decode("utf-8"))


#---------------------------------------------------------------
def get_from_bestprof(file_loc,
                      pointing_input=None):
    """Get info from a bestprof file

    Parameters
    ----------
    file_loc : `str`
        The path to the bestprof file.
    pointing_input : `str`
        If can not find a pointing in the bestprof file, will assume this is the pointing.
        |br| Default: `None`.

    Returns
    -------
    obsid : `int`
        The observation ID.
    pulsar : `str`
        The name of the pulsar.
    dm : `float`
        The dispersion measure of the pulsar.
    period : `float`
        The period of the pulsar in ms.
    period_uncer : `float`
        The uncertainty in the period measurement.
    sigma : `float`
        The sigma (similar to signal to noise ratio).
    obsstart : `int`
        The beginning time of the observation.
    obslength : `float`
        The length of the observation in seconds.
    profile: list
        A list of floats containing the profile data.
    bin_num : `int`
        The number of bins in the profile.
    """
    with open(file_loc, "r") as bestprof:
        lines = bestprof.readlines()
        # Find the obsid by finding a 10 digit int in the file name
        obsid = re.findall(r'(\d{10})', lines[0])[0]
        try:
            obsid = int(obsid)
        except ValueError:
            obsid = None

        # Get a pointing from the input fits file name
        ra, dec = lines[0].split("_ch")[0].split("_")[-2:]
        pointing = "{}_{}".format(ra, dec)
        if ":" not in pointing:
            if pointing_input is None:
                logger.error("No pointing found, please input one. Exiting.")
                sys.exit(0)
            else:
                pointing = pointing_input

        pulsar = str(lines[1].split("_")[-1][:-1])
        if not (pulsar.startswith('J') or pulsar.startswith('B')):
            pulsar = 'J{0}'.format(pulsar)

        dm = lines[14][22:-1]

        period = lines[15][22:-1]
        period, period_uncer = period.split('  +/- ')

        sigma = float(lines[13].split("~")[-1].split(" sigma")[0])

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

    return [obsid, pulsar, dm, period, period_uncer, sigma, obsstart, obslength, profile, bin_num, pointing]


def get_from_ascii(file_loc):
    """Retrieves the profile from an ascii file.

    Parameters
    ----------
    file_loc : `str`
        The location of the ascii file.

    Returns
    -------
    profile : `list`
        A list of floats containing the profile data.
    len(profile) : `int`
        The number of bins in the profile.
    """

    f = open(file_loc)
    lines = iter(f.readlines())
    next(lines)  # skip first line
    f.close()
    profile = []
    for line in lines:
        thisline = line.split()
        profile.append(float(thisline[3]))

    return [profile, len(profile)]


def get_stokes_from_ascii(file_loc):
    """Retrieves the all stokes components from an ascii file.

    Parameters
    ----------
    file_loc : `str`
        The location of the ascii file.

    Returns
    -------
    [I, Q, U, V, len(profile)] : `list`

        I : `list`
            Stokes I.
        Q : `list`
            Stokes Q.
        U : `list`
            Stokes U.
        V : `list`
            Stokes V.
        len(profile) : `int`
            The number of bins in the profile.
    """
    f = open(file_loc)
    lines = iter(f.readlines())
    f.close()
    next(lines)  # skip first line
    I = []
    Q = []
    U = []
    V = []
    for line in lines:
        thisline = line.split()
        I.append(float(thisline[3]))
        Q.append(float(thisline[4]))
        U.append(float(thisline[5]))
        V.append(float(thisline[6]))

    return [I, Q, U, V, len(I)]


def sigmaClip(data, alpha=3., tol=0.1, ntrials=10):
    """Sigma clipping operation:
    Compute the data's median, m, and its standard deviation, sigma.
    Keep only the data that falls in the range (m-alpha*sigma,m+alpha*sigma) for some value of alpha, and discard everything else.
    This operation is repeated ntrials number of times or until the tolerance level is hit.

    Parameters
    ----------
    data : `list`
        A list of floats - the data to clip.
    alpha : `float`, optional
        Determines the number of sigmas to use to determine the upper and lower limits. |br| Default: 3.
    tol : `float`, optional
        The fractional change in the standard deviation that determines when the tolerance is hit. |br| Default: 0.1.
    ntrials : `int`, optional
        The maximum number of times to apply the operation. |br| Default: 10.

    Returns
    -------
    oldstd : `float`
        The std of the clipped data.
    x : `list`
        The data list that contains only noise, with nans in place of 'real' data.
    """
    x = np.copy(data)
    oldstd = np.nanstd(x)
    # When the x[x<lolim] and x[x>hilim] commands encounter a nan it produces a
    # warning. This is expected because it is ignoring flagged data from a
    # previous trial so the warning is supressed.
    old_settings = np.seterr(all='ignore')
    for trial in range(ntrials):
        median = np.nanmedian(x)
        lolim = median - alpha * oldstd
        hilim = median + alpha * oldstd
        x[x < lolim] = np.nan
        x[x > hilim] = np.nan

        newstd = np.nanstd(x)
        tollvl = (oldstd - newstd) / newstd

        if tollvl <= tol:
            logger.debug(f"Took {trial+1} trials to reach tolerance")
            np.seterr(**old_settings)
            return oldstd, x

        if trial + 1 == ntrials:
            logger.debug(
                "Reached number of trials without reaching tolerance level")
            np.seterr(**old_settings)
            return oldstd, x

        oldstd = newstd


def sn_calc(profile, on_pulse_bool, noise_std, w_equiv_bins, sn_method="eqn_7.1"):
    """Calculates the flux desnity using the siganl to noise (The maximum/peak
    signal divided by the STD of the noise) and the radiometer equation (see eqn
    A 1.21 of the pulsar handbook)

    Parameters
    ----------
    profile : `list`
        The normalised profile of the pulsar.
    on_pulse_bool : `list`
        A list of booleans that are `True` if the profile is currently on the on
        pulse and `False` if the porifle is currently on the off pulse.
    noise_std : `float`
        The standard deviation of the noise in the normalised profile.
    w_equiv_bins : `int`
        The equivalent width (in bins) of a top-hat pulse with the same are and
        peak height as the observed profile.
    sn_method : `str`, optional
        The method of calculating the signal to noise ratio out of ["eqn_7.1", "simple"]. Default "eqn_7.1".

        "Simple" uses 1/noise_std.

        "eqn_7.1" uses equation 7.1 of the pulsar handbook.

    Returns
    -------
    sn : `float`
        The signal to noise ratio of the profile.
    u_sn : `float`
        The uncertainty of the signal to noise ratio.
    """
    proflen = len(profile)
    n_on_pulse = np.count_nonzero(on_pulse_bool)
    n_off_pulse = proflen - n_on_pulse

    u_noise_std = noise_std / np.sqrt(2 * n_off_pulse -2)

    # Estimate SN
    if sn_method == "eqn_7.1":
        # Estimate SN using Equation 7.1 in the handbook of pulsar astronomy
        sn = np.sum(profile) / (noise_std * np.sqrt(w_equiv_bins))
        # Standard uncertainty propigation
        # uncertainty of each point in the profile is the noise std
        u_sum = np.sqrt(proflen * noise_std**2)
        # assume w_equiv_bins uncertainty is 1 bin
        u_sn = np.sqrt( (u_sum / np.sum(profile))**2 + (u_noise_std/noise_std)**2 + (1/(2 * np.sqrt(w_equiv_bins)))**2 ) * sn
    elif sn_method == "simple":
        sn = 1 / noise_std
        # Standard uncertainty propigation
        u_sn = sn**2 * u_noise_std
    else:
        logger.error(f"sn_method {sn_method} not found. Exiting")
        sys.exit(1)
    logger.debug(f"sn: {sn} +/- {u_sn}")

    return sn, u_sn


def analyse_pulse_prof(prof_data, period, alpha=3, scattering_threshold=0.7):
    """Estimates the signal to noise ratio and many other properties from a pulse profile.

    Parameters
    ----------
    prof_data : `list`
        A list of floats that contains the pulse profile.
    period : `float`
        The pulsar's period in ms.
    alpha : `float`, optional
        The alpha value to use when clipping using :py:meth:`vcstools.prof_utils.sigmaClip`. |br| Default: 3.

    Returns
    -------
    prof_dict : `dict`
        The analysed profile's dictionary which contains the following keys:

        ``"sn"``
            The estimated signal to noise ratio (`float`).
        ``"sn_e"``
            The estimated signal to noise ratio's its uncertainty (`float`).
        ``"profile"``
            A list containting the noramlised pulse profile (`list`).
        ``"on_pulse_bool"``
            A list of booleans that are `True` if the profile is currently on the on
            pulse and `False` if the porifle is currently on the off pulse (`list`).
        ``"Weq"``
            The equivalent width of the profile measured in phase (`float`).
        ``"Weq_e"``
            The uncertaintiy in Weq_e (`float`).
        ``"w_equiv_bins"``
            The equivalent width of the profile measured in bins (`float`).
        ``"w_equiv_bins_e"``
            The uncertaintiy in w_equiv_bins (`float`).
        ``"w_equiv_ms"``
            The equivalent width of the profile measured in ms (`float`).
        ``"w_equiv_ms_e"``
            The uncertainty in w_equiv_ms (`float`).
        ``"scattering"``
            The scattering width in phase (`float`).
        ``"scattering_e"``
            The uncertainty in the scattering width in phase (`float`).
        ``"scattered"``
            When true, the profile is highly scattered (`boolean`).
        ``"flags"``
            A list of flagged data points (`list`).
    """
    prof_dict = {}
    nbins = len(prof_data)
    # centre the profile around the max and normalise
    shift = -int(np.argmax(prof_data))+int(nbins)//2
    prof_data = np.roll(prof_data, shift)
    prof_data = normamlise_prof(prof_data)

    # find sigma and check if profile is scattered
    sigma, flags = sigmaClip(prof_data, alpha=alpha, tol=0.01, ntrials=100)
    check_clip(flags)

    # Make a boolean list of if the part of the profile is on
    on_pulse_bool = []
    for opp in flags:
        if np.isnan(opp):
            on_pulse_bool.append(True)
        else:
            on_pulse_bool.append(False)

    pulse_width_bins = sum(on_pulse_bool)
    non_pulse_bins = pulse_width_bins - nbins

    # Check if scattered and calc SN
    if pulse_width_bins / nbins > scattering_threshold:
        logger.info("The profile is highly scattered. S/N estimate cannot be calculated")
        scattered = True
        sn = None
        u_sn = None
    else:
        scattered = False
        sn, u_sn = sn_calc(prof_data, on_pulse_bool, pulse_width_bins, sigma)

    # Centre profile on noise mean and renormalise
    #prof_data -= np.nanmean(flags)
    #prof_data /= max(prof_data)

    # Calculate equivelent widths
    # sum over entire profile to get are under pulse
    p_total = sum(prof_data)
    # uncertainty of each point in the profile is the noise std
    u_p_total = np.sqrt(nbins * sigma**2)
    # width equivelent is area/height
    w_equiv_bins = p_total # / 1
    u_w_equiv_bins = u_p_total
    # convert to phase
    w_equiv_phase = w_equiv_bins / nbins
    u_w_equiv_phase = u_w_equiv_bins / nbins
    # convert to ms
    w_equiv_ms = w_equiv_phase * period
    u_w_equiv_ms = u_w_equiv_phase * period

    # calc scattering
    scat_height = max(prof_data) / 2.71828
    scat_bins = 0
    for p in prof_data:
        if p > scat_height:
            scat_bins = scat_bins + 1
    scattering = float(scat_bins + 1)  / nbins
    # assumes the uncertainty is one bin
    u_scattering = 1. / nbins

    prof_dict["sn"] = sn
    prof_dict["sn_e"] = u_sn
    prof_dict["profile"] = prof_data
    prof_dict["on_pulse_bool"] = on_pulse_bool
    prof_dict["noise_std"] = sigma
    prof_dict["noise_mean"] = np.nanmean(flags)

    prof_dict["Weq"] = w_equiv_phase
    prof_dict["Weq_e"] = u_w_equiv_phase
    prof_dict["w_equiv_bins"] = w_equiv_bins
    prof_dict["w_equiv_bins_e"] = u_w_equiv_bins
    prof_dict["w_equiv_ms"] = w_equiv_ms
    prof_dict["w_equiv_ms_e"] = u_w_equiv_ms
    prof_dict["Wscat"] = scattering
    prof_dict["Wscat_e"] = u_scattering
    prof_dict["scattered"] = scattered
    prof_dict["flags"] = flags

    return prof_dict


def auto_analyse_pulse_prof(prof_data, period):
    """Automatically finds the best alpha value to use for analyse_pulse_prof() and returns the best resulting dictionary.

    Parameters
    ----------
    prof_data : `list`
        A list of floats that contains the pulse profile.
    period : `float`
        The period of the pulsar in ms.

    Returns
    -------
    fit_dict : `dict`
        The analysed profile's dictionary which contains the following keys:

        ``"sn"``
            The estimated signal to noise ratio (`float`).
        ``"sn_e"``
            The estimated signal to noise ratio's its uncertainty (`float`).
        ``"flags"``
            A list of flagged data points (`list`).
        ``"w_equiv_bins"``
            The equivalent width of the profile measured in bins (`float`).
        ``"w_equiv_bins_e"``
            The uncertaintiy in w_equiv_bins (`float`).
        ``"w_equiv_ms"``
            The equivalent width of the profile measured in ms (`float`).
        ``"w_equiv_ms_e"``
            The uncertainty in w_equiv_ms (`float`).
        ``"scattering"``
            The scattering width in ms (`float`).
        ``"scattering_e"``
            The uncertainty in the scattering width in ms (`float`).
        ``"scattered"``
            When true, the profile is highly scattered (`boolean`).
    """
    if not isinstance(period, float):
        period = float(period)

    alphas = np.linspace(1, 5, 9)
    attempts_dict = {}

    loglvl = logger.level
    logger.setLevel(logging.WARNING)  # squelch logging for the loop

    # loop over the gaussian evaluation fucntion, excepting in-built errors
    for alpha in alphas:
        try:
            prof_dict = analyse_pulse_prof(prof_data, period, alpha=alpha)
            attempts_dict[alpha] = prof_dict
        except(LittleClipError, LargeClipError, NoComponentsError, ProfileLengthError) as e:
            logger.setLevel(loglvl)
            logger.info(e)
            logger.info(f"Skipping alpha value: {alpha}")
            logger.setLevel(logging.WARNING)  # squelch logging for the loop
    logger.setLevel(loglvl)

    # Evaluate the best profile based on the SN error.
    sne = []
    sne_alphas = []
    scattered_trials = []
    if attempts_dict:
        for alpha_key in attempts_dict.keys():
            scattered_trials.append(attempts_dict[alpha_key]["scattered"])
            if not attempts_dict[alpha_key]["scattered"]:
                sne.append(attempts_dict[alpha_key]["sn_e"])
                sne_alphas.append(alpha_key)

    if sne:  # there is an SN estimate available. Only look through these analyses as candidates
        best_sne = min(sne)
        best_alpha = sne_alphas[sne.index(best_sne)]
        fit_dict = attempts_dict[best_alpha]
    elif scattered_trials:
        mse = []
        scattered_alphas = []
        for alpha_key in attempts_dict.keys():
            if attempts_dict[alpha_key]["scattered"]:
                # use mean square of width errors as a metric for the best fit
                mse.append(np.sqrt(
                    attempts_dict[alpha_key]["w_equiv_bins_e"]**2 + attempts_dict[alpha_key]["scattering_e"]**2))
                scattered_alphas.append(alpha_key)
        best_mse = min(mse)
        best_alpha = scattered_alphas[mse.index(best_mse)]
        fit_dict = attempts_dict[best_alpha]

    if not attempts_dict:  # sometimes things go wrong :/
        logger.error("Profile could not be fit. Returning empty dictionary!")
        return {}

    logger.info(f"Best profile analysis using an alpha value of {best_alpha}")

    return fit_dict


# SigmaClip isn't perfect. Use these next function to check for bad clips
def check_clip(prof_to_check, toomuch=0.9, toolittle_frac=0., toolittle_absolute=4):
    """Determines whether a clipped profile from :py:meth:`vcstools.prof_utils.sigmaClip` has been appropriately clipped by checking the number of nans.
    Raises a LittleClipError or a LargeClipError if too little or too much of the data has been clipped respectively.

    Parameters
    ----------
    prof_to_check : `list`
        The clipped profile from :py:meth:`vcstools.prof_utils.sigmaClip`
    toomuch : `float`, optional
        The fraction of the clipped profile beyond which is considered overclipped. |br| Default: 0.9
    toolittle : `float`, optional
        The fraction of the clipped profile below which is considered underclipped. |br| Default: 0.
    toolittle_absolute : `int`, optional
        If a profile has this many or less on-pulse bins, it is deemed not sufficient. |br| Default: 4
    """
    num_nans = 0
    for i in prof_to_check:
        if np.isnan(i):
            num_nans += 1
    if num_nans <= toolittle_frac*len(prof_to_check) or num_nans <= toolittle_absolute:
        raise LittleClipError(
            "Not enough data has been clipped. Condsier trying a smaller alpha value when clipping.")
    elif num_nans >= toomuch*len(prof_to_check):
        raise LargeClipError(
            "A large portion of the data has been clipped. Condsier trying a larger alpha value when clipping.")


def normamlise_prof(x):
    """Normalises a profile, x, to minimum 0 and maximum 1.

    Parameters
    ----------
    x : `list`
        A list of pulse profile values.

    Returns
    -------
    norm_x : `list`
        A normalised profile.
    """
    x = np.array(x)
    min_x = min(x)
    max_x = max(x)
    norm_x = (x - min_x)/(max_x - min_x)
    return norm_x


def error_in_x_pos(x, y, sigma_y, x_pos):
    """Finds error in a position, x_pos, given an error in y
    derivatives are found numerically.

    Parameters
    ----------
    x : `list`
        A list that covers the x range (ie. np.linspace(0, 1, len(x))).
    y : `list`
        The y components of the x range.
    sigma_y : `float`
        The absolute error in y.
    x_pos : `int`
        The location in x that we want to find the error for.

    Returns
    -------
    x_er : `float`
        The error in the x position.
    """
    # Spline, get derivative and roots
    spline = UnivariateSpline(x, y, s=0, k=4)
    dy_spline = spline.derivative()
    derivative_profile = dy_spline(x)
    dydx = derivative_profile[x_pos]

    # Double derivative
    spline = UnivariateSpline(x, y, s=0, k=5)
    d2y_spline = spline.derivative().derivative()
    double_derivative_profile = d2y_spline(x)
    d2ydx2 = double_derivative_profile[x_pos]

    # Derived from a taylor expansion and solved using quadratic formula:
    # x = ( dydx +/- sqrt( dydx^2 +/- 2*d2ydx2 * sigma_y ) / d2ydx2
    # We only care about two of the four solutions because the other two are identical but negative
    deriv_er_plus = (dydx + np.sqrt(dydx**2 + 2*d2ydx2 * sigma_y)) / d2ydx2
    deriv_er_minus = (
        dydx + np.sqrt(dydx**2 - 2*d2ydx2 * sigma_y)) / d2ydx2
    x_er = abs(np.nanmean((deriv_er_plus, deriv_er_minus)))

    return x_er


# Alternative method to find the location and width of peaks
def estimate_components_onpulse(profile, l=1e-5, plot_name=None):
    """ Works by using stickel with lambda=1e-5 (default) to get a generic profile outline, then uses UnivariateSpline to
    find the min/maxima and estimate which of these are real depending on their relative amplitude.
    NOTE: The on-pulse estimates are read form LEFT to RIGHT. Meaning onpulse[0] can be a greater index than onpulse[1].
    This indicates that the on-pulse is wrapped over the end of the phase.

    Parameters
    ----------
    profile : `list`
        The pulse profile.
    l : `float`, optional
        The lambda value to use for Stickel regularisation. Don't touch unless you know what you're doing. |br| Default: 1e-5.
    plot_name : `str`, optional
        If supplied, will make a plot of the profile, its splined regularisation, maxima and best estimate of on-pulse.

    Returns
    -------
    comp_est_dict : `dict`
        Component estimation dictionary with keys:

        ``"maxima"``
            The maximum points of the splined profile (`list`).
        ``"underest_on_pulse"``
            List of tuples containing the understimated beginning and end of each profile component (`list`).
        ``"overest_on_pulse"``
            List of tuples containing the overestimated beginning and end of each profile component (`list`).
        ``"underest_off_pulse"``
            List of tuples containing the beginning and end of each off-pulse range from the underestimated on pulse (`list`).
        ``"overest_off_pulse"``
            List of tuples containing the beginning and end of each off-pulse range from the overestimated on pulse (`list`).
        ``"noise"``
            The standard deviation of the noise taken as the off-pulse region (`float`).
        ``"alpha"``
            The alpha value def to :py:meth:`vcstools.prof_utils.sigmaClip` to attain the initial on-pulse and noise estimate (`float`).
    """
    x = np.linspace(0, len(profile)-1, len(profile))
    # Normalise profile
    norm_prof = normamlise_prof(profile)

    # Do only a very light sweep of sigmaClip to get a first order noise estimate
    alpha, clip_noise = _profile_noise_initial_estimate(norm_prof)
    _, y = sigmaClip(norm_prof, alpha=alpha)
    noise_mean = np.nanmean(y)

    # Take mean from normalised profile
    norm_prof = norm_prof - noise_mean

    # Stickel takes a numpy array in the form [[x0, y0], [x1, y1], ..., [xn, yn]]
    xy_prof = []
    for i, val in enumerate(norm_prof):
        xy_prof.append([x[i], val])
    xy_prof = np.array(xy_prof)

    # Perform the stickel smooth
    stickel = Stickel(xy_prof)
    stickel.smooth_y(l)
    reg_prof = stickel.yhat

    # Spline, get derivative and roots
    spline = UnivariateSpline(x, reg_prof, s=0, k=4)
    dy_spline = spline.derivative()
    dy_roots = dy_spline.roots()
    dy_roots = [int(round(i)) for i in dy_roots]  # round to integer
    spline = UnivariateSpline(x, reg_prof, s=0, k=5)
    dy2_spline = spline.derivative().derivative()
    dy2_roots = dy2_spline.roots()
    dy2_roots = [int(round(i)) for i in dy2_roots]  # round to integer
    dy2_profile = dy2_spline(x)

    # Profile with spline applied
    splined_profile = spline(x)

    # Find which dy roots resemble real signal
    noise_trials = (2.0, 1.5, 1.0, 0.5, 0.25, 0.0)
    for i in noise_trials: # This represents how generous we are with what is real signal
        real_max = []
        false_max = []
        for root in dy_roots:
            amp_at_root = splined_profile[root]
            if amp_at_root < i*clip_noise:  # If the amp of the regularised profile is too small, it's not signal
                false_max.append(root)
            # If this is a min, we don't care about it
            elif dy2_profile[root] > 0:
                false_max.append(root)
            else:  # Otherwise, this is a real max
                real_max.append(root)
        if i == noise_trials[-1]:
            logger.warn("This looks like a very noisy profile. Results may be unreliable.")
        if real_max:
            break  # Break if anything has populated real_max


    ##########################################################################################################
    # We will create an over and under-estimate...

    # The overestimate will mark the regions between the real max and the flanking false minima of the spline.
    # This range should include the entire on-pulse region but also some noise

    # The understimate will mark the regions between the real max and its flanking inflections.
    # This range should include no noise but will miss some of the on-pulse
    ##########################################################################################################

    overest_on_pulse = []
    underest_on_pulse = []
    for root in real_max:
        # Underestimated on-pulse
        # Find where the index of root would be
        dy2_roots = np.append(dy2_roots, root)
        dy2_roots = sorted(dy2_roots)
        i = dy2_roots.index(root)
        dy2_roots.remove(root)
        try:
            underestimate = [round(dy2_roots[i-1]), round(dy2_roots[i])]
        except IndexError as e:
            # The next inflection is on the other side of the profile
            if i == 0:  # Left side
                underestimate = [round(dy2_roots[-1]), round(dy2_roots[i])]
            else:  # Right side
                underestimate = [round(dy2_roots[i-1]), round(dy2_roots[0])]
        underest_on_pulse.append(underestimate)

        # Overestimated on-pulse
        false_max = np.append(false_max, root)
        false_max = sorted(false_max)
        j = false_max.index(root)
        false_max.remove(root)
        # Append the overestimate
        try:
            overestimate = [round(false_max[j-1]), round(false_max[j])]
        except IndexError as e:
            # The next inflection is on the other side of the profile
            if j == 0:  # Left side
                overestimate = [round(false_max[-1]), round(false_max[j])]
            else:  # Right side
                overestimate = [round(false_max[j-1]), round(false_max[0])]
        # If the bounds are touching
        if overest_on_pulse and overestimate[0] == overest_on_pulse[-1][1]:
            overest_on_pulse[-1][1] = overestimate[1]
        # If the last element wraps around to the first
        elif overest_on_pulse and overestimate[1] == overest_on_pulse[0][0]:
            overest_on_pulse[0][0] = overestimate[0]
        else:
            overest_on_pulse.append(overestimate)

    # Remove any duplicate items because somehow this happens sometimes
    overest_on_pulse.sort()
    overest_on_pulse = list(k for k, _ in itertools.groupby(overest_on_pulse))
    underest_on_pulse.sort()
    underest_on_pulse = list(
        k for k, _ in itertools.groupby(underest_on_pulse))

    # Fill the on-pulse regions
    overest_on_pulse = _fill_on_pulse_region(
        overest_on_pulse, splined_profile, norm_prof, clip_noise)
    underest_on_pulse = _fill_on_pulse_region(
        underest_on_pulse, splined_profile, norm_prof, clip_noise)
    # Get the off-pulse pairs as well
    overest_off_pulse = get_off_pulse(overest_on_pulse)
    underest_off_pulse = get_off_pulse(underest_on_pulse)

    # Work out the noise std and mean using the oversetimated on-pulse region
    off_pulse_noise_std = np.nanstd(
        profile_region_from_pairs(norm_prof, overest_off_pulse))
    off_pulse_noise_mean = np.nanmean(
        profile_region_from_pairs(normamlise_prof(profile), overest_off_pulse))

    # Sort the maxima based on amplitude
    maxima_amps = []
    for m in real_max:
        maxima_amps.append(splined_profile[m])
    real_max = [i for _, i in sorted(zip(maxima_amps, real_max), reverse=True)]

    # plot if necessary
    if plot_name:
        plt.figure(figsize=(14, 10))
        plt.plot(x, norm_prof, label="Profile", color="cornflowerblue")
        plt.plot(x, reg_prof, label="Smoothed Profile", color="orange")
        plt.xlim((0, len(x)-1))
        for i in real_max:
            plt.axvline(x=i, ls=":", lw=2, color="gray")
        for i in overest_on_pulse:
            plt.axvline(x=i[0], ls=":", lw=2, color="red")
            plt.axvline(x=i[1], ls=":", lw=2, color="red")
        plt.errorbar(0.03*len(norm_prof), 0.9, yerr=off_pulse_noise_std,
                     color="gray", capsize=10, markersize=0, label="Error")
        plt.text(0.04*len(norm_prof), 0.89,
                 f"Noise: {round(off_pulse_noise_std, 5)}")

        # Add some custom stuff to the legend
        colors = ['gray', 'red']
        lines = [Line2D([0], [0], color=c, linewidth=3, linestyle=':')
                 for c in colors]
        lines.append(
            Line2D([0], [0], color="cornflowerblue", linewidth=3, linestyle="-"))
        lines.append(Line2D([0], [0], color="orange",
                     linewidth=3, linestyle="-"))
        lines.append(Line2D([0], [0], color="gray",
                     linewidth=3, linestyle="-"))
        labels = ["Maxima", "Noise estimate bounds",
                  "Profile", "Smoothed Profile", "Error"]
        plt.legend(lines, labels)
        plt.xlabel("Bins")
        plt.ylabel("Amplitude")
        plt.title(plot_name[:-4])
        plt.savefig(plot_name, bbox_inches="tight")
        plt.close()

    comp_est_dict = {}
    comp_est_dict["maxima"] = real_max
    comp_est_dict["underest_on_pulse"] = underest_on_pulse
    comp_est_dict["overest_on_pulse"] = overest_on_pulse
    comp_est_dict["underest_off_pulse"] = underest_off_pulse
    comp_est_dict["overest_off_pulse"] = overest_off_pulse
    comp_est_dict["norm_noise_std"] = off_pulse_noise_std
    comp_est_dict["norm_noise_mean"] = off_pulse_noise_mean
    comp_est_dict["alpha"] = alpha

    return comp_est_dict


def get_off_pulse(on_pulse_pairs):
    """Works out the off pulse ranges from the on-pulse.

    Parameters
    ----------
    on_pulse_pairs : `list`
        List of sub-lists/tuples. Each sub-list is a pair of numbers that
        describes the upper and lower bounds of the ON PULSE region.

    Returns
    -------
    off_pulse : `list`
        A list of two-lists that decribes the off-pulse region.
    """
    off_pulse = []
    # Make cycling generators to deal with wrap-around
    current_pair_gen = itertools.cycle(on_pulse_pairs)
    next_pair_gen = itertools.cycle(on_pulse_pairs)
    next(next_pair_gen)
    # Figure out the off-pulse profile
    for i in range(len(on_pulse_pairs)):
        current_pair = next(current_pair_gen)
        next_pair = next(next_pair_gen)
        off_pulse.append([current_pair[1], next_pair[0]])
    return off_pulse


def profile_region_from_pairs(profile, pairs):
    """Acquires a list of the region described by a pair of indexes accounting for wrap-around.

    Parameters
    ----------
    profile : `list`
        The profile that we want to get the region of.
    pairs : `list`/`tuple`
        A list of lists/tuples of length two that describes the indexes of the region(s) to grab.

    Returns
    -------
    region : `list`
        The region of the profile described by the pairs.
    """
    profile = list(profile)
    region = []
    for pair in pairs:
        if pair[0] < pair[1]:
            region = profile[pair[0]:pair[1]]
        else:
            region = profile[pair[0]:]
            region.extend(profile[:pair[1]])
    return region


def filled_profile_region_between_pairs(profile, pairs, fill_value=0):
    """Fills the off-pair region of a profile with the fill value.

    Parameters
    ----------
    profile : `list`
        The profile we want to fill.
    pairs : `list`/`tuple`
        A list of lists/tupls of length two that describes the region we DON'T want to fill.
    fill_value : `float`, optional
        The value to fill the region with. |br| Default=0.
    """
    filled_prof = np.array(profile)
    # Create two cycling generators to avoid out of bounds errors for the last pair
    current_pair_gen = itertools.cycle(pairs)
    next_pair_gen = itertools.cycle(pairs)
    next(next_pair_gen)
    for _ in pairs:
        # Next pair from the generator
        current_pair = next(current_pair_gen)
        next_pair = next(next_pair_gen)
        if current_pair[1] < next_pair[0]:
            filled_prof[current_pair[1]:next_pair[0]] = fill_value
        else:
            filled_prof[current_pair[1]:] = fill_value
            filled_prof[:next_pair[0]] = fill_value
    return filled_prof


def _fill_on_pulse_region(on_pulse_pairs, smoothed_profile, real_profile, noise_est):
    """Attempts to fill small gaps in the on_pulse_pairs of the smoothed
    profile provided these gaps resemble real signal.

    Parameters
    ----------
    on_pulse_pairs : `list`
        A list of pairs of indexes that mark the on-pulse regions of the profile.
    smoothed_profile : `list`
        The profile that has been smoothed with the stickel regularisation algorithm.
    real_profile : `list`
        The profile that stickel was applied to.
    noise_est : `float`
        A first order estimate of the profile's noise. This will determine what looks like real signal.

    Return:
    -------
    true_on_pulse : `list`
        A list of pairs that represent the index ranges of the ammended on-pulse regions.
        Note that the first or last pairs may have idx[0] > idx[1] which represents that
        the on-pulse region wraps over the end of the phase.
    """
    # Nothing to do if the pairs are of length 1
    if len(on_pulse_pairs) <= 1:
        return on_pulse_pairs

    real_profile = list(real_profile)
    loop_pairs = on_pulse_pairs.copy()
    # Create two cycling generators to avoid out of bounds errors for the last pair
    on_pulse_ammendments = []
    current_pair_gen = itertools.cycle(loop_pairs)
    next_pair_gen = itertools.cycle(loop_pairs)
    next(next_pair_gen)
    for i in range(len(loop_pairs)):
        # Next pair from the generator
        current_pair = next(current_pair_gen)
        next_pair = next(next_pair_gen)

        # Figure out off pulse region
        off_pulse_region = profile_region_from_pairs(
            real_profile, [[current_pair[1], next_pair[0]]])

        # See if this looks like signal
        if min(off_pulse_region) >= 3*noise_est:
            # Add to ammendments that we need to make
            on_pulse_ammendments.append([current_pair[1], next_pair[0]])

    # Apply the ammendments to the on-pulse region
    true_on_pulse = loop_pairs.copy()
    for ammendment_pair in on_pulse_ammendments:
        # The beginning point of the on-pulse to be ammended
        begin = ammendment_pair[0]
        end = ammendment_pair[1]
        # Find where the beginning meets up with the endpoints
        current_begs = np.array([i[0] for i in true_on_pulse])
        current_ends = np.array([i[1] for i in true_on_pulse])
        # Find the indexes of the start and end pairs
        idx_begin = np.where(current_ends == begin)[0][0]
        idx_end = np.where(current_begs == end)[0][0]
        # Modify the on-pulse pair
        start_of_begin = true_on_pulse[idx_begin][0]
        end_of_end = true_on_pulse[idx_end][1]
        true_on_pulse[idx_begin] = [start_of_begin, end_of_end]
        # Remove the other pair as it's been absorbed
        del true_on_pulse[idx_end]

    return true_on_pulse


def _profile_noise_initial_estimate(profile):
    """Attempts to guess the best value for alpha for :py:meth:`vcstools.prof_utils.sigmaClip`
    in order to get a good estimate of the off-pulse noise level for a profile.

    Parameters
    ----------
    profile : `list`
        The profile to clip.

    Returns
    -------
    best_alpha : `float`
        The chosen alpha value for :py:meth:`vcstools.prof_utils.sigmaClip`.
    best_nois : `float`
        The corresponding standard deviatino of the noise.
    """
    alphas = np.linspace(3.0, 2.0, 5)
    valid_alphas = []
    valid_noise = []
    for alpha in alphas:
        noise_std, off_pulse = sigmaClip(profile, alpha=alpha)
        n_off_pulse = np.count_nonzero(~np.isnan(off_pulse))
        # At least 10% of the profile is still unclipped
        if n_off_pulse >= 0.2*len(profile):
            valid_alphas.append(alpha)
            valid_noise.append(noise_std)
    best_noise, best_alpha = sorted(zip(valid_noise, valid_alphas))[
        0]  # lowest noise is prioritised
    return best_alpha, best_noise