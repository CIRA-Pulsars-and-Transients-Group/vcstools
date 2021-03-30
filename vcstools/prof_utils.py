import numpy as np
import subprocess
import re
import logging
from astropy.time import Time

from vcstools.config import load_config_file

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
    f.close()
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
        logger.debug(f"max prof: {max(prof_data)} +/- {prof_e} ")
        logger.debug(f"sigma   : {sigma} +/- {sigma_e} ")
        logger.debug(f"sn      : {sn} +/- {sn_e} ")

    return [sn, sn_e, scattered]


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
            logger.info(f"Skipping alpha value: {alpha}")
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

    logger.info(f"Best profile analysis using an alpha value of {best_alpha}")

    return fit_dict


# SigmaClip isn't perfect. Use these next function to check for bad clips
def check_clip(prof_to_check, toomuch=0.8, toolittle_frac=0., toolittle_absolute=4):
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
    for i in prof_to_check:
        if np.isnan(i):
            num_nans += 1
    if num_nans <= toolittle_frac*len(prof_to_check) or num_nans <= toolittle_absolute:
        raise LittleClipError("Not enough data has been clipped. Condsier trying a smaller alpha value when clipping.")
    elif num_nans >= toomuch*len(prof_to_check):
        raise LargeClipError("A large portion of the data has been clipped. Condsier trying a larger alpha value when clipping.")
