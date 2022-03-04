"""Functions for preforming calculations related to the pulsar flux desnity radiometer equation
"""

import logging
import os
import sys
import numpy as np
import psrqpy
from matplotlib import pyplot as plt

#vcstools
from vcstools import data_load
from vcstools.config import load_config_file
from vcstools.beam_calc import get_beam_power_over_time, get_Trec,\
                               from_power_to_gain, source_beam_coverage,\
                               source_beam_coverage_and_times
from vcstools.metadb_utils import get_common_obs_metadata, obs_max_min,\
                                  mwa_alt_az_za
from vcstools.progress_bar import progress_bar
from vcstools.pulsar_spectra import flux_from_plaw, flux_from_spind,\
                                    find_spind, plot_flux_estimation
from vcstools.catalogue_utils import get_psrcat_ra_dec
from vcstools.gfit import gfit
from vcstools.beam_sim import read_sefd_file, launch_pabeam_sim
from vcstools.prof_utils import sn_calc
from vcstools import prof_utils

from pulsar_spectra.catalogues import flux_from_atnf

from mwa_pb import primarybeammap_tant as pbtant

logger = logging.getLogger(__name__)


# The two flux calculation methods

def flux_calc_radiometer_equation(profile, on_pulse_bool,
                                  noise_std, w_equiv_bins,
                                  sefd, u_sefd, t_int,
                                  bandwidth=30.72e6,
                                  sn_method="eqn_7.1"):
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
    sefd : `float`
        The System Equivalent Flux Density in Jy. Equivilent to system temperature divided by gain.
    u_sefd : `float`
        The uncertainty of the sefd.
    t_int : `float`
        The total observing time of the detection in seconds.
    bandwidth : `float`, optional
        The bandwidth of the observation in Hz. |br| Default: 30.72e6.
    sn_method : `str`, optional
        The method of calculating the signal to noise ratio out of ["eqn_7.1", "simple"]. Default "eqn_7.1".

        "Simple" uses 1/noise_std.

        "eqn_7.1" uses equation 7.1 of the pulsar handbook.

    Returns
    -------
    S_mean : `float`
        The mean flux density of the pulsar in mJy.
    u_S_mean : `float`
        The mean flux density uncertainty of the pulsar in mJy.
    """
    sn, u_sn = sn_calc(profile, on_pulse_bool, noise_std, w_equiv_bins, sn_method=sn_method)

    # Calculate S_mean
    S_mean = sn * sefd / (np.sqrt(2. * float(t_int) * bandwidth)) *\
             np.sqrt( w_equiv_bins / (len(profile) - w_equiv_bins)) * 1000.
    u_S_mean = np.sqrt( (u_sn / sn)**2 + (u_sefd / sefd)**2 ) * S_mean
    return S_mean, u_S_mean

def flux_calc_flux_profile(profile, noise_std,
                           sefd, u_sefd, t_int,
                           bandwidth=30.72e6):
    """Calculates the flux desnity using the siganl to noise (The maximum/peak
    signal divided by the STD of the noise) and the radiometer equation (see eqn
    A 1.21 of the pulsar handbook)

    Parameters
    ----------
    profile : `list`
        The normalised profile of the pulsar.
    noise_std : `float`
        The standard deviation of the noise in the normalised profile.
    sefd : `float`
        The System Equivalent Flux Density in Jy. Equivilent to system temperature divided by gain.
    u_sefd : `float`
        The uncertainty of the sefd.
    t_int : `float`
        The total observing time of the detection in seconds.
    bandwidth : `float`, optional
        The bandwidth of the observation in Hz. |br| Default: 30.72e6.

    Returns
    -------
    S_mean : `float`
        The mean flux density of the pulsar in mJy.
    u_S_mean : `float`
        The mean flux density uncertainty of the pulsar in mJy.
    """
    num_bins = len(profile)
    # Normalise to noise std
    profile = profile / noise_std
    # Work out how much observing time each bin represents (in sec)
    bin_time = t_int / num_bins

    # Put the profile into flux units with sefd
    flux_profile = profile * sefd / (np.sqrt(2. * bin_time * bandwidth))
    # Make an array of profile uncertainties. Profile uncertainty is assumed to be the std (1)
    u_flux_profile = np.sqrt( (1/profile)**2 + (u_sefd/sefd)**2 ) * flux_profile

    # Since we already took time into account we can just average and convert to mJy
    S_mean = np.mean(flux_profile) * 1000
    u_S_mean = np.sqrt(np.sum(u_flux_profile**2)) * 1000 / num_bins

    return S_mean, u_S_mean

def analyise_and_flux_cal(pulsar, bestprof_data,
                          flagged_tiles=None,
                          calid=None,
                          common_metadata=None,
                          trcvr=data_load.TRCVR_FILE,
                          simple_sefd=False, sefd_file=None,
                          vcstools_version='master',
                          args=None,
                          flux_method="radiometer"):
    """Analyise a pulse profile and calculates its flux density

    Parameters
    ----------
    pulsar : `str`
        The pulsar's Jname
    bestprof_data: list
        The output list from the function :py:meth:`vcstools.prof_utils.get_from_bestprof`

    Optional parameters:
    -------------------
    flagged_tiles : `str`
        The location of the flagged_tiles.txt file. If it's in the default location you can just supply the calid.
    calid : `int`
        The calibration ID of the detection. This is used to find the flagged_tiles.txt file.
    common_metadata: list
        The output of mwa_metadb_utils.get_common_obs_metadata(). If not supplied it will be downloaded.
    trcvr : `str`
        The file location of antena temperatures.
    simple_sefd : `boolean`
        If True perfom just a simple SEFD calculation instead of simulating the phased array beam response over the sky. Default: False.
    sefd_file : `str`
        The location of the pabeam.py's simulation of the phased array beam response over the sky output file. If not supplied will launch a pabeam.py simulation.
    vcstools_version : `str`
        The version of vcstools to use for the pabeam.py simulation. Default: master.
    args: Namespace
        The args from argparse to be used for job resubmission. Default: None.

    Returns
    -------
    det_kwargs: dict
    det_kwargs["flux"]: The mean flux density of the pulsar in mJy
    det_kwargs["flux_error"]: The flux desnity error in mJy
    det_kwargs["width"]: The equivalent width of the pulsar in ms
    det_kwargs["width_error"]: The error of the equivalent width in ms
    det_kwargs["scattering"]: The scattering width in s
    det_kwargs["scattering_error"]: The error of the scattering in s
    """
    # Load computer dependant config file
    comp_config = load_config_file()

    #unpack the bestprof_data
    obsid, prof_psr, _, period, _, sigma, beg, t_int, profile, num_bins, pointing = bestprof_data
    period=float(period)
    num_bins=int(num_bins)

    # Perform metadata calls
    if common_metadata is None:
        common_metadata = get_common_obs_metadata(obsid)
    obsid, ra, dec, dura, [xdelays, ydelays], centrefreq, channels = common_metadata
    # assume full bandwidth of 30.72 MHz
    bandwidth = 30.72e6

    # Find pulsar ra and dec
    _, pul_ra, pul_dec = get_psrcat_ra_dec(pulsar_list=[pulsar])[0]

    # Work out flagged tiles from calbration directory
    if not flagged_tiles:
        if calid:
            flagged_file = os.path.join(comp_config['base_data_dir'], obsid, "cal", calid, "rts", "flagged_tiles.txt")
            if os.path.exists(flagged_file):
                with open(flagged_file, "r") as ftf:
                    flagged_tiles = []
                    reader = csv.reader(ftf)
                    for row in reader:
                        flagged_tiles.append(row)
                    flagged_tiles = np.array(flagged_tiles).flatten()
            else:
                logger.warn("No flagged_tiles.txt file found so assuming no tiles have been flagged")
                flagged_tiles = []
        else:
            logger.warn("No flagged_tiles or calid provided so assuming no tiles have been flagged")
            flagged_tiles = []


    # Calc SEFD from the T_sys and gain
    if simple_sefd:
        t_sys, _, gain, u_gain = find_t_sys_gain(pulsar, obsid,
                                                 common_metadata=common_metadata,
                                                 beg=beg, end=(t_int + beg - 1))
        sefd = tsys / gain
    else:
        if sefd_file is None:
            launch_pabeam_sim(obsid, pointing, beg, t_int,
                              source_name=pulsar,
                              vcstools_version=vcstools_version,
                              flagged_tiles=flagged_tiles,
                              delays=xdelays,
                              args=args,
                              common_metadata=common_metadata)
            sys.exit(0)
        else:
            sefd_freq_time, sefd, u_sefd = read_sefd_file(sefd_file)

    #estimate S/N
    try:
        g_fitter = gfit(profile)
        g_fitter.plot_name = f"{obsid}_{pulsar}_{num_bins}_bins_gaussian_fit.png"
        g_fitter.component_plot_name = f"{obsid}_{pulsar}_{num_bins}_bins_gaussian_components.png"
        g_fitter.auto_fit()
        g_fitter.plot_fit()
        prof_dict = g_fitter.fit_dict
    except (prof_utils.ProfileLengthError, prof_utils.NoFitError):
        logger.info("Profile couldn't be fit. Using old style of profile analysis")
        prof_dict = prof_utils.auto_analyse_pulse_prof(profile, period)

    if not prof_dict:
        logger.warn("Profile could not be analysed using any methods")
        det_kwargs = {}
        det_kwargs["flux"]              = None
        det_kwargs["flux_error"]        = None
        det_kwargs["width"]             = None
        det_kwargs["width_error"]       = None
        det_kwargs["scattering"]        = None
        det_kwargs["scattering_error"]  = None
        return det_kwargs, None, None

    # Unpack dictionary
    sn = prof_dict["sn"]
    u_sn = prof_dict["sn_e"]
    profile = prof_dict["profile"]
    on_pulse_bool = prof_dict["on_pulse_bool"]
    noise_std = prof_dict["noise_std"]
    noise_mean = prof_dict["noise_mean"]
    w_equiv_phase = prof_dict["Weq"]
    u_w_equiv_phase =  prof_dict["Weq_e"]
    w_equiv_bins = w_equiv_phase * num_bins
    u_w_equiv_bins = w_equiv_phase * num_bins
    w_equiv_ms = period * w_equiv_phase
    u_w_equiv_ms = period * u_w_equiv_phase
    scattering = prof_dict["Wscat"]*period/1000 #convert to seconds
    u_scattering = prof_dict["Wscat_e"]*period/1000
    scattered = prof_dict["scattered"]

    logger.info("Profile scattered? {0}".format(scattered))
    logger.info("S/N: {0} +/- {1}".format(sn, u_sn))
    #logger.debug("Gain {0} K/Jy".format(gain))
    logger.debug("Equivalent width in ms: {0}".format(w_equiv_ms))
    #logger.debug("T_sys: {0} K".format(t_sys))
    logger.debug("Detection time: {0}".format(t_int))
    logger.debug("Number of bins: {0}".format(num_bins))


    # Renormalise around the noise mean
    noise = []
    noise_i = []
    on_pulse = []
    on_pulse_i = []
    for p, b, i in zip(profile, on_pulse_bool, range(len(profile))):
        if not b:
            noise.append(p)
            noise_i.append(i)
        else:
            on_pulse.append(p)
            on_pulse_i.append(i)
    plt.scatter(on_pulse_i, on_pulse, color="blue")
    plt.scatter(noise_i, noise, color="red")
    plt.savefig("test.png")
    noise_mean = np.mean(noise)
    print(f"Noise mean: {noise_mean}")
    profile = (profile - noise_mean) / max(profile - noise_mean)
    logger.debug(list(profile))
    logger.debug(on_pulse_bool)
    logger.debug(noise_std, w_equiv_bins, sefd, u_sefd, t_int)

    # Final calc of the mean flux density in mJy
    if flux_method == "radiometer":
        S_mean, u_S_mean = flux_calc_radiometer_equation(profile, on_pulse_bool,
                            noise_std, w_equiv_bins,
                            sefd, u_sefd, t_int, bandwidth=bandwidth)
    elif flux_method == "flux_profile":
        S_mean, u_S_mean = flux_calc_flux_profile(profile,
                            noise_std,
                            sefd, u_sefd, t_int, bandwidth=bandwidth)


    logger.info('Smean {0:.3f} +/- {1:.3f} mJy'.format(S_mean, u_S_mean))

    #prevent TypeError caused by trying to format Nones given to fluxes for highly scattered pulsars
    S_mean = float("{0:.3f}".format(S_mean))
    u_S_mean = float("{0:.3f}".format(u_S_mean))

    # Plot flux comparisons for ANTF
    freq_all, flux_all, flux_err_all, _ = flux_from_atnf(pulsar)
    logger.debug("Freqs: {0}".format(freq_all))
    logger.debug("Fluxes: {0}".format(flux_all))
    logger.debug("Flux Errors: {0}".format(flux_err_all))
    logger.debug("{0} there are {1} flux values available on the ATNF database"\
                .format(pulsar, len(flux_all)))

    # Check if there is enough data to estimate the flux
    #if len(flux_all) == 0:
    #    logger.debug("{} no flux values on archive. Cannot estimate flux.".format(pulsar))
    #elif ( len(flux_all) == 1 ) and ( ( not spind ) or ( not spind_err ) ):
    #    logger.debug("{} has only a single flux value and no spectral index. Cannot estimate flux. Will return Nones".format(pulsar))
    #else:
    #    spind, spind_err, K, covar_mat = find_spind(pulsar, freq_all, flux_all, flux_err_all)
    #    plot_flux_estimation(pulsar, freq_all, flux_all, flux_err_all, spind,
    #                            my_nu=centrefreq, my_S=S_mean, my_S_e=u_S_mean, obsid=obsid,
    #                            a_err=spind_err,  K=K, covar_mat=covar_mat)

    #format data for uploading
    w_equiv_ms   = float("{0:.2f}".format(w_equiv_ms))
    u_w_equiv_ms = float("{0:.2f}".format(u_w_equiv_ms))
    scattering   = float("{0:.5f}".format(scattering))
    u_scattering = float("{0:.5f}".format(u_scattering))

    det_kwargs = {}
    det_kwargs["flux"]              = S_mean
    det_kwargs["flux_error"]        = u_S_mean
    det_kwargs["width"]             = w_equiv_ms
    det_kwargs["width_error"]       = u_w_equiv_ms
    det_kwargs["scattering"]        = scattering
    det_kwargs["scattering_error"]  = u_scattering
    return det_kwargs, sn, u_sn


def est_pulsar_flux(pulsar, obsid, plot_flux=False, common_metadata=None, query=None):
    """Estimates a pulsar's flux from archival data by assuming a power law relation between flux and frequency

    Parameters
    ----------
    pulsar : `str`
        The Jname of the pulsar.
    obsid : `int`
        The MWA observation ID
    plot_flux : `boolean`, optional
        Whether or not to make a plot of the flux estimation. |br| Default: `False`
    common_metadata : `list`, optional
        The list of common metadata generated from :py:meth:`vcstools.metadb_utils.get_common_obs_metadata`
    query : psrqpy object, optional
        A previous psrqpy.QueryATNF query. Can be supplied to prevent performing a new query.

    Returns
    -------
    flux : `float`
        The estimated flux in Jy.
    flux_err : `float`
        The estimated flux's uncertainty in Jy.
    """
    if query is None:
        query = psrqpy.QueryATNF(psrs=[pulsar], loadfromdb=data_load.ATNF_LOC).pandas
    query_id = list(query['PSRJ']).index(pulsar)
    if common_metadata is None:
        logger.debug("Obtaining mean freq from obs metadata")
        common_metadata = get_common_obs_metadata(obsid)
    f_mean = common_metadata[5]*1e6

    #freq_all, flux_all, flux_err_all, spind, spind_err = flux_from_atnf(pulsar, query=query)
    freq_all, flux_all, flux_err_all, ref_all = flux_from_atnf(pulsar, query=query)
    # convert to Hz
    freq_all     = [f*1e6 for f in freq_all]
    # convert to Jy
    flux_all     = [f*1e-3 for f in flux_all]
    flux_err_all = [f*1e-3 for f in flux_err_all]

    spind = query["SPINDX"][query_id]
    spind_err = query["SPINDX_ERR"][query_id]

    logger.debug("Freqs: {0}".format(freq_all))
    logger.debug("Fluxes: {0}".format(flux_all))
    logger.debug("Flux Errors: {0}".format(flux_err_all))
    logger.debug("{0} there are {1} flux values available on the ATNF database"\
                .format(pulsar, len(flux_all)))

    # Check if there is enough data to estimate the flux
    if len(flux_all) == 0:
        logger.debug("{} no flux values on archive. Cannot estimate flux. Will return Nones".format(pulsar))
        return None, None
    elif ( len(flux_all) == 1 ) and ( ( not spind ) or ( not spind_err ) ):
        logger.debug("{} has only a single flux value and no spectral index. Cannot estimate flux. Will return Nones".format(pulsar))
        return None, None

    if ( not spind ) or ( not spind_err ):
        # If no spind on ATNF fit our own
        spind, spind_err, K, covar_mat = find_spind(pulsar, freq_all, flux_all, flux_err_all)
    else:
        K = covar_mat = None

    if K and covar_mat is not None:
        # Use the spind power law we fit
        flux_est, flux_est_err = flux_from_plaw(f_mean, K, spind, covar_mat)
    elif spind and spind_err and flux_all:
        # Use ATNF spind
        flux_est, flux_est_err = flux_from_spind(f_mean, freq_all[0], flux_all[0], flux_err_all[0],\
                                                 spind, spind_err)
    logger.debug("Finished estimating flux")

    if plot_flux == True:
        plot_flux_estimation(pulsar, freq_all, flux_all, flux_err_all, spind,
                             my_nu=f_mean, my_S=flux_est, my_S_e=flux_est_err, obsid=obsid,
                             a_err=spind_err,  K=K, covar_mat=covar_mat)

    return flux_est, flux_est_err


def find_pulsar_w50(pulsar, query=None):
    """Attempts to find a pulsar's W50 from the ATNF catalogue. If unavailable, will estimate

    Parameters
    ----------
    pulsar : `str`
        The Jname of the pulsar.
    query : psrqpy object, optional
        A previous psrqpy.QueryATNF query. Can be supplied to prevent performing a new query.

    Returns
    -------
    w50 : `float`
        The pulsar's w50 in s.
    w50_err : `float`
        The W50's uncertainty.
    """
    if query is None:
        #returns W_50 and error for a pulsar from the ATNF archive IN SECONDS
        logger.debug("Accessing ATNF database")
        query = psrqpy.QueryATNF(psrs=[pulsar], loadfromdb=data_load.ATNF_LOC).pandas
    query_id = list(query['PSRJ']).index(pulsar)

    W_50     = query["W50"][query_id]
    W_50_err = query["W50_ERR"][query_id]
    if np.isnan(W_50):
        W_50 = np.nan
        W_50_err = np.nan
    else:
        #convert to seconds
        W_50 = W_50 / 1000.

    if np.isnan(W_50_err) and not np.isnan(W_50):
        logger.debug("{} W_50 error not on archive. returning standard 5% error".format(pulsar))
        W_50_err = W_50 * 0.05
    else:
        #convert to seconds
        W_50_err = W_50_err / 1000.

    if np.isnan(W_50):
        logger.debug("{} applying estimated W_50. Uncertainty will be inflated".format(pulsar))
        #Rankin1993 - W = x*P^0.5 where x=4.8+/-0.5 degrees of rotation at 1GHz
        #We will nflate this error due to differing frequencies and pulsar behaviour. W_50_err=1. degrees
        coeff = 4.8
        coeff_err = 2.
        period = float(query["P0"][query_id])

        #This estimation is worse for msps, add extra uncetainty if period < 50ms
        if period <0.05:
            coeff_err = 4.

        #calculate
        W_50 = coeff*period**0.5
        W_50_err = coeff_err*period**0.5 #error from variance calculation

        #W_50 is now in degrees of rotation phase. Convert to seconds
        W_50 = (W_50/360.)*period
        W_50_err = (W_50_err/360.)*period

    if np.isnan(W_50):
        W_50=None
    if np.isnan(W_50_err):
        W_50_err=None
    return W_50, W_50_err


def find_t_sys_gain(pulsar, obsid,
                    p_ra=None, p_dec=None,
                    dect_beg=None, dect_end=None,
                    obs_beg=None, obs_end=None,
                    common_metadata=None, full_metadata=None,
                    query=None, min_z_power=0.3, trcvr=data_load.TRCVR_FILE):
    """Finds the system temperature and gain for an observation.

    Parameters
    ----------
    pulsar : `str`
        The Jname of the pulsar.
    obsid : `int`
        The MWA Observation ID.
    p_ra, p_dec : `str`, optional
        The target's right ascension and declination in sexidecimals. If not supplied will use the values from the ANTF.
    dect_beg, dect_end : `int`, optional
        The beg and end GPS time of the detection to calculate over.
        If not supplied will estimate beam enter and exit.
    obs_beg, obs_end : `int`, optional
        Beginning and end GPS time of the observation.
        If not supplied will use :py:meth:`vcstools.metadb_utils.obs_max_min` to find it.
    common_metadata : `list`, optional
        The list of common metadata generated from :py:meth:`vcstools.metadb_utils.get_common_obs_metadata`.
    full_metadata : `dict`, optional
        The dictionary of metadata generated from :py:meth:`vcstools.metadb_utils.getmeta`.
    query : psrqpy object, optional
        A previous psrqpy.QueryATNF query. Can be supplied to prevent performing a new query.
    min_z_power : `float`, optional
        Zenith normalised power cut off. |br| Default: 0.3.
    trcvr : `str`
        The location of the MWA receiver temp csv file. |br| Default: <vcstools_data_dir>MWA_Trcvr_tile_56.csv

    Returns
    -------
    t_sys : `float`
        The system temperature in K.
    t_sys_err : `float`
        The system temperature's uncertainty.
    gain : `float`
        The system gain in K/Jy.
    gain_err : `float`
        The gain's uncertainty.
    """
    # get ra and dec if not supplied
    if query is None:
        logger.debug("Obtaining pulsar RA and Dec from ATNF")
        query = psrqpy.QueryATNF(psrs=[pulsar], loadfromdb=data_load.ATNF_LOC).pandas
    query_id = list(query['PSRJ']).index(pulsar)
    if not p_ra or not p_dec:
        p_ra = query["RAJ"][query_id]
        p_dec= query["DECJ"][query_id]

    # get metadata if not supplied
    if not common_metadata:
        logger.debug("Obtaining obs metadata")
        common_metadata = get_common_obs_metadata(obsid)

    obsid, obs_ra, obs_dec, _, delays, centrefreq, channels = common_metadata

    if not dect_beg or not dect_end:
        # Estimate integration time from when the source enters and exits the beam
        dect_beg, dect_end = source_beam_coverage_and_times(obsid, pulsar,
                                p_ra=p_ra, p_dec=p_dec,
                                obs_beg=obs_beg, obs_end=obs_end,
                                min_z_power=min_z_power,
                                common_metadata=common_metadata,
                                query=query)[:2]
    start_time = dect_end - int(obsid)
    t_int = dect_end - dect_beg + 1

    #Get important info
    ntiles = 128 #TODO actually we excluded some tiles during beamforming, so we'll need to account for that here

    beam_power = get_beam_power_over_time(np.array([[pulsar, p_ra, p_dec]]),
                                          common_metadata=[obsid, obs_ra, obs_dec, t_int, delays,
                                                           centrefreq, channels],
                                          dt=100, start_time=start_time)
    mean_beam_power = np.mean(beam_power)

    # Usa a primary beam function to convolve the sky temperature with the primary beam
    # prints suppressed
    sys.stdout = open(os.devnull, 'w')
    _, _, Tsky_XX, _, _, _, Tsky_YY, _ = pbtant.make_primarybeammap(int(obsid), delays, centrefreq*1e6, 'analytic', plottype='None')
    sys.stdout = sys.__stdout__

    #TODO can be inaccurate for coherent but is too difficult to simulate
    t_sky = (Tsky_XX + Tsky_YY) / 2.
    # Get T_sys by adding Trec and Tsky (other temperatures are assumed to be negligible
    t_sys_table = t_sky + get_Trec(centrefreq, trcvr_file=trcvr)
    t_sys = np.mean(t_sys_table)
    t_sys_err = t_sys*0.02 #TODO: figure out what t_sys error is

    logger.debug("pul_ra: {} pul_dec: {}".format(p_ra, p_dec))
    _, _, zas = mwa_alt_az_za(obsid, ra=p_ra, dec=p_dec)
    theta = np.radians(zas)
    gain = from_power_to_gain(mean_beam_power, centrefreq*1e6, ntiles, coh=True)
    logger.debug("mean_beam_power: {} theta: {} pi: {}".format(mean_beam_power, theta, np.pi))
    gain_err = gain * ((1. - mean_beam_power)*0.12 + 2.*(theta/(0.5*np.pi))**2. + 0.1)

    return t_sys, t_sys_err, gain, gain_err


def est_pulsar_sn(pulsar, obsid,
                  p_ra=None, p_dec=None,
                  dect_beg=None, dect_end=None,
                  obs_beg=None, obs_end=None,
                  common_metadata=None, full_metadata=None,
                  query=None, plot_flux=False,
                  min_z_power=0.3, trcvr=data_load.TRCVR_FILE):
    """Estimates the signal to noise ratio for a pulsar in a given observation using the radiometer equation

    .. math:: S/N = \frac{s_{mean}  gain}{t_{sys}}  \sqrt{n_p  t_{int} df \frac{period - W_{50}}{W_{50}}

    Note that W_50 should be W_equiv but we can't figure that out so we're estimating

    Parameters
    ----------
    pulsar : `str`
        The Jname of the pulsar.
    obsid : `int`
        The MWA Observation ID.
    p_ra, p_dec : `str`, optional
        The target's right ascension and declination in sexidecimals. If not supplied will use the values from the ANTF.
    dect_beg, dect_end : `int`, optional
        The beg and end GPS time of the detection to calculate over.
        If not supplied will estimate beam enter and exit.
    obs_beg, obs_end : `int`, optional
        Beginning and end GPS time of the observation.
        If not supplied will use :py:meth:`vcstools.metadb_utils.obs_max_min` to find it.
    common_metadata : `list`, optional
        The list of common metadata generated from :py:meth:`vcstools.metadb_utils.get_common_obs_metadata`
    full_metadata : `dict`, optional
        The dictionary of metadata generated from :py:meth:`vcstools.metadb_utils.getmeta`
    query : psrqpy object, optional
        A previous psrqpy.QueryATNF query. Can be supplied to prevent performing a new query.
    plot_flux : `boolean`
        OPTIONAL - whether or not to produce a plot of the flux estimation. |br| Default: `False`.
    min_z_power : `float`, optional
        Zenith normalised power cut off. |br| Default: 0.3.
    trcvr : `str`
        The location of the MWA receiver temp csv file. |br| Default: <vcstools_data_dir>/MWA_Trcvr_tile_56.csv

    Returns
    -------
    sn : `float`
        The expected signal to noise ratio for the given inputs.
    sn_err : `float`
        The uncertainty in the signal to noise ratio.
    """
    # We will attain uncertainties for s_mean, gain, t_sys and W_50.
    # other uncertainties are considered negligible
    if query is None:
        query = psrqpy.QueryATNF(psrs=pulsar, loadfromdb=data_load.ATNF_LOC).pandas
    query_id = list(query['PSRJ']).index(pulsar)

    if p_ra is None or p_dec is None:
        # Get some basic pulsar and obs info info
        p_ra = query["RAJ"][query_id]
        p_dec = query["DECJ"][query_id]

    # Get metadata if not supplied
    if not common_metadata or not full_metadata:
        logger.debug("Obtaining obs metadata")
        common_metadata, full_metadata = get_common_obs_metadata(obsid, return_all=True, full_metadata=full_metadata)

    n_p = 2 #constant
    df = 30.72e6 #(24*1.28e6)

    # Estimate flux
    s_mean, s_mean_err = est_pulsar_flux(pulsar, obsid, plot_flux=plot_flux,
                                         common_metadata=common_metadata, query=query)
    # Fluxes may be Nones. If so, return None
    if s_mean is None and s_mean_err is None:
        return None, None, None, None

    if not dect_beg or not dect_end:
        # Estimate integration time from when the source enters and exits the beam
        dect_beg, dect_end = source_beam_coverage_and_times(obsid, pulsar,
                                p_ra=p_ra, p_dec=p_dec,
                                obs_beg=obs_beg, obs_end=obs_end,
                                min_z_power=min_z_power,
                                common_metadata=common_metadata,
                                query=query)[:2]
    t_int = dect_end - dect_beg + 1
    if t_int<=0.:
        logger.warning("{} not in beam for obs files or specificed beginning and end times"\
                    .format(pulsar))
        return None, None, None, None

    # Find system temp and gain
    t_sys, t_sys_err, gain, gain_err = find_t_sys_gain(pulsar, obsid,
                p_ra=p_ra, p_dec=p_dec,
                dect_beg=dect_beg, dect_end=dect_end,
                obs_beg=obs_beg, obs_end=obs_end,
                common_metadata=common_metadata, full_metadata=full_metadata,
                trcvr=trcvr, min_z_power=min_z_power, query=query)

    #Find W_50
    W_50, W_50_err = find_pulsar_w50(pulsar, query=query)

    # Calculate SN
    period = float(query["P0"][query_id])
    SN = ((s_mean * gain)/t_sys) * np.sqrt(n_p * t_int * df * (period - W_50)/W_50)

    #Calculate SN uncertainty using variance formula. Assuming error from period, df and t_int is zero
    dc_expr    = np.sqrt((period-W_50)/W_50)
    var_s_mean = (gain * np.sqrt(n_p * t_int * df)) * dc_expr / t_sys
    var_gain   = s_mean * np.sqrt(n_p * t_int * df) * dc_expr / t_sys
    var_W_50   = s_mean * gain * np.sqrt(n_p * t_int * df)/t_sys * (period/(-2.*W_50**2.)) * dc_expr**-1
    var_t_sys  = -s_mean * gain * np.sqrt(n_p * t_int * df) * dc_expr / t_sys**2.
    var_s_mean = var_s_mean**2. * s_mean_err**2.
    var_gain   = var_gain**2.   * gain_err**2.
    var_W_50   = var_W_50**2.   * W_50_err**2.
    var_t_sys  = var_t_sys**2.  * t_sys_err**2.

    logger.debug("variance estimates for s_mean: {0}, gain: {1}, W_50: {2}, t_sys: {3}"\
                 .format(var_s_mean, var_gain, var_W_50, var_t_sys))
    SN_err = np.sqrt(var_s_mean + var_gain + var_W_50 + var_t_sys)

    logger.debug("Gain: {0} +/- {1}".format(gain, gain_err))
    logger.debug("t_int: {0}".format(t_int))
    logger.debug("df: {0}".format(df))
    logger.debug("period: {0}".format(period))
    logger.debug("W_50: {0} +/- {1}".format(W_50, W_50_err))
    logger.debug("t_sys: {0} +/- {1}".format(t_sys, t_sys_err))

    return SN, SN_err, s_mean, s_mean_err


def multi_psr_snfe(pulsar_list, obsid, obs_beg, obs_end,
                   common_metadata=None, full_metadata=None,
                   query=None, plot_flux=False,
                   min_z_power=0.3, trcvr=data_load.TRCVR_FILE):
    """Runs :py:meth:`vcstools.sn_flux_utils.est_pulsar_sn` for multiple pulsars in the same MWA observation.

    Parameters
    ----------
    pulsar : `list`
        A list of the pulsar Jnames.
    obsid : `int`
        The MWA Observation ID.
    obs_beg, obs_end : `int`
        Beginning and end GPS time of the observation.
        If not supplied will use :py:meth:`vcstools.metadb_utils.obs_max_min` to find it.
    common_metadata : `list`, optional
        The list of common metadata generated from :py:meth:`vcstools.metadb_utils.get_common_obs_metadata`
    full_metadata : `dict`, optional
        The dictionary of metadata generated from :py:meth:`vcstools.metadb_utils.getmeta`
    query : psrqpy object, optional
        A previous psrqpy.QueryATNF query. Can be supplied to prevent performing a new query.
    plot_flux : `boolean`, optional
        If `True` will produce a plot of the flux estimation. |br| Default: False
    min_z_power : `float`, optional
        Zenith normalised power cut off. |br| Default: 0.3.
    trcvr : `str`
        The location of the MWA receiver temp csv file. |br| Default: <vcstools_data_dir>MWA_Trcvr_tile_56.csv

    Returns
    -------
    sn_dict : `dict`
        A dictionary where eacy key is the pulsar jname and contains a list of the following
        sn_dict[pulsar]=[sn, sn_e, s, s_e]
        sn : `float`
            The expected signal to noise ratio for the given inputs.
        sn_err : `float`
            The uncertainty in the signal to noise ratio.
        s : `float`
            The expected flux density of the pulsar.
        s_e : `float`
            The uncertainty expected flux density of the pulsar.
    """
    logger.info("""This script may use estimations where data is missing.
    For full verbosity, use the DEBUG logger (ie. -L DEBUG)""")

    if common_metadata is None or full_metadata is None:
        logger.debug("Obtaining obs metadata")
        common_metadata, full_metadata = get_common_obs_metadata(obsid, return_all=True, full_metadata=full_metadata)

    if obs_beg is None or obs_end is None:
        obs_beg, obs_end = obs_max_min(obsid)

    mega_query = psrqpy.QueryATNF(psrs=pulsar_list, loadfromdb=data_load.ATNF_LOC).pandas
    sn_dict = {}
    for i, pulsar in enumerate(progress_bar(mega_query["PSRJ"], "Calculating pulsar SN: ")):
        psr_query = {}
        for key in mega_query.keys():
            psr_query[key] = [mega_query[key][i]]

        sn, sn_e, s, s_e = est_pulsar_sn(pulsar, obsid,
                                         obs_beg=obs_beg, obs_end=obs_end,
                                         common_metadata=common_metadata, full_metadata=full_metadata,
                                         plot_flux=plot_flux, query=psr_query,
                                         min_z_power=min_z_power, trcvr=trcvr)

        sn_dict[pulsar]=[sn, sn_e, s, s_e]

    return sn_dict