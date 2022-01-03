"""Functions for preforming calculations related to the pulsar flux desnity radiometer equation
"""

import logging
import os
import sys
import numpy as np
import psrqpy

#vcstools
from vcstools import data_load
from vcstools.beam_calc import get_beam_power_over_time, get_Trec,\
                               from_power_to_gain, source_beam_coverage,\
                               source_beam_coverage_and_times
from vcstools.metadb_utils import get_common_obs_metadata, obs_max_min,\
                                  mwa_alt_az_za
from vcstools.progress_bar import progress_bar
from vcstools.pulsar_spectra import flux_from_plaw, flux_from_spind,\
                                    find_spind, plot_flux_estimation

from pulsar_spectra.catalogues import flux_from_atnf

from mwa_pb import primarybeammap_tant as pbtant

logger = logging.getLogger(__name__)


def est_pulsar_flux(pulsar, obsid, plot_flux=False, common_metadata=None, query=None):
    """
    Estimates a pulsar's flux from archival data by assuming a power law relation between flux and frequency

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