#other
import logging
import os
import sys
import numpy as np
import psrqpy
import glob

#matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#Astropy
from astropy.table import Table

#vcstools
from vcstools import data_load
from vcstools.catalogue_utils import grab_source_alog
from vcstools.beam_calc import get_beam_power_over_time, find_sources_in_obs,\
                               from_power_to_gain, get_Trec
from vcstools.metadb_utils import get_common_obs_metadata, obs_max_min,\
                                  mwa_alt_az_za, get_channels
from vcstools.progress_bar import progress_bar

from mwa_pb import primarybeammap_tant as pbtant

logger = logging.getLogger(__name__)


def find_combined_beg_end(obsid, base_path="/astro/mwavcs/vcs/", channels=None):
    """
    looks through the comined files of the input obsid and returns the max and min in gps time

    Input:
        obsid: The MWA observation ID
    Optional Input:
        base_path: the direct path the base vcs directory. Default: /astro/mwavcs/vcs/\
        channels: a list of the frequency channel ids. Default None which then gets the
                  from the mwa metadata
    """
    #TODO have some sort of check to look for gaps
    if glob.glob("{0}/{1}/combined/{1}*_ics.dat".format(base_path, obsid)):
        combined_files = glob.glob("{0}/{1}/combined/{1}*_ics.dat".format(base_path, obsid))
    else:
        channels = get_channels(obsid, channels)
        combined_files = glob.glob("{0}/{1}/combined/{1}*_ch{2}.dat".\
                                   format(base_path, obsid, channels[-1]))
    if len(combined_files) > 0:
        comb_times = []
        for comb in combined_files:
            comb_times.append(int(comb.split("_")[1]))
        beg = min(comb_times)
        end = max(comb_times)
    else:
        logger.warn("No combined files on disk for {0}".format(obsid))
        beg = None
        end = None

    return beg, end


def plot_flux_estimation(pulsar, nu_atnf, S_atnf, S_atnf_e, a,\
                        my_nu=None, my_S=None, my_S_e=None, obsid=None,\
                        K=None, covar_mat=None, a_err=None):
    """
    Used for plotting the estimated flux density against the known flux values.
    Can plot against either a least-sqaures fit flux or a spectral-index calculated flux.
    For the former, supply the covariance matrix and the K value.
    For the latter supply the spectral index error.

    Parameters:
    -----------
    pulsar: string
        The name of the pulsar
    nu_atnf: list
        The frequencies in which the known flux values correspond to (Hz)
    S_atnf: list
        The known flux values (Jy)
    S_atnf_e: list
        The uncertainties correspodning to the known fluxes
    a: float
        The spectral index
    my_nu: float
        OPTIONAL - The frequency you're estimating the flux at (Hz). Default: None
    my_S: float
        OPTIONAL - The estimated flux (Jy). Default: None
    my_S_e: float
        OPTIONAL - The uncertainty in the estimated flux (Jy). Default: None
    obsid: int
        OPTIONAL - The observation ID. Default: None
    K: float
        OPTIONAL - The K value of the least-squares fit. Use only when a least-sqaures fit has been done. Default: None
    covar_mat: numpy.matrix object
        OPTIONAL - The covariance matrix from the least-squares fit. Use only when a least-squares fit has been done. Default: None
    a_err: float
        OPTIONAL - The error in the spectral index. Use only when the flux has been estimated without least-squares. Default: None
    """
    #making title and .png name
    title ="{0} Flux Estimate".format(pulsar)
    save_name = "flux_density_{0}".format(pulsar)
    if obsid:
        title += " {}".format(obsid)
        save_name += "_{}".format(obsid)
    else:
        title += " ATNF"
        save_name += "_ATNF"
    save_name += ".png"

    #input data
    nu_range = list(nu_atnf)
    if my_nu:
        nu_range.append(my_nu)
    S_range = list(S_atnf)
    if my_S:
        S_range.append(my_S)

    #making errors
    nu_cont = np.logspace(np.log10(min(nu_range)), np.log10(max(nu_range)), num=500)
    if covar_mat is not None and K is not None:
        S_cont, S_cont_e = flux_from_plaw(nu_cont, K, a, covar_mat)
        a_err = np.sqrt(abs(covar_mat.item(3)))
    elif a_err is not None:
        S_cont, S_cont_e = flux_from_spind(nu_cont, nu_atnf[0], S_atnf[0], S_atnf_e[0], a, a_err)
    else:
        logger.warn("Requires more information to plot. Please refer to docs for more info.")
        return

    #x ticks
    possible_ticks = [1e7, 3e7, 1e8, 3e8, 1e9, 3e9, 1e10, 3e10, 1e11]
    lower = min(possible_ticks, key=lambda x:abs(x-min(nu_range)))
    upper = min(possible_ticks, key=lambda x:abs(x-max(nu_range)))
    lower_idx = possible_ticks.index(lower)
    upper_idx = possible_ticks.index(upper)
    if min(nu_range) < possible_ticks[lower_idx]:
        lower_idx = lower_idx - 1
    if max(nu_range) > possible_ticks[upper_idx]:
        upper_idx = upper_idx + 1
    xticks = []
    xtick_labels = []
    for i in range(lower_idx, (upper_idx + 1)):
        xticks.append(possible_ticks[i])
        xtick_labels.append(possible_ticks[i]/1e9)

    #y ticks
    possible_ticks = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1e0, 3e0, 1e1, 3e1, 1e2, 3e2, 1e3]
    lower = min(possible_ticks, key=lambda x:abs(x-min(S_range)))
    upper = min(possible_ticks, key=lambda x:abs(x-max(S_range)))
    lower_idx = possible_ticks.index(lower)
    upper_idx = possible_ticks.index(upper)
    if min(S_range) < possible_ticks[lower_idx]:
        lower_idx = lower_idx - 1
    if max(S_range) > possible_ticks[upper_idx]:
        upper_idx = upper_idx + 1
    yticks = []
    for i in range(lower_idx, (upper_idx + 1)):
        yticks.append(possible_ticks[i])

    #Plotting
    _, ax = plt.subplots(figsize=(12, 8))
    ax.grid()
    plt.text(0.05, 0.1, "Derived spectral index = {0} +/- {1}".format(round(a, 2), round(a_err, 4)),\
            fontsize=10, color="black", transform=ax.transAxes)
    if my_nu and my_S and my_S_e:
        plt.text(0.05, 0.05, "Flux est at {0}MHz: {1} +/- {2}mJy".format(round(my_nu/1e6, 2), round(my_S*1000, 2), round(my_S_e*1000, 2)),\
                fontsize=10, color="black", transform=ax.transAxes)
        plt.errorbar(my_nu, my_S, yerr=my_S_e, fmt='o', label="Extrapolated data points", color="orange") # Extrapolated data point

    plt.fill_between(nu_cont, S_cont - S_cont_e, S_cont + S_cont_e, facecolor='gray') # Model errors
    plt.plot(nu_cont, S_cont, 'k--', label="model") # Modelled line
    plt.errorbar(nu_atnf, S_atnf, yerr=S_atnf_e, fmt='o', label="ATNF data points", color="blue") # ATNF data points
    plt.yscale('log')
    plt.xscale('log')
    plt.xticks(xticks, xtick_labels)
    plt.yticks(yticks, yticks)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Flux (Jy)')
    plt.title(title)
    plt.legend()
    plt.savefig(save_name)
    plt.close()


def pulsar_beam_coverage(obsid, pulsar, beg, end, metadata=None, full_meta=None, ondisk=False, min_z_power=0.3, query=None):
    """
    Finds the normalised time that a pulsar is in the beam for a given obsid
    If pulsar is not in beam, returns None, None

    Parameters:
    -----------
    obsid: int
        The observation ID
    pulsar: string
        The pulsar's J name
    beg: int
        The beginning of the processed observing time in gps time
    end: int
        The end of the processed observing time in gps time
    obs_beg: int
        OPTIONAL - The beginning of the observation in gps time
    obs_end: int
        OPTIONAL - The end of the observation in gps time
    ondisk: boolean
        Whether to use files that are on-disk for beginning and end times. Default=False

    Returns:
    --------
    enter_files: float
        A float between 0 and 1 that describes the normalised time that the pulsar enters the beam
    exit_files: float
         a float between 0 and 1 that describes the normalised time that the pulsar exits the beam
    """
    if not metadata or not full_meta:
        metadata, full_meta = get_common_obs_metadata(obsid, return_all=True)

    #Find the beginning and end of obs
    obs_beg, obs_end = files_beg, files_end = obs_max_min(obsid)
    obs_dur = obs_end-obs_beg + 1

    #Logic loop:
    if ondisk==True:
        #find the beginning and end time of the observation FILES you have on disk
        files_beg, files_end = find_combined_beg_end(obsid)
        files_duration = files_end - files_beg + 1
    else:
        #uses manually input beginning and end times to find beam coverage
        files_beg = beg
        files_end = end
        files_duration = files_end - files_beg + 1

    #find the enter and exit times of pulsar normalized with the observing time
    names_ra_dec = grab_source_alog(pulsar_list=[pulsar], query=query)
    beam_source_data, _ = find_sources_in_obs([obsid], names_ra_dec, min_power=min_z_power, metadata_list=[[metadata, full_meta]])
    if beam_source_data[obsid]:
        enter_obs_norm = beam_source_data[obsid][0][1]
        exit_obs_norm = beam_source_data[obsid][0][2]
    else:
        logger.warn("{} not in beam".format(pulsar))
        return None, None

    #times the source enters and exits beam
    time_enter = obs_beg + obs_dur*enter_obs_norm -1
    time_exit = obs_beg + obs_dur*exit_obs_norm -1

    #normalised time the source enters/exits the beam in the files
    enter_files = (time_enter-files_beg)/files_duration
    exit_files = (time_exit-files_beg)/files_duration

    if enter_files<0.:
        enter_files=0.
    if exit_files>1.:
        exit_files=1.
    if enter_files>1. or exit_files<0.:
        logger.debug("source {0} is not in the beam for the files on disk".format(pulsar))
        enter_files = None
        exit_files = None

    return enter_files, exit_files


def ATNF_spectral_data_plot(pulsar_list):
    """
    Given a list of pulsars, plots the available spectral energy distribution for each one using data from ATNF

    Parameters:
    -----------
    pulsar_list: list
        A list of the J names of pulsars to plot
    """
    K=None
    covar_mat=None
    a_err=None
    for pulsar in pulsar_list:
        nu_atnf, S_atnf, S_atnf_e, a, a_err = flux_from_atnf(pulsar)
        if not nu_atnf:
            logger.warn("{}: No data available on ATNF database. Cannot plot.".format(pulsar))
            continue
        if not a or a_err:
            a, a_err, K, covar_mat = find_spind(pulsar, nu_atnf, S_atnf, S_atnf_e)

        plot_flux_estimation(pulsar, nu_atnf, S_atnf, S_atnf_e, a,\
                        K=K, covar_mat=covar_mat, a_err=a_err)


def least_squares_fit_plaw(x_data, y_data, y_err):
    """
    Used primarily by est_pulsar_flux() to to attain a power law function. Intended for use with pulsar flux densities

    Parameters:
    -----------
    x_data: list
        A list of frequencies in Hz
    y_data: list
        A list of fluxes in Jy correspodning to the input x_data frequencies
    y_err: list
        A list containing the corresponding error values for y_data in Hz

    Returns:
    --------
    a: float
        The fit spectral index
    c: float
        The fit y-intercept
    covar_mat: np.matrix
        The covariance matrix of the fit. Contains the information required to for uncertainty calculations

    """
    #convert everything to numpy arrays
    x_data = np.array(x_data)
    y_data = np.array(y_data)
    y_err = np.array(y_err)

    #Set up matrices in log space
    Y = np.log(y_data)
    X = np.vstack((np.ones(len(x_data)), np.log(x_data))).T

    #Set up errors. We will use the avg. length of the errors in logspace
    Y_err = 0.5 * np.log((y_data + y_err)/(y_data - y_err))

    #Convert the errors to weights
    W = np.diag(1/Y_err**2)

    #for reference: https://en.wikipedia.org/wiki/Weighted_least_squares
    # B =(X'*W*X)' * X'WY
    XW     = np.matmul(X.T, W)
    XWX    = np.matmul(XW, X)
    XWY    = np.matmul(XW, Y)
    XWXinv = np.linalg.pinv(XWX)
    b      = np.matmul(XWXinv, XWY)

    c = b[0]
    a = b[1]
    #convert c from log tospace to linear space
    K = np.exp(c)

    #The covariance matrix
    covar_mat = XWXinv

    logger.debug("a: {0}".format(a))
    logger.debug("K: {0}".format(K))
    logger.debug("Covariance Matrix: {0}".format(covar_mat))

    return a, K, covar_mat


def flux_from_plaw(freq, K, a, covar_mat):
    """
    Calculates the flux and error from a power law fit by extrapolating to the desired frequency.
    The power law is of the form S = c * nu**a

    Parameters:
    -----------
    freq: float
        The frequency for which we want to calculate a flux for (nu)
    K: float
        The value of K from the power law function
    a: float
        The value of a (spectral index) from the power law function
    covar_matrix: numpy matrix
        The covariance matrix from our power law fit. The main diagonal elements are sigma_c^2, sigma_a^2 respectively

    Returns:
    --------
    flux: float
        The calculated flux
    flux_err: float
        The uncertainty of the calculated flux
    """

    def plaw_func(nu, K, a):
        #Power law function
        return K*nu**a

    flux = plaw_func(freq, K, a)

    #Calculate the error. For reference: https://en.wikipedia.org/wiki/Propagation_of_uncertainty
    log_freq = np.log(freq)

    dc2 = covar_mat[0, 0]
    da2 = covar_mat[1, 1]
    dac = covar_mat[0, 1]

    flux_err_log = np.sqrt(dc2 + da2*log_freq**2 + 2*dac*log_freq)

    #convert the error to linear space. We will use the 'average' logspace error
    z = np.exp(2*flux_err_log)
    flux_err = flux*(z-1)/(z+1)

    return flux, flux_err


def flux_from_spind(nu_1, nu_2, s_2, s_2_err, a, a_err):
    """
    Calculates a flux value based on the spectral index, a using the formula:
    S_1 = nu_1^a * nu_2 ^-a * S_2
    Unvcertainty in nu is negligable

    Parameters:
    -----------
    nu_1: float
        The frequency at the desired flux estimation (Hz)
    nu_2: float
        The frequency at the known flux value (Hz)
    s_2: float
        The known flux value (Jy)
    S_2_err: float
        The uncertainty in the known flux value (Jy)
    a: float
        The spectral index
    a_err: float
        The uncertainty in the spectral index

    Returns:
    --------
    flux_est: float
        The estimated flux (Jy)
    flux_est_err: float
        The uncertainty in the estimated flux - calculated using the variance formula (Jy)
    """

    flux_est = nu_1**a * nu_2**(-a) * s_2
    #variance formula error est
    s_2_var = nu_1**a * nu_2**(-a)
    s_2_var = s_2_var**2 * s_2_err**2
    a_var = s_2 * nu_1**a * nu_2**(-a) * (np.log(nu_1)-np.log(nu_2))
    a_var = a_var**2 * a_err**2
    flux_est_err = np.sqrt(s_2_var + a_var)

    return flux_est, flux_est_err


def flux_from_atnf(pulsar, query=None):
    """
    Queries the ATNF database for flux and spectral index info on a particular pulsar at all frequencies

    Parameters:
    -----------
    pulsar: string
        The J name of the pulsar
    query: object
        OPTIONAL - The return from psrqpy.QueryATNF for this pulsar. Default: None

    Returns:
    --------
    freq_all: list
        All frequencies in Hz with flux values on ATNF
    flux_all: list
        The flux values corresponding to the freq_all list in Jy
    flux_err_all: list
        The uncertainty in the flux_all values
    spind: float
        The spectral index from ATNF, will be None if not available
    spind_err: float
        The ucnertainty in spind from ATNF, will be None if not available
    """
    if query is None:
        query = psrqpy.QueryATNF(psrs=[pulsar], loadfromdb=data_load.ATNF_LOC).pandas

    flux_queries = ["S40", "S50", "S60", "S80", "S100", "S150", "S200",\
                    "S300", "S400", "S600", "S700", "S800", "S900",\
                    "S1400", "S1600", "S2000", "S3000", "S4000", "S5000",\
                    "S6000", "S8000"]
    freq_all=[]
    flux_all=[]
    flux_err_all=[]
    #Get all available data from dataframe and check for missing values
    for flux_query in flux_queries:
        flux = query[flux_query][0]
        if not np.isnan(flux):
            #sometimes error values don't exist, causing a key error in pandas
            try:
                flux_err = query[flux_query+"_ERR"][0]
                if flux_err == 0.0:
                    logger.debug("{0} flux error for query: {1}, is zero. Assuming 20% uncertainty"\
                            .format(pulsar, flux_query))
                    flux_err = flux*0.2
            except KeyError:
                logger.debug("{0} flux error value {1}, not available. assuming 20% uncertainty"\
                            .format(pulsar, flux_query))
                flux_err = flux*0.2

            if np.isnan(flux_err):
                logger.debug("{0} flux error value for {1} not available. assuming 20% uncertainty"\
                            .format(pulsar, flux_query))
                flux_err = flux*0.2

            freq_all.append(int(flux_query.split()[0][1:])*1e6) #convert to Hz
            flux_all.append(flux*1e-3) #convert to Jy
            flux_err_all.append(flux_err*1e-3) #convert to Jy

    #Also get spectral index if it exists
    spind = query["SPINDX"][0]
    spind_err = query["SPINDX_ERR"][0]

    return freq_all, flux_all, flux_err_all, spind, spind_err


def find_spind(pulsar, freq_all, flux_all, flux_err_all):
    """
    Tries to attain a spectral index from input data.
    If this fails, uses an index of -1.4 with uncertainty of 1 as per Bates 2013.

    Parameters:
    -----------
    pulsar: string
        The J name of the pulsar
    freq_all: list
        The frequencies in Hz
    flux_all:list
        The fluxes in Jy
    flux_err_all: list
        The uncertainty in the fluxes

    Returns:
    --------
    spind: float
        The spectral index
    spind_err: float
        The uncertainty in the spectral index. Will be None if power law fit was attained
    K: float
        The K value of the power law fit. Will be None if not attained
    covar_mat: numpy.matrix
        The covariance matrix of the power law fit. Will be None if not attained
    """
    spind = None
    K = None
    covar_mat = None
    spind_err = None
    #Attempt to estimate flux
    if len(flux_all) > 1:
        logger.debug("{0} calculating power law".format(pulsar))
        for i, _ in enumerate(flux_all):
            flux_all[i] = flux_all[i]
            flux_err_all[i] = flux_err_all[i]

        #Find params from least squares fit
        spind, K, covar_mat = least_squares_fit_plaw(freq_all, flux_all, flux_err_all)
        logger.debug("{0} derived spectral index: {1} +/- {2}".format(pulsar, spind, np.sqrt(abs(covar_mat.item(3)))))

    #Do something different if there is only one flux value in archive
    elif len(flux_all) == 1:
        logger.debug("{} Only a single flux value available".format(pulsar))
        if spind and not spind_err:
            logger.debug("{} spectral index error not available. Assuming 20% error".format(pulsar))
            spind_err = spind*0.2
        if not spind:
            logger.debug("{} insufficient data to estimate spectral index. Using alpha=-1.4 +/- 1.0 as per Bates2013".format(pulsar))
            spind = -1.4
            spind_err = 1.

    elif len(flux_all) < 1:
        logger.debug("{} no flux values. Cannot estimate flux. Will return Nones".format(pulsar))

    return spind, spind_err, K, covar_mat


def est_pulsar_flux(pulsar, obsid, plot_flux=False, metadata=None, query=None):
    """
    Estimates a pulsar's flux from archival data by assuming a power law relation between flux and frequency

    Parameters:
    -----------
    pulsar: string
        The puslar's name. e.g. 'J2241-5236'
    obsid: int
        The observation ID
    plot_flux: boolean
        OPTIONAL - Whether or not to make a plot of the flux estimation. Default = False
    metadata: list
        OPTIONAL - The metadata call for this obsid
    query: object
        OPTIONAL - The return from psrqpy.QueryATNF for this pulsar

    Returns:
    -------
    flux: float
        The estimated flux in Jy
    flux_err: float
        The estimated flux's uncertainty in Jy
    """

    if metadata is None:
        logger.debug("obtaining mean freq from obs metadata")
        metadata = get_common_obs_metadata(obsid)
    f_mean = metadata[5]*1e6

    freq_all, flux_all, flux_err_all, spind, spind_err = flux_from_atnf(pulsar, query=query)

    logger.debug("Freqs: {0}".format(freq_all))
    logger.debug("Fluxes: {0}".format(flux_all))
    logger.debug("Flux Errors: {0}".format(flux_err_all))
    logger.debug("{0} there are {1} flux values available on the ATNF database"\
                .format(pulsar, len(flux_all)))

    if not spind or spind_err:
        spind, spind_err, K, covar_mat = find_spind(pulsar, freq_all, flux_all, flux_err_all)

    if K and covar_mat is not None and spind:
        flux_est, flux_est_err = flux_from_plaw(f_mean, K, spind, covar_mat)
    elif spind and spind_err:
        flux_est, flux_est_err = flux_from_spind(f_mean, freq_all[0], flux_all[0], flux_err_all[0],\
                                                spind, spind_err)
    else:
        logger.debug("{} no flux values on archive. Cannot estimate flux. Will return Nones".format(pulsar))
        return None, None

    if plot_flux==True:
            plot_flux_estimation(pulsar, freq_all, flux_all, flux_err_all, spind,
                                my_nu=f_mean, my_S=flux_est, my_S_e=flux_est_err, obsid=obsid,
                                a_err=spind_err,  K=K, covar_mat=covar_mat)

    return flux_est, flux_est_err


def find_pulsar_w50(pulsar, query=None):
    """
    Attempts to find a pulsar's W50 from the ATNF catalogue. If unavailable, will estimate

    Paramters:
    ---------
    pulsar: string
        The J-name of the pulsar. e.g. 'J2241-5236'
    query: object
        OPTIONAL - The query for 'pulsar' returned from psrqpy.QueryATNF
    Returns:
    --------
    w50: float
        The pulsar's W50
    w50_err: float
        The W50's uncertainty
    """
    if query is None:
        #returns W_50 and error for a pulsar from the ATNF archive IN SECONDS
        logger.debug("Accessing ATNF database")
        query = psrqpy.QueryATNF(psrs=[pulsar], loadfromdb=data_load.ATNF_LOC).pandas

    W_50 = query["W50"][0]
    W_50_err = query["W50_ERR"][0]
    if np.isnan(W_50):
        W_50=np.nan
        W_50_err=np.nan
    else:
        #convert to seconds
        W_50 = W_50/1000.

    if np.isnan(W_50_err) and not np.isnan(W_50):
        logger.debug("{} W_50 error not on archive. returning standard 5% error".format(pulsar))
        W_50_err = W_50*0.05
    else:
        #convert to seconds
        W_50_err = W_50_err/1000.

    if np.isnan(W_50):
        logger.debug("{} applying estimated W_50. Uncertainty will be inflated".format(pulsar))
        #Rankin1993 - W = x*P^0.5 where x=4.8+/-0.5 degrees of rotation at 1GHz
        #We will nflate this error due to differing frequencies and pulsar behaviour. W_50_err=1. degrees
        coeff = 4.8
        coeff_err = 2.
        period = float(query["P0"][0])

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


def find_times(obsid, pulsar, beg, end, metadata=None, full_meta=None, min_z_power=0.3, query=None):
    """
    Find the total integration time of a pulsar in the primary beam of an obsid

    Parameters:
    ----------
    obsid: int
        The observation ID
    pulsar: string
        The J name of the pulsar
    beg: int
        The beginning of the processed observing time
    end: int
        The end of the processed observing time

    Returns:
    -------
    enter_time: int
        The time when the pulsar enters the beam in gps
    exit_time: int
        The time when the pulsar exits the beam in gps
    t_int: int
        The total time that the pulsar is in the beam in seconds
    """
    #type assurances
    obsid = int(obsid)

    if not metadata or not full_meta:
        metadata, full_meta = get_common_obs_metadata(obsid, return_all=True)
    obs_beg, obs_end = obs_max_min(obsid)
    enter_norm, exit_norm = pulsar_beam_coverage(obsid, pulsar, beg, end, metadata=metadata, full_meta=full_meta, min_z_power=min_z_power, query=query)
    if beg is not None and end is not None:
        dur = end-beg
    else: #use entire obs duration
        beg = obs_beg
        end = obs_end
        dur = end - beg + 1
    if enter_norm is not None and exit_norm is not None:
        enter_time = beg + enter_norm * dur
        exit_time = beg + exit_norm * dur
        t_int = dur*(exit_norm-enter_norm)
    else:
        logger.debug("Integration time calculation failed")
        enter_time = 0
        exit_time = 0
        t_int=0

    return enter_time, exit_time, t_int


def find_t_sys_gain(pulsar, obsid, beg, end, p_ra=None, p_dec=None, enter=None, t_int=None,\
                    obs_metadata=None, full_meta=None, query=None, min_z_power=0.3, trcvr=data_load.TRCVR_FILE):

    """
    Finds the system temperature and gain for an observation.
    A function snippet originally written by Nick Swainston - adapted for general VCS use.

    Parameters:
    -----------
    pulsar: str
        the J name of the pulsar. e.g. J2241-5236
    obsid: int
        The observation ID. e.g. 1226406800
    beg: int
        The beginning of the observing time
    end: float
        The end of the observing time
    enter: float
        OPTIONAL - The fractional time that the pulsar enters the beam wrt beg and end
    t_int: float
        OPTINOAL - The frctional time that the pulsar is in the beam wrt beg and end
    p_ra: str
        OPTIONAL - the target's right ascension
    p_dec: str
        OPTIONAL - the target's declination
    obs_metadata: list
        OPTIONAL - the array generated from get_common_obs_metadata(obsid)
    query: object
        OPTIONAL - The return of the psrqpy function for this pulsar
    trcvr: str
        The location of the MWA receiver temp csv file. Default = <vcstools_data_dir>MWA_Trcvr_tile_56.csv

    Returns:
    --------
    t_sys: float
        The system temperature
    t_sys_err: float
        The system temperature's uncertainty
    gain: float
        The system gain
    gain_err: float
        The gain's uncertainty
    """

    #get ra and dec if not supplied
    if p_ra is None or p_dec is None and query is None:
        logger.debug("Obtaining pulsar RA and Dec from ATNF")
        query = psrqpy.QueryATNF(psrs=[pulsar], loadfromdb=data_load.ATNF_LOC).pandas
        p_ra = query["RAJ"][0]
        p_dec = query["DECJ"][0]
    elif p_ra is None and p_dec is None and query is not None:
        p_ra = query["RAJ"][0]
        p_dec= query["DECJ"][0]

    #get metadata if not supplied
    if not obs_metadata or not full_meta:
        logger.debug("Obtaining obs metadata")
        obs_metadata, full_meta = get_common_obs_metadata(obsid, return_all=True)

    obsid, obs_ra, obs_dec, _, delays, centrefreq, channels = obs_metadata

    #get enter time
    if not enter or not t_int:
        logger.debug("Calculating beginning time for pulsar coverage")
        enter, _, t_int = find_times(obsid, pulsar, beg, end, metadata=obs_metadata, full_meta=full_meta, min_z_power=min_z_power, query=query)

    start_time =  enter-int(obsid)

    #Get important info
    trec_table = Table.read(trcvr,format="csv")
    ntiles = 128 #TODO actually we excluded some tiles during beamforming, so we'll need to account for that here

    beam_power = get_beam_power_over_time([obsid, obs_ra, obs_dec, t_int, delays,\
                                           centrefreq, channels],\
                                           np.array([[pulsar, p_ra, p_dec]]),\
                                           dt=100, start_time=start_time)
    beam_power = np.mean(beam_power)

    # Usa a primary beam function to convolve the sky temperature with the primary beam
    # prints suppressed
    sys.stdout = open(os.devnull, 'w')
    _, _, Tsky_XX, _, _, _, Tsky_YY, _ = pbtant.make_primarybeammap(int(obsid), delays, centrefreq*1e6, 'analytic', plottype='None')
    sys.stdout = sys.__stdout__

    #TODO can be inaccurate for coherent but is too difficult to simulate
    t_sky = (Tsky_XX + Tsky_YY) / 2.
    # Get T_sys by adding Trec and Tsky (other temperatures are assumed to be negligible
    t_sys_table = t_sky + get_Trec(trec_table, centrefreq)
    t_sys = np.mean(t_sys_table)
    t_sys_err = t_sys*0.02 #TODO: figure out what t_sys error is

    logger.debug("pul_ra: {} pul_dec: {}".format(p_ra, p_dec))
    _, _, zas = mwa_alt_az_za(obsid, ra=p_ra, dec=p_dec)
    theta = np.radians(zas)
    gain = from_power_to_gain(beam_power, centrefreq*1e6, ntiles, coh=True)
    logger.debug("beam_power: {} theta: {} pi: {}".format(beam_power, theta, np.pi))
    gain_err = gain * ((1. - beam_power)*0.12 + 2.*(theta/(0.5*np.pi))**2. + 0.1)

    return t_sys, t_sys_err, gain, gain_err


def est_pulsar_sn(pulsar, obsid, beg, end,
                 p_ra=None, p_dec=None, obs_metadata=None, full_meta=None, plot_flux=False,\
                 query=None, min_z_power=0.3, trcvr=data_load.TRCVR_FILE):

    """
    Estimates the signal to noise ratio for a pulsar in a given observation using the radiometer equation
    S/N = (s_mean * gain * sqrt(n_p * t_int * df * (period - W_50)/W_50)) / t_sys
    Note that W_50 should be W_equiv but we can't figure that out so we're estimating

    Parameters:
    ----------
    pulsar: string
        Name of the pulsar e.g. J2241-5236
    obsid: int
        Observation ID e.g. 1226406800
    beg: int
        beginning of the observing time
    end: int
        end of the observing time
    p_ra: str
        OPTIONAL - the target's right ascension
    p_dec: str
        OPTIONAL - the target's declination
    obs_metadata: list
        OPTIONAL - the array generated from get_common_obs_metadata(obsid)
    plot_flux: boolean
        OPTIONAL - whether or not to produce a plot of the flux estimation. Default = False

    Returns:
    --------
    sn: float
        The expected signal to noise ratio for the given inputs
    sn_err: float
        The uncertainty in the signal to noise ratio
    """
    #We will attain uncertainties for s_mean, gain, t_sys and W_50.
    # other uncertainties are considered negligible
    logger.info("hello")
    if query is None:
        query = psrqpy.QueryATNF(psrs=pulsar, loadfromdb=data_load.ATNF_LOC).pandas

    if p_ra is None or p_dec is None:
        #Get some basic pulsar and obs info info
        p_ra = query["RAJ"][0]
        p_dec = query["DECJ"][0]

    #get metadata if not supplied
    if not obs_metadata or not full_meta:
        logger.debug("Obtaining obs metadata")
        obs_metadata, full_meta = get_common_obs_metadata(obsid, return_all=True)

    n_p = 2 #constant
    df = 30.72e6 #(24*1.28e6)

    #estimate flux
    s_mean, s_mean_err = est_pulsar_flux(pulsar, obsid, plot_flux=plot_flux,\
                         metadata=obs_metadata, query=query)
    #fluxes may be Nones. If so, return None
    if s_mean is None and s_mean_err is None:
        return None, None, None, None

    #find integration time
    enter, _, t_int = find_times(obsid, pulsar, beg=beg, end=end, metadata=obs_metadata, full_meta=full_meta, min_z_power=min_z_power, query=query)
    if t_int<=0.:
        logger.warning("{} not in beam for obs files or specificed beginning and end times"\
                    .format(pulsar))
        return 0., 0., 0., 0.

    #find system temp and gain
    t_sys, t_sys_err, gain, gain_err = find_t_sys_gain(pulsar, obsid, enter=enter, t_int=t_int,
                                beg=beg, end=end, p_ra=p_ra, p_dec=p_dec, query=query,\
                                obs_metadata=obs_metadata, full_meta=full_meta, trcvr=trcvr, min_z_power=min_z_power)

    #Find W_50
    W_50, W_50_err = find_pulsar_w50(pulsar, query=query)

    #calculate SN
    period = float(query["P0"][0])
    SN = ((s_mean * gain)/t_sys) * np.sqrt(n_p * t_int * df * (period - W_50)/W_50)

    #Calculate SN uncertainty using variance formula. Assuming error from period, df and t_int is zero
    dc_expr = np.sqrt((period-W_50)/W_50)
    var_s_mean = (gain * np.sqrt(n_p * t_int * df)) * dc_expr / t_sys
    var_gain = s_mean * np.sqrt(n_p * t_int * df) * dc_expr / t_sys
    var_W_50 = s_mean * gain * np.sqrt(n_p * t_int * df)/t_sys * (period/(-2.*W_50**2.)) * dc_expr**-1
    var_t_sys = -s_mean * gain * np.sqrt(n_p * t_int * df) * dc_expr / t_sys**2.
    var_s_mean      = var_s_mean**2.    * s_mean_err**2.
    var_gain        = var_gain**2.      * gain_err**2.
    var_W_50        = var_W_50**2.      * W_50_err**2.
    var_t_sys       = var_t_sys**2.     * t_sys_err**2.

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


def multi_psr_snfe(pulsar_list, obsid, beg, end,
                    obs_metadata=None, full_meta=None, plot_flux=False,\
                    query=None, min_z_power=0.3, trcvr=data_load.TRCVR_FILE):

    logger.info("""This script may use estimations where data is missing.
    For full verbosity, use the DEBUG logger (ie. -L DEBUG)""")

    if obs_metadata is None or full_meta is None:
        logger.debug("Obtaining obs metadata")
        obs_metadata, full_meta = get_common_obs_metadata(obsid, return_all=True)

    obs_beg, obs_end = obs_max_min(obsid)
    if beg is None:
        beg = obs_beg
    if end is None:
        end = obs_end

    mega_query = psrqpy.QueryATNF(psrs=pulsar_list, loadfromdb=data_load.ATNF_LOC).pandas
    sn_dict = {}
    for i, pulsar in enumerate(progress_bar(mega_query["PSRJ"], "Calculating pulsar SN: ")):
        psr_query = {}
        for key in mega_query.keys():
            psr_query[key] = [mega_query[key][i]]

        sn, sn_e, s, s_e = est_pulsar_sn(pulsar, obsid,\
                                 beg=beg, end=end, obs_metadata=obs_metadata, full_meta=full_meta, plot_flux=plot_flux,\
                                 query=psr_query, min_z_power=min_z_power, trcvr=trcvr)

        sn_dict[pulsar]=[sn, sn_e, s, s_e]

    return sn_dict