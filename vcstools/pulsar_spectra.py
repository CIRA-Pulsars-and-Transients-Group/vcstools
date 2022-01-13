"""Functions to download pulsar fluxes and fit their spectra
"""

import logging
import os
import numpy as np
import psrqpy

#matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#vcstools
from vcstools import data_load

logger = logging.getLogger(__name__)


def plot_flux_estimation(pulsar, nu_atnf, S_atnf, S_atnf_e, a,
                         my_nu=None, my_S=None, my_S_e=None, obsid=None,
                         K=None, covar_mat=None, a_err=None):
    """Used for plotting the estimated flux density against the known flux values.
    Can plot against either a least-sqaures fit flux or a spectral-index calculated flux.
    For the former, supply the covariance matrix and the K value.
    For the latter supply the spectral index error.

    Parameters
    ----------
    pulsar : `str`
        The name of the pulsar.
    nu_atnf : `list`
        The frequencies in which the known flux values correspond to in Hz.
    S_atnf : `list`
        The known flux values in Jy.
    S_atnf_e : `list`
        The uncertainties correspodning to the known fluxes.
    a : `float`
        The spectral index.
    my_nu : `float`, optional
        The frequency you're estimating the flux at (Hz). |br| Default: `None`
    my_S : `float`, optional
        The estimated flux (Jy). |br| Default: `None`
    my_S_e : `float`, optional
        The uncertainty in the estimated flux (Jy). |br| Default: `None`
    obsid : `int`, optional
        The MWA observation ID. |br| Default: `None`
    K : `float`, optional
        The K value of the least-squares fit. Use only when a least-sqaures fit has been done. |br| Default: `None`
    covar_mat: numpy.matrix object, optional
        The covariance matrix from the least-squares fit. Use only when a least-squares fit has been done. |br| Default: `None`
    a_err : `float`, optional
        The error in the spectral index. Use only when the flux has been estimated without least-squares. |br| Default: `None`
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



def ATNF_spectral_data_plot(pulsar_list):
    """Given a list of pulsars, plots the available spectral energy distribution for each one using data from ATNF

    Parameters
    ----------
    pulsar_list : `list`
        A list of the J names of pulsars to plot
    """
    K = None
    covar_mat = None
    a_err = None
    for pulsar in pulsar_list:
        nu_atnf, S_atnf, S_atnf_e, a, a_err = flux_from_atnf(pulsar)
        if not nu_atnf:
            logger.warn("{}: No data available on ATNF database. Cannot plot.".format(pulsar))
            continue
        if not a or a_err:
            a, a_err, K, covar_mat = find_spind(pulsar, nu_atnf, S_atnf, S_atnf_e)

        plot_flux_estimation(pulsar, nu_atnf, S_atnf, S_atnf_e, a,
                             K=K, covar_mat=covar_mat, a_err=a_err)


def least_squares_fit_plaw(x_data, y_data, y_err):
    """Used primarily by est_pulsar_flux() to to attain a power law function. Intended for use with pulsar flux densities.

    Parameters
    ----------
    x_data : `list`
        A list of frequencies in Hz.
    y_data : `list`
        A list of fluxes in Jy correspodning to the input x_data frequencies.
    y_err : `list`
        A list containing the corresponding error values for y_data in Hz.

    Returns
    -------
    a : `float`
        The fit spectral index
    c : `float`
        The fit y-intercept
    covar_mat : np.matrix
        The covariance matrix of the fit. Contains the information required to for uncertainty calculations.
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
    """Calculates the flux and error from a power law fit by extrapolating to the desired frequency. The power law is of the form :math:`S = c \\nu^a`.

    Parameters
    ----------
    freq : `float`
        The frequency for which we want to calculate a flux for (:math:`\\nu`).
    K : `float`
        The value of K from the power law function.
    a : `float`
        The value of a (spectral index) from the power law function.
    covar_matrix : numpy.matrix
        The covariance matrix from our power law fit. The main diagonal elements are :math:`\\sigma_c^2, \\sigma_a^2` respectively.

    Returns
    -------
    flux : `float`
        The calculated flux (mJy).
    flux_err : `float`
        The uncertainty of the calculated flux (mJy).
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
    """Calculates a flux value based on the spectral index, a using the formula:
    :math:`S_1 = \\nu_1^a  \\nu_2 ^{-a}  S_2`
    Uncertainty in :math:`\\nu` is negligable.

    Parameters
    ----------
    nu_1 : `float`
        The frequency at the desired flux estimation (Hz).
    nu_2 : `float`
        The frequency at the known flux value (Hz).
    s_2 : `float`
        The known flux value (Jy).
    S_2_err : `float`
        The uncertainty in the known flux value (Jy).
    a : `float`
        The spectral index.
    a_err : `float`
        The uncertainty in the spectral index.

    Returns
    -------
    flux_est : `float`
        The estimated flux (Jy).
    flux_est_err : `float`
        The uncertainty in the estimated flux - calculated using the variance formula (Jy).
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
    """Queries the ATNF database for flux and spectral index info on a particular pulsar at all frequencies

    Parameters
    ----------
    pulsar : `str`
        The Jname of the pulsar.
    query : psrqpy object, optional
        A previous psrqpy.QueryATNF query. Can be supplied to prevent performing a new query.

    Returns
    -------
    freq_all : `list`
        All frequencies in Hz with flux values on ATNF.
    flux_all : `list`
        The flux values corresponding to the freq_all list in Jy.
    flux_err_all : `list`
        The uncertainty in the flux_all values.
    spind : `float`
        The spectral index from ATNF, will be None if not available.
    spind_err : `float`
        The ucnertainty in spind from ATNF, will be None if not available.
    """
    if query is None:
        query = psrqpy.QueryATNF(psrs=[pulsar], loadfromdb=data_load.ATNF_LOC).pandas
    query_id = list(query['PSRJ']).index(pulsar)

    flux_queries = ["S40", "S50", "S60", "S80", "S100", "S150", "S200",\
                    "S300", "S400", "S600", "S700", "S800", "S900",\
                    "S1400", "S1600", "S2000", "S3000", "S4000", "S5000",\
                    "S6000", "S8000"]
    freq_all=[]
    flux_all=[]
    flux_err_all=[]
    #Get all available data from dataframe and check for missing values
    for flux_query in flux_queries:
        flux = query[flux_query][query_id]
        if not np.isnan(flux):
            #sometimes error values don't exist, causing a key error in pandas
            try:
                flux_err = query[flux_query+"_ERR"][query_id]
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
    spind = query["SPINDX"][query_id]
    spind_err = query["SPINDX_ERR"][query_id]

    return freq_all, flux_all, flux_err_all, spind, spind_err


def find_spind(pulsar, freq_all, flux_all, flux_err_all):
    """Tries to attain a spectral index from input data.
    If this fails, uses an index of -1.4 with uncertainty of 1 as per Bates 2013.

    Parameters
    ----------
    pulsar : `str`
        The Jname of the pulsar.
    freq_all : `list`
        The frequencies in Hz.
    flux_all : `list`
        The fluxes in Jy.
    flux_err_all : `list`
        The uncertainty in the fluxes.

    Returns
    -------
    spind : `float`
        The spectral index.
    spind_err : `float`
        The uncertainty in the spectral index. Will be None if power law fit was attained.
    K : `float`
        The K value of the power law fit. Will be None if not attained.
    covar_mat : numpy.matrix
        The covariance matrix of the power law fit. Will be None if not attained.
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