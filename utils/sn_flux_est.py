#! /usr/bin/env python3

#other
import logging
import argparse
import os
import glob
import config
import sys
import numpy as np
import psrqpy
import pandas as pd

#matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#Astropy
from astropy.table import Table
from astropy.time import Time

#scipy
from scipy.optimize import curve_fit

#vcstools and mwa_search
from mwa_pb import primarybeammap_tant as pbtant
import file_maxmin
import find_pulsar_in_obs as fpio
import mwa_metadb_utils
import submit_to_database
import process_vcs

logger = logging.getLogger(__name__)

try:
    ATNF_LOC = os.environ['PSRCAT_FILE']
except:
    logger.warn("ATNF database could not be found on disk.")
    ATNF_LOC = None

#---------------------------------------------------------------
def plot_flux_estimation(freqs, fluxes, flux_errors, pulsar, obsid, alpha=None, c=None):

    plt.errorbar(freqs, fluxes, yerr=flux_errors, fmt="o", label="Data")
    plt.yscale("log")
    plt.xscale("log")
    plt.title("Flux Density Spectrum for {0}".format(pulsar))
    plt.xlabel("Frequency MHz")
    plt.ylabel("Flux Jy")

    if alpha is not None and c is not None:
        minfreq = min(freqs)
        maxfreq = max(freqs)
        x = np.linspace(minfreq/2, maxfreq*2, 200)
        y = np.exp(alpha * np.log(x) + c)
        plt.plot(x, y, "--r", label="Fitted Plaw")
        plt.legend()

    plt.savefig("flux_density_{0}_{1}.png".format(pulsar, obsid))

#---------------------------------------------------------------
def analyse_pulse_prof(prof_path=None, prof_data=None, period=None, verbose=False):
    """
    Gets an estimate of S/N with error from a pulse profile
    Based on code oringally writted by Nick Swainston

    Inputs:
    user must supply EITHER a bestprof path OR prof_data and period of your pulse profile
    the profile data should be a list of floats that makes the pulse profile
    verbose: a boolean that determines whether to return more detailed information. Discussed below

    Return:
    [sn, u_sn]: the signal to noise and its uncertainty
    verbose=True will return much more information in the form: 
    [sn, u_sn, flags, w_equiv_bins, u_w_equiv_bins w_equiv_ms, u_w_equiv_ms, scattered]
    """
    if prof_path is None and (prof_data is None or period is None):
        logger.warn("Insufficient information to attain SN estimate from profile. Returning Nones")
        return None, None

    if prof_data is None:
        _, _, _, period, _, _, _, prof_data, nbins = submit_to_database.get_from_bestprof(prof_path)
        nbins = float(nbins)
        period = float(period)
    else:
        nbins = len(prof_data)

    #centre the profile around the max
    shift = -int(np.argmax(prof_data))+int(nbins)//2
    prof_data = np.roll(prof_data, shift)
    
    #find sigma and check if profile is scattered 
    sigma, flags = submit_to_database.sigmaClip(prof_data)    
    bot_prof_min = (max(prof_data) - min(prof_data)) * .1 + min(prof_data)
    scattered=False
    if (np.nanmin(flags) > bot_prof_min) or ( not np.isnan(flags).any() ):
        logger.warn("The profile is highly scattered. S/N estimate cannot be calculated")
        if verbose == True:
            scattered=True
            #making a new profile with the only bin being the lowest point
            prof_min_i = np.argmin(prof_data)
            flags = []
            for fi in range(len(prof_data)):
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

    if verbose==False:
        return [sn, u_sn]

    elif scattered==False:
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
        
    return [sn, u_sn, flags, w_equiv_bins, u_w_equiv_bins, w_equiv_ms, u_w_equiv_ms, scattered]

#---------------------------------------------------------------
def pulsar_beam_coverage(obsid, pulsar, beg=None, end=None, ondisk=False):
    """
    Finds the normalised time that a pulsar is in the beam for a given obsidtime that a pulsar is in the primary beam for the obsid files
    #beg and end should only be supplied if the files are not present on the system
    #if pulsar is not in beam, returns -1, -1

    #find the enter and exit times of pulsar normalized with the observing time
    names_ra_dec = fpio.grab_source_alog(pulsar_list=[pulsar])
    beam_source_data, _ = fpio.find_sources_in_obs([obsid], names_ra_dec)

    enter_obs_norm = beam_source_data[obsid][0][1]
    exit_obs_norm = beam_source_data[obsid][0][2]

    #handle negative numbers
    if beg is not None and end is not None:
        if beg < 0. or end < 0.:
            logger.warn("Negative beg/end supplied")
            beg = None
            end = None

    if beg is None and end is None:
        #find the beginning and end time of the observation FILES you have on disk
        files_beg, files_end = process_vcs.find_combined_beg_end(obsid)
        if files_beg is None or files_end is None:
            #use entire obs
            logger.warn("Using entire obs duration")
            files_beg, files_end = mwa_metadb_utils.obs_max_min(obsid)
        files_duration = files_end - files_beg + 1

    else:
        #uses manually input beginning and end times to find beam coverage
        files_beg = beg
        files_end = end
        files_duration = files_end - files_beg + 1

    #find how long the total observation is (because that's what enter and exit uses)
    obs_beg, obs_end, obs_dur = file_maxmin.print_minmax(obsid)
    obs_dur = obs_end-obs_beg

    #times the source enters and exits beam
    time_enter = obs_beg + obs_dur*enter_obs_norm
    time_exit = obs_beg + obs_dur*exit_obs_norm

    #normalised time the source enters/exits the beam in the files
    enter_files = (time_enter-files_beg)/files_duration
    exit_files = (time_exit-files_beg)/files_duration

    if enter_files<0.:
        enter_files=0.
    if exit_files>1.:
        exit_files=1.
    if enter_files>1. or exit_files<0.:
        logger.warn("source {0} is not in the beam for the files on disk".format(pulsar))
        enter_files = -1
        exit_files = -1

    return enter_files, exit_files

#---------------------------------------------------------------
def fit_plaw_psr(x_data, y_data, alpha_initial=-1.5, c_initial = 30., alpha_bound=[-3., 0.],\
                c_bound=[0., 50.]):
    """
    Used primarily by est_pulsar_flux() to fit a power law function to input data. Intended for use with pulsar flux densities

    Parameters:
    -----------
    x_data: list
        A list of frequencies in Hz
    y_data: list
        A list of fluxes in Jy correspodning to the input x_data frequencies
    alpha_initial: float
        OPTIONAL - An initial estimate for the spectral index. Default = -1.5
    c_initial: float
        OPTIONAL - An intiial estimate for the value of the y-intercept. Default = 30.
    alpha_bound: list
        OPTIONAL - Boundary conditions for the spectral index. Default = [-3., 0.]
    c_bound: list
        OPTIONAL - Boundary conditions for the y-intercept. Default = [0., 50.] 

    Returns:
    --------
    a: float
        The fit spectral index
    c: float
        The fit y-intercept
    covar_matrix: np.matrix
        The covariance matrix of the fit. Contains the information required to for uncertainty calculations
    """
    def log_plaw_func(nu, a, c):
        #pass the log values of nu
        return np.exp(a * np.log(nu) + c)

    def check_fit(sol):
        #checks that the fit parameters are sensible 
        covar_matrix = np.matrix(sol[1])
        alpha = sol[0][0]
        c = sol[0][1]
        
        covar_0_0 = covar_matrix.item(0)
        covar_0_1 = covar_matrix.item(1)
        covar_1_0 = covar_matrix.item(2)
        covar_1_1 = covar_matrix.item(3)

        test = True
        #check the following conditions
        #alpha
        if alpha>10 or alpha < -10:
            test =  False
        #c
        elif c>1000 or c<-1000:
            test = False
        #covariance matrix
        elif covar_0_0==0. or covar_0_1==0. or covar_1_0==0.or covar_1_1==0.:
            test = False
        elif covar_0_0>1e6 or covar_0_0<-1e6 or covar_0_1>1e6 or covar_0_1<-1e6\
            or covar_1_0>1e6 or covar_1_0<-1e6 or covar_1_1>1e6 or covar_1_1<-1e6:
            test = False
        return test 

    initial = [alpha_initial, c_initial]
    
    #make the bounds tuple
    lower_bound = [alpha_bound[0], c_bound[0]]
    upper_bound = [alpha_bound[1], c_bound[1]]
    bounds = (lower_bound, upper_bound)

    #force the data into numpy float arrays because scipy likes it
    x_data = np.array(x_data, dtype="float64")
    y_data = np.array(y_data, dtype="float64")
    initial = np.array(initial, dtype="float64")
    logger.debug("x_data: {0}".format(x_data))
    logger.debug("y_data: {0}".format(y_data))
    
    #fit a line
    #sol = curve_fit(function, x_data, y_data, p0=initial)
    function=log_plaw_func
    sol = curve_fit(function, x_data, y_data, p0=initial, bounds=bounds)
    covar_matrix = np.matrix(sol[1])
    a = sol[0][0]
    c = sol[0][1]
    
    #Check the fit
    

    test_check = check_fit(sol)
    if test_check==False:
        logger.warn("Initial parameters could not fit a power law. Trying without bounds...")
        sol = curve_fit(function, x_data, y_data, p0=initial)

    test_check = check_fit(sol)
    if test_check==False:
        logger.warn("Trying without bounds and initial conditions...")
        sol = curve_fit(function, x_data, y_data)
        covar_matrix = np.matrix(sol[1])
        a = sol[0][0]
        c = sol[0][1]

    test_check = check_fit(sol)    
    if test_check==False:
        covar_matrix = np.matrix(sol[1])
        a = sol[0][0]
        c = sol[0][1]
        logger.warn("Bad power law fit. Results may be unphysical")
        logger.warn("a: {0}".format(a))
        logger.warn("c: {0}".format(c))
        logger.warn("Covariance Matrix: {0}".format(covar_matrix))
        return a, c, covar_matrix
    else:
        return a, c, covar_matrix


    covar_matrix = np.matrix(sol[1])
    a = sol[0][0]
    c = sol[0][1]
    logger.debug("a: {0}".format(a))
    logger.debug("c: {0}".format(c))
    logger.debug("Solution: {0}".format(sol))
    logger.debug("Covariance Matrix: {0}".format(covar_matrix))

    return a, c, covar_matrix 

#---------------------------------------------------------------
def est_pulsar_flux(pulsar, obsid, plot_flux=False, f_mean=None):
    """
    Estimates a pulsar's flux from archival data by assuming a power law relation between flux and frequency. Frist tries to attain a apectral index from the ATNF database. If this fails, try to work out a spectrla index. If this fails, uses an index of -1.4 with uncertainty of 1.

    Parameters:
    -----------
    pulsar: string
        The puslar's name. e.g. 'J2241-5236'
    obsid: int
        The observation ID
    plot_flux: boolean
        OPTIONAL - Whether or not to make a plot of the flux estimation. Default = False
    f_mean: float
        OPTIONAL - The frequency for which to estimate the flux in Hz. If not specified, will use the obsid mean frequqncy.
    
    Returns:
    ------- 
    flux: float
        The estimated flux in Jy
    flux_err: float
        The estimated flux's uncertainty in Jy
    """
    
    if f_mean==None:
        logger.debug("obtaining mean freq from obs metadata")
        f_mean = mwa_metadb_utils.get_common_obs_metadata(obsid)[5]*1e6
    
    #query psrcat for flux values
    flux_queries = ["S40", "S50", "S60", "S80", "S100", "S150", "S200",\
                    "S300", "S400", "S600", "S700", "S800", "S900",\
                    "S1400", "S1600", "S2000", "S3000", "S4000", "S5000",\
                    "S6000", "S8000"] 
    flux_error_queries = []
    for i in flux_queries:
            flux_error_queries.append(i+"_ERR")
    all_queries = flux_queries
    all_queries.append("SPINDX")

    #for some reason i have to do this
    if flux_queries[-1] == "SPINDX":
        flux_queries = flux_queries[:-1]

    df = psrqpy.QueryATNF(params=all_queries, psrs=[pulsar], loadfromdb=ATNF_LOC).pandas
    freq_all=[]
    flux_all=[]
    flux_err_all=[]
    #Get all available data from dataframe and check for missing values
    for query in flux_queries:
        flux = df[query][0]
        if not np.isnan(flux):
            #sometimes error values don't exist, causing a key error in pandas
            try:
                flux_err = df[query+"_ERR"][0]
                if flux_err == 0.0:
                    logger.warn("Flux error for query: {0}, pulsar {1}, is zero. Assuming 20% uncertainty".format(query, pulsar))
                    flux_err = flux*0.2
            except KeyError:
                logger.warn("flux error value for {0}, pulsar {1}, not available. assuming 20% uncertainty".format(query, pulsar))
                flux_err = flux*0.2
           
            if np.isnan(flux_err):
                logger.warn("flux error value for {0}, pulsar {1}, not available. assuming 20% uncertainty".format(query, pulsar))
                flux_err = flux*0.2

            freq_all.append(int(query.split()[0][1:])*1e6) #convert to Hz
            flux_all.append(flux*1e-3) #convert to Jy
            flux_err_all.append(flux_err*1e-3) #convert to Jy

    #Also get spectral index if it exists
    spind = df["SPINDX"]
    spind_err = df["SPINDX_ERR"]    

    logger.debug("Freqs: {0}".format(freq_all))
    logger.debug("Fluxes: {0}".format(flux_all))
    logger.debug("Flux Errors: {0}".format(flux_err_all))

    #query psrcat for spectral index
    spind_query = psrqpy.QueryATNF(params=["SPINDX"], psrs=[pulsar], loadfromdb=ATNF_LOC).pandas
    spind = spind_query["SPINDX"][0]
    spind_err = spind_query["SPINDX_ERR"][0]

    logger.info("There are {0} flux values available on the ATNF database for {1}".format(len(flux_all), pulsar))
    #Attempt to estimate flux
    if len(flux_all) > 1:
        logger.info("Fitting power law to archive data")
        for i in range(len(flux_all)):
            flux_all[i] = flux_all[i]
            flux_err_all[i] = flux_err_all[i]

        #apply spectral index bounds if an error already exists
        initial_spind = -1.5
        spind_bounds = None
        if not np.isnan(spind):
            initial_spind = spind
            if not np.isnan(spind_err):
                #if spind_error exists on cat, use it to create 5 sigma bounds for fitting
                spind_bounds = [spind-spind_err*5, spind+spind_err*5]
         
    
        spind, c, covar_matrix = fit_plaw_psr(freq_all, flux_all, alpha_initial=initial_spind, alpha_bound=spind_bounds)
        logger.info("Derived spectral index: {0} +/- {1}".format(spind, covar_matrix.item(0)))
        #plot if in debug mode
        if plot_flux==True:
            plot_flux_estimation(freq_all, flux_all, flux_err_all, pulsar, obsid, alpha=spind, c=c)        

        #flux calc.  
        flux_est = np.exp(spind*np.log(f_mean)+c) 
        #calculate error from covariance matrix
        a_mat = np.matrix([np.log(f_mean), 1])
        log_flux_err = np.sqrt( a_mat * covar_matrix * a_mat.T )    
        #to find the error, we take the average log error in linear space
        b = np.exp(log_flux_err)
        flux_est_err = flux_est/2. * (b - (1/b))
        flux_est_err = flux_est_err.item(0)

    #Do something different if there is only one flux value in archive
    elif len(flux_all) == 1:
        logger.warn("Only a single flux value available on the archive")

        if not np.isnan(spind) and np.isnan(spind_err):
            logger.warn("Spectral index error not available. Assuming 20% error")
            spind_err = spind*0.2
        if np.isnan(spind):
            logger.warn("Insufficient archival data to estimate spectral index. Using alpha=-1.4 +/- 1.0 as per Bates2013")
            spind = -1.4
            spind_err = 1.
       
        #formula for flux: 
        #S_1 = nu_1^a * nu_2 ^-a * S_2     
        nu_1 = f_mean #in Hz
        nu_2 = freq_all[0] #in Hz 
        s_2 = flux_all[0] #in Jy
        s_2_err = flux_err_all[0]  
        a = spind
        a_err = spind_err
        
        logger.debug("nu1 {0}".format(nu_1))        
        logger.debug("nu2 {0}".format(nu_2))        
        logger.debug("s2 {0}".format(s_2))        
        logger.info("calculating flux using spectral index: {0} and error: {1}".format(spind, spind_err))
        
        flux_est = nu_1**a * nu_2**(-a) * s_2 
        #variance formula error est
        s_2_var = nu_1**a * nu_2**(-a)
        a_var = s_2 * nu_1**a * nu_2**(-a) * (np.log(nu_1)-np.log(nu_2))  
        a_var = a_var**2 * a_err**2
        s_2_var = s_2_var**2 * s_2_err**2
        flux_est_err = np.sqrt(s_2_var + a_var)

    elif len(flux_all) < 1:
        logger.warn("No flux values on archive for {0}. Cannot estimate flux. Will return Nones".format(pulsar))
        return None, None

    
    logger.info("Source flux estimate at {0} MHz: {1} +/- {2} Jy".format(f_mean/1e6, flux_est, flux_est_err))
    
    return flux_est, flux_est_err

#---------------------------------------------------------------
def find_pulsar_w50(pulsar):
    """
    Attempts to find a pulsar's W50 from the ATNF catalogue. If unavailable, will estimate

    Paramters:
    ---------
    pulsar: string
        The J-name of the pulsar. e.g. 'J2241-5236'

    Returns:
    --------
    w50: float
        The pulsar's W50
    w50_err: float
        The W50's uncertainty
    """
    #returns W_50 and error for a pulsar from the ATNF archive IN SECONDS
    logger.debug("Accessing ATNF database")
    query = psrqpy.QueryATNF(params=["W50"], psrs=[pulsar], loadfromdb=ATNF_LOC).pandas 
    W_50 = query["W50"][0]
    W_50_err = query["W50_ERR"][0]
    if np.isnan(W_50):
        logger.warn("W_50 is not on archive")
        W_50=None
        W_50_err=None
    else:
        #convert to seconds
        W_50 = W_50/1000.

    if np.isnan(W_50_err) and not np.isnan(W_50):
        logger.warn("W_50 error not on archive for {0}. returning standard 5% error".format(pulsar))
        W_50_err = W_50*0.05   
    else:
        #convert to seconds
        W_50_err = W_50_err/1000.

    if W_50 is None:
        logger.warn("Applying estimated W_50 for {0}. Uncertainty will be inflated".format(pulsar))
        #Rankin1993 - W = x*P^0.5 where x=4.8+/-0.5 degrees of rotation at 1GHz
        #We will nflate this error due to differing frequencies and pulsar behaviour. W_50_err=1. degrees
        coeff = 4.8
        coeff_err = 2.
        period_query = prsqpy.QueryATNF(params=["P0"], psrs=[pulsar], loadfromdb=ATNF_LOC).pandas
        period = float(preiod_query["P0"][0])
        
        #This estimation is worse for msps, add extra uncetainty if period < 50ms
        if period <0.05:
            coeff_err = 4.
    
        #calculate
        W_50 = coeff*period**0.5
        W_50_err = coeff_err*period**0.5 #error from variance calculation

        #W_50 is now in degrees of rotation phase. Convert to seconds
        W_50 = (W_50/360.)*period
        W_50_err = (W_50_err/360.)*period

    return W_50, W_50_err

#---------------------------------------------------------------
def find_times(obsid, pulsar, beg=None, end=None, base_path="/group/mwaops/vcs/"):

    """
    Find the total integration time of a pulsar in the primary beam of an obsid

    Parameters:
    ----------
    obsid: int
        The observation ID
    pulsar: string
        The J name of the pulsar
    beg: int
        OPTIONAL - The beginning of the observing time
    end: int
        OPTINAL - The end of the observing time
    base_path: string
        OPTIONAL - The location of the system's observations. Default = '/group/mwaops/vcs/'

    Returns:
    -------
    beg: int
        The time when the pulsar enters the beam in gps
    end: int
        The time when the pulsar exits thebeam in gps
    t_int: int
        The total time that the pulsar is in the beam
    """
    if beg is None or end is None:
        logger.info("Using duration for entire observation")
        beg, end = mwa_metadb_utils.obs_max_min(obsid)
        t_int = end-beg
        enter, exit = pulsar_beam_coverage(obsid, pulsar, beg=beg, end=end)
        beg = beg+enter*t_int
        end = beg+exit*t_int
    
    if t_int is None:
        enter, exit = pulsar_beam_coverage(obsid, pulsar, beg=beg, end=end)
        if beg > 0. and end > 0.:
            dur = end-beg
        else: #use entire obs duration
            beg, end = mwa_metadb_utils.obs_max_min(obsid)
            dur = end - beg
        t_int = dur*(exit-enter)

    return beg, end, t_int

#---------------------------------------------------------------
def find_t_sys_gain(pulsar, obsid, beg=None, p_ra=None, p_dec=None,\
                    obs_metadata=None, trcvr="/group/mwaops/PULSAR/MWA_Trcvr_tile_56.csv"):

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
    p_ra: str
        OPTIONAL - the target's right ascension
    p_dec: str
        OPTIONAL - the target's declination
    obs_metadata: list
        OPTIONAL - the array generated from mwa_metadb_utils.get_common_obs_metadata(obsid)
    trcvr: str
        The location of the MWA receiver temp csv file. Default = '/group/mwaops/PULSAR/MWA_Trcvr_tile_56.csv'

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
    if p_ra is None or p_dec is None:
        logger.debug("Obtaining pulsar RA and Dec from ATNF")
        ra_dec_q = psrqpy.QueryATNF(params=["RAJ", "DECJ"], psrs=[pulsar], loadfromdb=ATNF_LOC).pandas 
        p_ra = ra_dec_q["RAJ"]
        p_dec = ra_dec_q["DECJ"]
    
    #get metadata if not supplied
    if obs_metadata is None:
        logger.debug("Obtaining obs metadata")
        obs_metadata = mwa_metadb_utils.get_common_obs_metadata(obsid) 
   
    obsid, obs_ra, obs_dec, obs_dur, delays, centrefreq, channels = obs_metadata

    #get beg if not supplied
    if beg is None:
        logger.debug("Calculating beginning time for pulsar coverage")
        comp_config=config.load_config_file()
        base_path = comp_config["base_product_dir"]
        beg, _, _ = find_times(obsid, pulsar, beg=beg, base_path=base_path)
        
    #Find 'start_time' for fpio
    obs_start, obs_end = mwa_metadb_utils.obs_max_min(obsid)
    start_time = beg-obs_start

    #Get important info    
    trec_table = Table.read(trcvr,format="csv")
    ntiles = 128 #TODO actually we excluded some tiles during beamforming, so we'll need to account for that here
    beam_power = fpio.get_beam_power_over_time(obs_metadata,[[pulsar, p_ra, p_dec]],\
                                               dt=100, start_time=start_time)
    beam_power = np.mean(beam_power)

    # Usa a primary beam function to convolve the sky temperature with the primary beam
    # (prints suppressed)
    sys.stdout = open(os.devnull, 'w')
    _, _, Tsky_XX, _, _, _, Tsky_YY, _ = pbtant.make_primarybeammap(int(obsid), delays, centrefreq*1e6, 'analytic', plottype='None')
    sys.stdout = sys.__stdout__


    #TODO can be inaccurate for coherent but is too difficult to simulate
    t_sky = (Tsky_XX + Tsky_YY) / 2.
    # Get T_sys by adding Trec and Tsky (other temperatures are assumed to be negligible
    t_sys_table = t_sky + submit_to_database.get_Trec(trec_table, centrefreq)
    t_sys = np.mean(t_sys_table)
    t_sys_err = t_sys*0.02 #TODO: figure out what t_sys error is

    alts, azs, zas = mwa_metadb_utils.mwa_alt_az_za(obsid, p_ra, p_dec)
    theta = np.radians(zas)
    logger.debug("Theta: {0}".format(theta))
    gain = submit_to_database.from_power_to_gain(beam_power, centrefreq*1e6, ntiles, coh=True)
    gain_err = gain * ((1. - beam_power)*0.12 + 2.*(theta/(0.5*np.pi))**2. + 0.1)

    #sometimes gain_err is a numpy array and sometimes it isnt so i have to to this...
    try:
        gain_err.shape
        gain_err = gain_err[0]
    except:
        pass   
 
    return t_sys, t_sys_err, gain, gain_err

#---------------------------------------------------------------
def est_pulsar_sn(pulsar, obsid,\
                 beg=None, end=None, p_ra=None, p_dec=None, obs_metadata=None, plot_flux=False,\
                enter=None, exit=None):
    
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
        OPTIONAL - beginning of the observing time
    end: int
        OPTIONAL - end of the observing time
    p_ra: str
        OPTIONAL - the target's right ascension
    p_dec: str
        OPTIONAL - the target's declination
    obs_metadata: list
        OPTIONAL - the array generated from mwa_metadb_utils.get_common_obs_metadata(obsid)
    plot_flux: boolean
        OPTIONAL - whether or not to produce a plot of the flux estimation. Default = False
    enter: float
        OPTIONAL - The normalised enter time of the pulsar's coverage in the beam (between 0 and 1)
    exit: float
        OPTIONAL - The normalised exit time of the pulsar's covreage in the beam (between 0 and 1)

    Returns:
    --------
    sn: float
        The expected signal to noise ratio for the given inputs
    sn_err: float
        The uncertainty in the signal to noise ratio
    """
    #We will attain uncertainties for s_mean, gain, t_sys and W_50.
    # other uncertainties are considered negligible

    #Get some basic pulsar and obs info info
    if p_ra is None or p_dec is None:
        logger.debug("Obtaining pulsar RA and Dec from ATNF")
        query = psrqpy.QueryATNF(params=["RAJ", "DECJ"], psrs=[pulsar], loadfromdb=ATNF_LOC).pandas
        p_ra = query["RAJ"] 
        p_dec = query["DECJ"]
    
    #get metadata if not supplied
    if obs_metadata is None:
        logger.debug("Obtaining obs metadata")
        obs_metadata = mwa_metadb_utils.get_common_obs_metadata(obsid)    

    n_p = 2 #constant
    df = 30.72e6 #(24*1.28e6)
        t_int = dur*(exit-enter)

    return beg, end, t_int

#---------------------------------------------------------------
def find_t_sys_gain(pulsar, obsid, beg=None, p_ra=None, p_dec=None,\
                    obs_metadata=None, trcvr="/group/mwaops/PULSAR/MWA_Trcvr_tile_56.csv"):

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
    p_ra: str
        OPTIONAL - the target's right ascension
    p_dec: str
        OPTIONAL - the target's declination
    obs_metadata: list
        OPTIONAL - the array generated from mwa_metadb_utils.get_common_obs_metadata(obsid)
    trcvr: str
        The location of the MWA receiver temp csv file. Default = '/group/mwaops/PULSAR/MWA_Trcvr_tile_56.csv'

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
    if p_ra is None or p_dec is None:
        logger.debug("Obtaining pulsar RA and Dec from ATNF")
        ra_dec_q = psrqpy.QueryATNF(params=["RAJ", "DECJ"], psrs=[pulsar], loadfromdb=ATNF_LOC).pandas 
        p_ra = ra_dec_q["RAJ"]
        p_dec = ra_dec_q["DECJ"]
    
    #get metadata if not supplied
    if obs_metadata is None:
        logger.debug("Obtaining obs metadata")
        obs_metadata = mwa_metadb_utils.get_common_obs_metadata(obsid) 
   
    obsid, obs_ra, obs_dec, obs_dur, delays, centrefreq, channels = obs_metadata

    #get beg if not supplied
    if beg is None:
        logger.debug("Calculating beginning time for pulsar coverage")
        comp_config=config.load_config_file()
        base_path = comp_config["base_product_dir"]
        beg, _, _ = find_times(obsid, pulsar, beg=beg, base_path=base_path)
        
    #Find 'start_time' for fpio
    obs_start, obs_end = mwa_metadb_utils.obs_max_min(obsid)
    start_time = beg-obs_start

    #Get important info    
    trec_table = Table.read(trcvr,format="csv")
    ntiles = 128 #TODO actually we excluded some tiles during beamforming, so we'll need to account for that here
    beam_power = fpio.get_beam_power_over_time(obs_metadata,[[pulsar, p_ra, p_dec]],\
                                               dt=100, start_time=start_time)
    beam_power = np.mean(beam_power)

    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.propagate = False

    if args.obsid==None or args.pulsar==None:
        logger.error("Obsid and Pulsar name must be supplied. Exiting...")
        sys.exit(1)

    #Decide what to use as ra and dec
    if args.pointing==None:
        query = psrqpy.QueryATNF(params=["RAJ", "DECJ"], psrs=[args.pulsar], loadfromdb=ATNF_LOC).pandas
        raj = query["RAJ"]
        decj = query["DECJ"]
    else:
        raj = args.pointing.split("_")[0]
        decj = args.pointing.split("_")[1]

    if args.loglvl=="DEBUG":
        plot=True
    else:
        plot=False
    
    SN, SN_err = est_pulsar_sn(args.pulsar, args.obsid,\
             beg=args.beg, end=args.end, p_ra=raj, p_dec=decj, plot_flux=plot)
    logger.info("Pulsar S/N: {0} +/- {1}".format(SN, SN_err))
