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

def get_sn_from_prof(prof_path):
    #Gets an estimate of S/N with error from a bestprof
    #Based on code oringally writted by Nick Swainston
    
    _, _, _, _, _ , _, prof_data, nbins = submit_to_database.get_from_bestprof(prof_path)
    #centre the profile around the max
    shift = -int(np.argmax(prof_data))+int(nbins)//2
    prof_data = np.roll(prof_data, shift)
    
    #find sigma and check if profile is scattered 
    sigma, flags = submit_to_database.sigmaClip(prof_data)    
    bot_prof_min = (max(prof_data) - min(prof_data)) * .1 + min(prof_data)
    if (np.nanmin(flags) > bot_prof_min) or ( not np.isnan(flags).any() ):
        logger.warn("The profile is highly scattered. S/N estimate may be incorrect")

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

    return sn, u_sn

def pulsar_beam_coverage(obsid, pulsar, beg=None, end=None):
    #returns the beginning and end time as a fraction that a pulsar is in the primary beam for the obsid files
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

def fit_plaw_psr(x_data, y_data, alpha_initial=-1.5, c_initial = 30., alpha_bound=None, c_bound=None):

    initial = [alpha_initial, c_initial]
    
    #make the bounds tuple
    if alpha_bound is None:
        alpha_bound = [-np.inf, np.inf]
    if c_bound is None:
        c_bound = [-np.inf, np.inf]
    
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


    a = sol[0][0]
    c = sol[0][1]
    logger.debug("a: {0}".format(a))
    logger.debug("c: {0}".format(c))
    logger.debug("Solution: {0}".format(sol))
    covar_matrix = np.matrix(sol[1])
    logger.debug("Covariance Matrix: {0}".format(covar_matrix))

    return a, c, covar_matrix 

def log_plaw_func(nu, a, c):
    #pass the log values of nu
    return np.exp(a * np.log(nu) + c)

def est_pulsar_flux(pulsar, obsid=None, f_mean=None):

    """
    Estimates a pulsar's flux from archival data by assuming a power law relation between flux and frequency. Frist tries to attain a apectral index from the ATNF database. If this fails, try to work out a spectrla index. If this fails, uses an index of -1.4 with uncertainty of 1.

    INPPUT:
    pulsar: the puslar's name
    OPTIONAL:
    obsid: the Observation ID
    f_mean: forces a calculation of the flux using this frequency, rather than the obs central freq
        NOTE - should be supplied in Hz. One of either obsid or f_mean must be supplied. if both, uses f_mean
    output: 
    flux, flux_err: the expected pulsar flux and error at the central frequency of the obs in Jy
    """
    if obsid is None and f_mean is None:
        logger.error("Please supply an obsid or mean frequency")
        sys.exit(1)

    if f_mean is None:
        logger.info("obtaining mean freq from obs metadata")
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

    df = psrqpy.QueryATNF(params=all_queries, psrs=[pulsar]).pandas
    freq_all=[]
    flux_all=[]
    flux_err_all=[]
    #Get all available data from dataframe and check for missing values
    for query in flux_queries:
        flux = df[query][0]
        #sometimes error values don't exist, causing a key error in pandas
        try:
            flux_err = df[query+"_ERR"][0]
            if flux_err == 0.0:
                logger.warn("Flux erorr for query: {0} is zero. Assuming 20% uncertainty".format(query))
                flux_err = flux*0.2
        except KeyError:
            logger.warn("flux error value for {0} not available. Assuming 20% uncertainty".format(query))
            flux_err = flux*0.2
        
        if not np.isnan(flux) and not np.isnan(flux_err):
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
    spind_query = psrqpy.QueryATNF(params=["SPINDX"], psrs=[pulsar]).pandas
    spind = spind_query["SPINDX"][0]
    spind_err = spind_query["SPINDX_ERR"][0]

    logger.info("There are {0} flux values available on the ATNF database for this pulsar".format(len(flux_all)))
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
                #if spind_error exists on cat, use it to create 5 sigma bounds
                spind_bounds = [spind-spind_err*5, spind+spind_err*5]
         
    
        spind, c, covar_matrix = fit_plaw_psr(freq_all, flux_all, alpha_initial=initial_spind, alpha_bound=spind_bounds)
        logger.info(np.shape(covar_matrix))
        logger.info(np.item(0))
        logger.info("Derived spectral index: {0} +/- {1}".format(spind, covar_matrix[0][0][0]))
        
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
            logger.warn("Spectral index error not available. Assuming 10% error")
            spind_err = spind*0.1
        if np.isnan(spind):
            logger.warn("Insufficient archival data to estimate spectral index. Using alpha=-1.4 +/- 1.0 as per Bates2013")
            spind = -1.4
            spind_err = 1.
       
        #formula for flux: 
        #S_1 = nu_1^a * nu_2 ^-a * S_2     
        nu_1 = f_mean
        nu_2 = freq_all[0] * 1e6 #convert to Hz from MHz
        s_2 = flux_all[0] * 1e-3 #convert to Jy from mJy
        s_2_err = flux_err_all[0] * 1e-3 
        a = spind
        a_err = spind_err        
        
        logger.info("calculating flux using spectral index: {0} and error: {1}".format(spind, spind_err))
        
        flux_est = nu_1**a * nu_2**(-a) * s_2 
        #variance formula error est
        s_2_var = nu_1**a * nu_2**(-a)
        a_var = s_2 * nu_1**a * nu_2**(-a) * (np.log(nu_1)-np.log(nu_2))  
        a_var = a_var**2 * a_err**2
        s_2_var = s_2_var**2 * s_2_err**2
        flux_est_err = np.sqrt(s_2_var + a_var)

    elif len(flux_all) < 1:
        logger.error("No flux values on archive. Cannot estimate flux. Terminating")
        sys.exit(1)

    
    logger.info("Source flux estimate at {0} MHz: {1} +/- {2} Jy".format(f_mean/1e6, flux_est, flux_est_err))
    
    return flux_est, flux_est_err

def find_pulsar_w50(pulsar):
    
    #returns W_50 and error for a pulsar from the ATNF archive IN SECONDS
    logger.debug("Accessing ATNF database")
    query = psrqpy.QueryATNF(params=["W50"], psrs=[pulsar]).pandas 
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
        logger.warn("W_50 error not on archive. returning standard 5% error")
        W_50_err = W_50*0.05   
    else:
        #convert to seconds
        W_50_err = W_50_err/1000.

    if W_50 is None:
        logger.warn("Applying estimated W_50. Uncertainty will be inflated")
        #Rankin1993 - W = x*P^0.5 where x=4.8+/-0.5 degrees of rotation at 1GHz
        #We will nflate this error due to differing frequencies and pulsar behaviour. W_50_err=1. degrees
        coeff = 4.8
        coeff_err = 2.
        period_query = prsqpy.QueryATNF(params=["P0"], psrs=[pulsar]).pandas
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

def find_times(obsid, pulsar, beg=None, end=None, base_path="/group/mwaops/vcs/"):

    #Find the total integration time of a pulsar in an obsid
    t_int=None
    if beg is None or end is None:
        #see if files are on disk
        logger.info("Checking if files for {0} are on disk".format(obsid))
        comb_path = base_path + "/{0}/combined".format(obsid)
        if os.path.exists(comb_path):
            #see if the combined folder is not empty
            if len(glob.glob("{0}/{1}*_ics.dat".format(comb_path, obsid)))>0:
                logger.info("Files found on disk")
                beg, end = process_vcs.find_combined_beg_end(obsid, base_path)
                dur = end-beg
                enter, exit = pulsar_beam_coverage(obsid, pulsar)
                if (enter>=0. and exit <=1.) and (enter>=0 and exit >=0):
                    t_int = dur*(exit-enter)
                else:
                    t_int=0.
            else:
                logger.warn("Combined folder is empty")
    
    if (beg is None or end is None) and t_int is None:
        logger.info("Files not on disk. Using duration for entire observation")
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

def find_t_sys_gain(pulsar, obsid, beg=None, p_ra=None, p_dec=None,\
                    obs_metadata=None, trcvr="/group/mwaops/PULSAR/MWA_Trcvr_tile_56.csv"):

    """
    A function snippet originally written by Nick Swainston - adapted for general VCS use.
    REQUIRED:
    pulsar: the J name of the pulsar. e.g. J2241-5236
    obsid: the observation ID. e.g. 1226406800

    OPTIONAL:
    beg, end: the beginning and end times of the observation. If unsupplied, will check disk files or use entire obs
    t_int: the integration time - the time the pulsar is in the beam
    p_ra, p_dec: the pulsar's RA and Dec. If unsupplied, will try to find it from metadata
    obs_metadata: metadata list from mwa_meatadb_utils. Will find if unsupplied
    trcvr: the location of the MWA receiver temp csv file

    OUTPUT:
    t_sys, t_sys_err: the system temperature and its uncertainty
    gain, gain_err: the gain and its uncertainty
    """

    #get ra and dec if not supplied
    if p_ra is None or p_dec is None:
        logger.info("Obtaining pulsar RA and Dec from ATNF")
        ra_dec_q = psrqpy.QueryATNF(params=["RAJ", "DECJ"], psrs=[pulsar]).pandas 
        p_ra = ra_dec_q["RAJ"]
        p_dec = ra_dec_q["DECJ"]
    
    #get metadata if not supplied
    if obs_metadata is None:
        logger.info("Obtaining obs metadata")
        obs_metadata = mwa_metadb_utils.get_common_obs_metadata(obsid) 
   
    obsid, obs_ra, obs_dec, obs_dur, delays, centrefreq, channels = obs_metadata

    #get beg if not supplied
    if beg is None:
        logger.info("Calculating beginning time for pulsar coverage")
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
    
def est_pulsar_sn(pulsar, obsid, beg=None, end=None, p_ra=None, p_dec=None, obs_metadata=None):
    
    """
    S/N = (s_mean * gain * sqrt(n_p * t_int * df * (period - W_50)/W_50)) / t_sys
    Note that W_50 should be W_equiv but we can't figure that out so we're estimating    
    
    required:
    Pulsar: Name of the pulsar e.g. J2241-5236
    Obsid: Observation ID e.g. 1226406800

    optional:
    beg, end: beginning of integration time. NOTE - leve beg and end blank if files are on disk
    p_ra, p_dec: the source's RA and Dec. Will try to find from the pulsar name if unsupplied
    """
    #We will attain uncertainties for s_mean, gain, t_sys and W_50.
    # other uncertainties are considered negligible


    #Get some basic pulsar and obs info info
    if p_ra is None or p_dec is None:
        logger.info("Obtaining pulsar RA and Dec from ATNF")
        name_ra_dec = fpio.get_psrcat_ra_dec([pulsar])
        p_ra = name_ra_dec[0][1]
        p_dec = name_ra_dec[0][2]
    
    #get metadata if not supplied
    if obs_metadata is None:
        logger.info("Obtaining obs metadata")
        obs_metadata = mwa_metadb_utils.get_common_obs_metadata(obsid)    
    
    n_p = 2 #constant
    df = 30.72e6 #(24*1.28e6)
    f_mean = obs_metadata[5]*1e6
    period = float(psrqpy.QueryATNF(params=["P0"], psrs=[pulsar]).pandas["P0"])

    #find integration time
    beg, end, t_int = find_times(obsid, pulsar, beg, end)
    if t_int<=0.:
        logger.warn("Pulsar not in beam for obs files or specificed beginning and end times")
        return 0., 0.

    #find system temp and gain
    t_sys, t_sys_err, gain, gain_err = find_t_sys_gain(pulsar, obsid,\
                                 beg=beg, p_ra=p_ra, p_dec=p_dec, obs_metadata=obs_metadata)
 
    #estimate flux
    s_mean, s_mean_err = est_pulsar_flux(pulsar, obsid, f_mean=f_mean)
    
    #Find W_50
    W_50, W_50_err = find_pulsar_w50(pulsar)

    #calculate SN
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

    logger.debug("S_mean: {0} +/- {1}".format(s_mean, s_mean_err))
    logger.debug("Gain: {0} +/- {1}".format(gain, gain_err))
    logger.debug("t_int: {0}".format(t_int))
    logger.debug("df: {0}".format(df))
    logger.debug("period: {0}".format(period))
    logger.debug("W_50: {0} +/- {1}".format(W_50, W_50_err))
    logger.debug("t_sys: {0} +/- {1}".format(t_sys, t_sys_err))

    return SN, SN_err

if __name__ == "__main__":

    loglevels = dict(DEBUG=logging.DEBUG,\
                    INFO=logging.INFO,\
                    WARNING=logging.WARNING,\
                    ERROR=logging.ERROR)

    parser = argparse.ArgumentParser(description="""A utility file for estimating the S/N of a pulsar in an obsid""")
    parser.add_argument("-o", "--obsid", type=int, help="The Observation ID (e.g. 1221399680)")
    parser.add_argument("-p", "--pulsar", type=str, help="The pulsar's name (e.g. J2241-5236)")
    parser.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbostity level. Default: INFO")
    parser.add_argument("-b", "--beg", type=int, default=None, help="The beginning time of observation. If None, will use beginning given by a metadata call")
    parser.add_argument("-e", "--end", type=int, default=None, help="The end time of observation. If None, will use the end given by a metadata call")
    parser.add_argument("--ondisk", action="store_true", help="Use this tag to calculate S/N for files on disk")
    parser.add_argument("--raj", type=str, default=None, help="The RA of the target. If None, will obtain from a call to ATNF")
    parser.add_argument("--decj", type=str, default=None, help="The Dec of the target. If None, will obtain from a call to ATNF")
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

    if args.ondisk==False:
        #Decide what to use as beg and end
        if args.beg==None or args.end==None:
            obsbeg, obsend = mwa_metadb_utils.obs_max_min(args.obsid)
            if args.beg==None:
                beg = obsbeg
            if args.end==None:
                end = obsend        
        else:
            beg = args.beg
            end = args.end

    #Decide what to use as ra and dec
    if args.raj==None or ags.decj==None:
        query = psrqpy.QueryATNF(params=["RAJ", "DECJ"], psrs=[args.pulsar]).pandas
        if args.raj==None:
            raj = query["RAJ"]
        if args.decj==None:
            decj = query["DECJ"]
    else:
        raj = args.raj
        decj = args.decj

    #metadata call to speed up script
    obs_metadata = mwa_metadb_utils.get_common_obs_metadata(args.obsid)

    SN, SN_err = est_pulsar_sn(args.pulsar, args.obsid,\
             beg=beg, end=end, p_ra=raj, p_dec=decj, obs_metadata=obs_metadata)
    logger.info("Pulsar S/N: {0} +/- {1}".format(SN, SN_err))
