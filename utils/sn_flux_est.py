#! /usr/bin/env python3

#vcstools and mwa_search
import file_maxmin
import find_pulsar_in_obs as fpio
import mwa_metadb_utils
from mwa_pb import primarybeammap_tant as pbtant
import submit_to_database
import process_vcs

#Astropy
from astropy.table import Table
from astropy.time import Time

#scipy
from scipy.optimize import curve_fit

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

logger = logging.getLogger(__name__)

def pulsar_beam_coverage(obsid, pulsar, beg=None, end=None):
    #returns the beginning and end time as a fraction that a pulsar is in the primary beam for the obsid files
    #beg and end should only be supplied if the files are not present on the system

    #find the enter and exit times of pulsar normalized with the observing time
    names_ra_dec = fpio.grab_source_alog(pulsar_list=[pulsar])
    beam_source_data, _ = fpio.find_sources_in_obs([obsid], names_ra_dec)

    enter_obs_norm = beam_source_data[obsid][0][1]
    exit_obs_norm = beam_source_data[obsid][0][2]

    if beg is None and end is None:
        #find the beginning and end time of the observation FILES you have on disk
        files_beg, files_end = process_vcs.find_beg_end(obsid)
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
    if enter_files>1.:
        logger.warn("source {0} is not in the beam for the files on disk".format(pulsar))
    if exit_files<0.:
        logger.warn("source {0} is not in the beam for the files on the disk".format(pulsar))

    return enter_files, exit_files

def fit_plaw_psr(x_data, y_data, y_err=None, initial=-1.5):

    function=plaw_func
    #convert to logs
    log_x = []
    log_y = []
    for x in x_data:
        log_x.append(np.log(x))
    for y in y_data:
        log_y.append(np.log(y))
 
    if y_err is None:
        sol = curve_fit(function, log_x, log_y, p0=[initial])
    else:
        log_y_err = []
        for y in y_err:
            log_y_err.append(np.log(y))
        sol = curve_fit(function, log_x, log_y, sigma=log_y_err, p0=[initial])

    spind = sol[0][0]
    spind_err = sol[1][0][0]
    logger.debug("solution: {0}".format(sol))
    return spind, spind_err

def plaw_func(nu, a):
    #pass the log values of nu
    return a*nu

def est_pulsar_flux(pulsar, obsid=None, f_mean=None):

    """
    Estimates a pulsar's flux from archival data by assuming a power law relation between flux and frequency. Frist tries to attain a apectral index from the ATNF database. If this fails, try to work out a spectrla index. If this fails, uses an index of -1.4 with uncertainty of 1.

    INPPUT:
    pulsar: the puslar's name
    OPTIONAL:
    obsid: the Observation ID
    f_mean: forces a calculation of the flux using this frequency, rather than the obs central freq
        NOTE - one of either obsid or f_mean must be supplied. if both, uses f_mean
    output: 
    flux, flux_err: the expected pulsar flux and error at the central frequency of the obs in Jy
    """
    #query psrcat for spectral index
    spind_query = psrqpy.QueryATNF(params=["SPINDX"], psrs=[pulsar]).pandas
    spind = spind_query["SPINDX"][0]
    spind_err = spind_query["SPINDX_ERR"][0]

    #If spind is not available, we will have to estimate it
    if np.isnan(spind) or np.isnan(spind_err):
        logger.info("Pulsar {0} spectral index not available. Attempting to estimate...".format(pulsar))
        flux_queries = ["S40", "S50", "S60", "S80", "S100", "S150", "S200",\
                        "S300", "S400", "S600", "S700", "S800", "S900",\
                        "S1400", "S1600", "S2000", "S3000", "S4000", "S5000",\
                        "S6000", "S8000"] 
        error_queries = []
        for i in flux_queries:
                error_queries.append(i+"_ERR")

        #query psrcat for flux values
        flux_df = psrqpy.QueryATNF(params=flux_queries, psrs=[pulsar]).pandas
        freq_all=[]
        flux_all=[]
        flux_err_all=[]
        #Get all available data from dataframe
        for query in flux_queries:
            flux = flux_df[query][0]
            #sometimes error values don't exist, causing a key error in pandas
            try:
                flux_err = flux_df[query+"_ERR"][0]
            except KeyError:
                logger.warn("flux error value for {0} not available. Assuming 5% uncertainty".format(query))
                flux_err = flux*0.05

            if not np.isnan(flux) and not np.isnan(flux_err):
                freq_all.append(int(query.split()[0][1:])*1e6)
                flux_all.append(flux*1e-29)
                flux_err_all.append(flux_err*1e-29)

        logger.debug("Number of available freqs, fluxes and flux uncertainties: {0}, {1}, {2}"\
                    .format(len(freq_all), len(flux_all), len(flux_err_all)))
        logger.debug("Freqs: {0}".format(freq_all))
        logger.debug("Fluxes: {0}".format(flux_all))
        logger.debug("Flux Errors: {0}".format(flux_err_all))

        if len(flux_all) >1:
            logger.info("Fitting power law to archive data")
            spind, spind_err = fit_plaw_psr(freq_all, flux_all, flux_err_all)                    
        else:
            logger.warn("Insufficient archival data to estimate spectral index. Using alpha=-1.4 +/- 1.0 as per Bates2013")
            spind=-1.4
            spind_err=1.
        
        if f_mean is None:
            f_mean = mwa_metadb_utils.get_common_obs_metadata(obsid)[5]
        
        #flux calc. Err uses standard error propogation variance formula 
        logger.info("calculating flux using spectral index: {0} and error: {1}".format(spind, spind_err))
        flux_est = (f_mean*1e6)**spind
        flux_est_err = flux_est*spind_err 
        logger.info("Source flux estimate: {0} +/- {1} Jy".format(flux_est*10**26, flux_est_err*10**26))
        
    return flux_est, flux_est_err

def get_psr_w50(pulsar):
    
    #returns W50 and error for a pulsar from the ATNF archive IN SECONDS
    logger.debug("Accessing ATNF database")
    query = psrqpy.QueryATNF(params=["W50"], psrs=[pulsar]).pandas 
    W50 = query["W50"][0]
    W50_err = query["W50_ERR"][0]
    if np.isnan(W50):
        logger.warn("W50 is not on archive")
        W50=None
        W50_err=None
    else:
        #convert to seconds
        W50 = W50/1000.

    if np.isnan(W50_err) and not np.isnan(W50):
        logger.warn("W50 error not on archive. returning standard 2% error")
        W50_err = W50*0.02   
    else:
        #convert to seconds
        W50_err = W50_err/1000.

    if W50 is None:
        logger.warn("Applying estimated W50. Uncertainty will be inflated")
        #Rankin1993 - W = x*P^0.5 where x=4.8+/-0.5 degrees of rotation at 1GHz
        #We will nflate this error due to differing frequencies and pulsar behaviour. W50_err=1. degrees
        coeff = 4.8
        coeff_err = 2.
        period_query = prsqpy.QueryATNF(params=["P0"], psrs=[pulsar]).pandas
        period = float(preiod_query["P0"][0])
        
        #This estimation is worse for msps, add extra uncetainty if period < 50ms
        if period <0.05:
            coeff_err = 4.
    
        #calculate
        W50 = coeff*period**0.5
        W50_err = coeff_err*period**0.5 #error from variance calculation

        #W50 is now in degrees of rotation phase. Convert to seconds
        W50 = (W50/360.)*period
        W50_err = (W50_err/360.)*period

    return W50, W50_err

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
                beg, end = process_vcs.find_beg_end(obsid, base_path)
                dur = end-beg
                enter, exit = pulsar_beam_coverage(obsid, pulsar)
                if enter>=0. and exit <=1.:
                    t_int = dur*(exit-enter)
                else:
                    t_int=0.
            else:
                logger.warn("Combined folder is empty")
    
    if (beg is None or end is None) and t_int is None:
        logger.info("Files not on disk. Using duration for entire observation")
        beg, end, t_int = file_maxmin.print_minmax(obsid)
        enter, exit = pulsar_beam_coverage(obsid, pulsar, beg=beg, end=end)
        beg = beg+enter*t_int
        end = beg+exit*t_int
    
    if t_int is None:
        enter, exit = pulsar_beam_coverage(obsid, pulsar, beg=beg, end=end)
        dur = end-beg
        t_int = dur*(exit-enter)
 
    return beg, end, t_int

def find_t_sys_gain(pulsar, obsid, beg=None, end=None, t_int=None,\
                    p_ra=None, p_dec=None, obs_metadata=None,\
                    trcvr="/group/mwaops/PULSAR/MWA_Trcvr_tile_56.csv"):

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

    #get t_int, beg, end if not supplied
    if t_int is None or beg is None or end is None:
        logger.info("Calculating beginning, end and integration times for pulsar coverage")
        comp_config=config.load_config_file()
        base_path = comp_config["base_product_dir"]
        beg, end, t_int = find_times(pulsar, obsid, beg=beg, end=end, base_path=base_path)
        
    #Find 'start_time' for fpio
    obs_start, obs_end, obs_dur = file_maxmin.print_minmax(obsid)
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
    t_sys_err = t_sys*0.05 #TODO: figure out what t_sys error is


    gain = submit_to_database.from_power_to_gain(beam_power, centrefreq*1e6, ntiles, coh)
    gain_err = gain * ((1. - avg_power)*0.12 + (theta/(0.5*np.pi))*(theta/(0.5*np.pi))*2. + 0.1)
    
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
        logger.info("obtaining pulsar RA and Dec")
        name_ra_dec = fpio.get_psrcat_ra_dec([pulsar])
        ra = name_ra_dec[0][1]
        dec = name_ra_dec[0][2]
    
    #get metadata if not supplied
    if obs_metadata is None:
        logger.info("Obtaining obs metadata")
        obs_metadata = mwa_metadb_utils.get_common_obs_metadata(obsid)    
    
    n_p = 2 #constant
    df = 30.72e6 #(24*1.28e6)
    f_mean = obs_metadata[5]*1e6
    period = float(psrqpy.QueryATNF(params=["P0"], psrs=[pulsar]).pandas["P0"])

    #find integration time
    t_int, beg, end = find_times(obsid, pulsar, beg, end)
    if t_int>0.:
        s_mean = est_pulsar_flux(pulsar, obsid)
        t_sys, gain = find_t_sys_gain(pulsar, obsid)
    else:
        logger.warn("Pulsar not in beam for obs files or specificed beginning and end times")
        return 0.

    #find system temp and gain
    t_sys, t_sys_err, gain, gain_err = find_t_sys_gain(pulsar, obsid, beg=beg, end=end, t_int=t_int, obs_metadata=obs_metadata)

    #estimate flux
    s_mean, s_mean_err = est_psr_flux(pulsar, f_mean=obs_metadata[5])
    
    #Find W50
    W50, W50_err = find_psr_w50(pulsar)

    #calculate SN
    SN = (s_mean * gain * np.sqrt(n_p * t_int * df * (period - W_50)/W_50)) / t_sys
     
    #Calculate SN uncertainty using variance formula. Assuming error from period, df and t_int is zero
    var_s_mean = (gain * np.sqrt(n_p * t_int * df) * np.sqrt((period - W50)/W50))/t_sys
    var_s_mean = var_s_mean**2. * s_mean_err**2.
    var_gain = (s_mean * np.sqrt(n_p * t_int * df) * np.sqrt((period - W50)/W50))/t_sys
    var_gain = var_gain**2. * gain_err**2.
    var_W50 = (s_mean * gain * np.sqrt(n_p * t_int * df)/t_sys) * (period/(-2.*W50**2.)) * np.sqrt((period-W50)/W50)
    var_W50 = var_W50**2. * W50_err**2.
    var_t_sys = (-s_mean * gain * np.sqrt(n_p * t_int * df)) * np.sqrt((period-W50)/W50)/t_sys**2.
    var_t_sys = var_t_sys**2. * t_sys_err*2.
    SN_err = np.sqrt(var_s_mean + var_gain + var_W50 + var_t_sys)

    return SN, SN_err

if __name__ == "__main__":

    loglevels = dict(DEBUG=logging.DEBUG,\
                    INFO=logging.INFO,\
                    WARNING=logging.WARNING,\
                    ERROR=logging.ERROR)

    parser = argparse.ArgumentParser(description="""A utility file for estimating the S/N of a pulsar in an obsid""")
    parser.add_argument("-o", "--obsid", type=str, help="The Observation ID (e.g. 1221399680)")
    parser.add_argument("-p", "--pulsar", type=str, help="The pulsar's name (e.g. J2241-5236)")
    parser.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbostity level. Default: INFO")

    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.propagate = False

    SN = est_pulsar_sn(args.pulsar, args.obsid)
    logger.info("Pulsar S/N: {0}".format(SN))
