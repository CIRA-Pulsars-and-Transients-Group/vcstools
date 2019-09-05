#! /usr/bin/env python3

#vcstools and mwa_search
import check_known_pulsars
import binfinder
import file_maxmin
import find_pulsar_in_obs as fpio
import mwa_metadb_utils
from mwa_pb import primarybeammap_tant as pbtant
import submit_to_database
import process_vcs

#Astropy
from astropy.table import Table
from astropy.time import Time

#pythony stuff
import subprocess
import os
import glob
import config
import sys
import numpy

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

def est_pulsar_flux(pulsar, obsid):
    

    return flux_est

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
                beg, end = check_known_pulsars.find_beg_end(obsid, base_path)
                dur = end-beg
                enter, exit = binfinder.pulsar_beam_coverage(obsid, pulsar)
                if enter>=0. and exit <=1.:
                    t_int = dur*(exit-enter)
                else:
                    t_int=0.
            else:
                logger.warn("Combined folder is empty")
    
    if (beg is None or end is None) and t_int is None:
        logger.info("Files not on disk. Using duration for entire observation")
        #TODO use fpio instead
        beg, end, t_int = file_maxmin.print_minmax(obsid)
        enter, exit = pulsar_beam_coverage(obsid, pulsar, beg=beg, end=end)
        beg = beg+enter*t_int
        end = beg+exit*t_int
    
    if t_int is None:
        t_int = beg-end
        #TODO: add pulsar beam coverage
 
    return beg, end, t_int

def find_t_sys_gain(pulsar, obsid, beg=None, end=None, t_int=None,\
                    p_ra=None, p_dec=None, obs_metadata=None,\
                    trcvr="/group/mwaops/PULSAR/MWA_Trcvr_tile_56.csv"):

    #get ra and dec if not supplied
    if p_ra is None or p_dec is None:
        logger.info("obtaining pulsar RA and Dec")
        p_name_ra_dec = fpio.get_psrcat_ra_dec([pulsar])
        p_ra = name_ra_dec[1]
        p_dec = name_ra_dec[2]
    
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
    bandpowers = fpio.get_beam_power_over_time(obs_metadata,[[pulsar, p_ra, p_dec]],\
                                               dt=100, start_time=start_time)
    bandpowers = np.mean(bandpowers)

    #make the primary beam map function (prints suppressed)
    sys.stdout = open(os.devnull, 'w')
    beamsky_sum_XX,beam_sum_XX,Tant_XX,beam_dOMEGA_sum_XX,\
     beamsky_sum_YY,beam_sum_YY,Tant_YY,beam_dOMEGA_sum_YY =\
     pbtant.make_primarybeammap(int(obsid), delays, centrefreq*1e6, 'analytic', plottype='None')
    sys.stdout = sys.__stdout__

    #TODO can be inaccurate for coherent but is too difficult to simulate
    tant = (Tant_XX + Tant_YY) / 2.
    t_sys_table = tant + submit_to_database.get_Trec(trec_table,centrefreq)
    gain = submit_to_database.from_power_to_gain(bandpowers,centrefreq*1e6,ntiles,coh)
    t_sys = np.mean(t_sys_table)


    return t_sys, gain
    
def est_pulsar_SN(pulsar, obsid, beg=None, end=None):
    
    """
    S/N = (s_mean * gain * sqrt(n_p * t_int * df * (period - W_50)/W_50)) / T_sys
    Note that W_50 should be W_equiv but we can't figure that out so we're estimating    
    
    required:
    Pulsar: Name of the pulsar
    Obsid: Observation ID

    optional:
    beg: beginning of integration time
    end: end of integration time
        NOTE - leve beg and end blank if files are on disk
    """
    #Get some basic pulsar and obs info info
    name_ra_dec = fpio.get_psrcat_ra_dec([pulsar])
    ra = name_ra_dec[1]
    dec = name_ra_dec[2]
    n_p = 2
    t_int, beg, end = find_t_int(obsid, pulsar, beg, end, on_disk)
    if t_int>0.:
        s_mean = est_pulsar_flux(pulsar, obsid)
        t_sys, gain = find_t_sys(pulsar, obsid)
    else:
        SN=0.

    return SN

if __name__ == "__main__":

    loglevels = dict(DEBUG=logging.DEBUG,\
                    INFO=logging.INFO,\
                    WARNING=logging.WARNING,\
                    ERROR=logging.ERROR)

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.propagate = False

    parser = argparse.ArugmentParser(Description="""A utility file for estimating the S/N of a pulsar in an obsid""")
    parser.add_argument("-o", "--obsid", type=str, help="The Observation ID (e.g. 1221399680)")
    parser.add_argument("-p", "--pulsar", type=str, help="The pulsar's name (e.g. J2241-5236)")
    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbostity level. Default: INFO")
