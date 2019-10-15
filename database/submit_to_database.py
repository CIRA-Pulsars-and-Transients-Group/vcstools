#! /usr/bin/env python3

"""
Author: Nicholas Swainston
Creation Date: /05/2016

The MWA Pulsar Database was created by David Pallot and he wrote the database call functions.

This code is used to submit observations to the MWA Pulsar Database. The goal is to calculate all need values (flux density, width, scattering) for each observation and  submit them to the pulsar database without having to manually input them.
"""
__author__ = 'Nicholas Swainston'
__date__ = '2016-05-12'


import os
import argparse
import numpy as np
import subprocess
import sys
from shutil import copyfile as cp
import math
from scipy.interpolate import InterpolatedUnivariateSpline
import glob
import textwrap as _textwrap

#Astropy imports
from astropy.table import Table
from astropy.time import Time

#MWA software imports
from mwa_pb import primarybeammap_tant as pbtant
from mwa_pb import primary_beam
import mwa_pb.primarybeammap as pbl
from mwa_pulsar_client import client
from mwa_metadb_utils import get_common_obs_metadata, mwa_alt_az_za
import find_pulsar_in_obs as fpio
import sn_flux_est as snfe

import matplotlib.pyplot as plt
import logging
logger = logging.getLogger(__name__)

web_address = 'https://mwa-pawsey-volt01.pawsey.org.au'

class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _split_lines(self, text, width):
        text = _textwrap.dedent(self._whitespace_matcher.sub(' ', text).strip())
        return _textwrap.wrap(text, width)


def send_cmd(cmd): 
    output = subprocess.Popen(cmd.split(' '), stdin=subprocess.PIPE, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.STDOUT).communicate()[0].decode() 
    return output


def get_pulsar_dm_p(pulsar):
    #Gets the ra and dec from the output of PSRCAT
    cmd = 'psrcat -c dm {}'.format(pulsar)
    output = send_cmd(cmd)
    lines = output.split('\n')
    for l in lines[4:-1]:
        columns = l.split()

        if len(columns) > 1:
            dm = columns[1]
    cmd = 'psrcat -c p0 {}'.format(pulsar)
    output = send_cmd(cmd)
    lines = output.split('\n')
    for l in lines[4:-1]:
        columns = l.split()
        if len(columns) > 1:
            p = columns[1]
    return [dm, p]


def sex2deg( ra, dec):
    """
    sex2deg( ra, dec)
    ra - the right ascension in HH:MM:SS
    dec - the declination in DD:MM:SS
    
    Convert sexagesimal coordinates to degrees.
    """ 
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    c = SkyCoord( ra, dec, frame='icrs', unit=(u.hourangle,u.deg))
    
    # return RA and DEC in degrees in degrees
    return [c.ra.deg, c.dec.deg]
    
def get_from_bestprof(file_loc):
    """
    Get info from a bestprof file
    returns: [obsid, pulsar, dm, period, period_uncer, obsstart, obslength, profile, bin_num]
    """
    
    import re
    with open(file_loc,"r") as bestprof:
        lines = bestprof.readlines()
        # Find the obsid by finding a 10 digit int in the file name
        obsid = re.findall(r'(\d{10})', lines[0])[0]
        print(obsid)
        try:
            obsid = int(obsid)
        except:
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
        for p in range(len(orig_profile)):
            profile[p] = orig_profile[p] - min_prof
    return [obsid, pulsar, dm, period, period_uncer, obsstart, obslength, profile, bin_num]

def get_from_ascii(file_loc):
    with open(file_loc,"r") as f_ascii:
        lines = f_ascii.readlines()

        orig_profile = []
        for l in lines[1:]:
            orig_profile.append(float(l.split(" ")[-1]))
        bin_num = len(orig_profile)
        profile = np.zeros(bin_num)
        min_prof = min(orig_profile)
        for p in range(len(orig_profile)):
            profile[p] = orig_profile[p] - min_prof
        #maybe centre it around the pulse later        
    return [profile, bin_num]




def from_power_to_gain(powers,cfreq,n,coh=True):
    from astropy.constants import c,k_B
    from math import sqrt

    obswl = c.value/cfreq
    #for coherent
    if coh:
        coeff = obswl**2*16*n/(4*np.pi*k_B.value)
    else:
        coeff = obswl**2*16*sqrt(n)/(4*np.pi*k_B.value)
    logger.debug("Wavelength",obswl,"m")
    logger.debug("Gain coefficient:",coeff)
    SI_to_Jy = 1e-26
    return (powers*coeff)*SI_to_Jy


def get_Trec(tab,obsfreq):
    Trec = 0.0
    for r in range(len(tab)-1):
        if tab[r][0]==obsfreq:
            Trec = tab[r][1]
        elif tab[r][0] < obsfreq < tab[r+1][0]:
            Trec = ((tab[r][1] + tab[r+1][1])/2)
    if Trec == 0.0:
        logger.debug("ERROR getting Trec")
    return Trec


def sigmaClip(data, alpha=3, tol=0.1, ntrials=10):
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
            logger.info("Took {0} trials to reach tolerance".format(trial+1))
            np.seterr(**old_settings)
            return oldstd, x

        if trial + 1 == ntrials:
            logger.warning("Reached number of trials without reaching tolerance level")
            np.seterr(**old_settings)
            return oldstd, x

        oldstd = newstd


# Removed this function as this will always be calculated from the bestprof
#def enter_exit_calc(time_detection, time_obs, metadata, start=None, stop=None):


def zip_calibration_files(base_dir, cal_obsid, source_file):
    """
    Checkes that all the expected calibration files are where they should be 
    and returns the file location of the zipped file
    """
    import tarfile

    #Bandpass calibrations
    bandpass_files = glob.glob("{0}/BandpassCalibration_node0*.dat".format(base_dir))
    if len(bandpass_files) != 24:
        logger.error("Bandpass files not found. Exiting")
        quit()

    #DI Jones matricies
    DIJ_files = glob.glob("{0}/DI_JonesMatrices_node0*.dat".format(base_dir))
    if len(DIJ_files) != 24:
        logger.error("DI Jones matricies not found. Exiting")
        quit()

    #flagged txt files
    if not ( os.path.isfile("{0}/flagged_channels.txt".format(base_dir))
            and os.path.isfile("{0}/flagged_tiles.txt".format(base_dir)) ):
        logger.error("Flag files not found. Exiting")
        quit()

    #rts.in
    rts_files = glob.glob("{0}/rts_{1}*.in".format(base_dir, cal_obsid))
    if len(rts_files) == 0:
        logger.error("No rts in file. Exiting")
        quit()

    #source file
    if not os.path.isfile(source_file):
        logger.error("No source file. Exiting")
        quit()

    #zip the files
    zip_file_location = "{0}/{1}_rts_calibrator.zip".format(base_dir, cal_obsid)
    out = tarfile.open(zip_file_location, mode='w')
    for bf in bandpass_files:
        out.add(bf, arcname=bf.split("/")[-1])
    for DIJ in DIJ_files:
        out.add(DIJ, arcname=DIJ.split("/")[-1])
    for rts in rts_files:
        out.add(rts, arcname=rts.split("/")[-1])
    out.add("{0}/flagged_channels.txt".format(base_dir), arcname="flagged_channels.txt")
    out.add("{0}/flagged_tiles.txt".format(base_dir), arcname="flagged_tiles.txt")
    out.add(source_file, arcname=source_file.split("/")[-1])
    out.close()

    return zip_file_location


def flux_cal_and_submit(time_detection, time_obs, metadata, bestprof_data,
                        pul_ra, pul_dec, coh, auth,
                        trcvr="/group/mwaops/PULSAR/MWA_Trcvr_tile_56.csv"):
    """
    time_detection: the time in seconds of the dectection from the bestprof file
    time_obs: the time in seconds of the dectection from the metadata
    metadata: list from the function get_obs_metadata
    bestprof_data: list from the function get_from_bestprof
    trcvr: the file location of antena temperatures
    """

    #unpack the bestprof_data
    #[obsid, pulsar, dm, period, period_uncer, obsstart, obslength, profile, bin_num]
    obsid, pulsar, _, period, _, _, _, profile, num_bins = bestprof_data 
    period=float(period)
    num_bins=int(num_bins)

    #get r_sys and gain
    t_sys, u_t_sys, gain, u_gain = snfe.find_t_sys_gain(pulsar, obsid)
   
    #estimate S/N
    sn, u_sn, flagged_profile, w_equiv_bins, u_w_equiv_bins, w_equiv_ms, u_w_equiv_ms, scattered = snfe.analyse_pulse_prof(prof_data=profile, period=period, verbose=True)
  
    logger.debug("Profile scattered? {0}".format(scattered))
    logger.debug("S/N: {0} +/- {1}".format(sn, u_sn)) 
    if scattered == False:

  
        #final calc of the mean fluxdesnity in mJy
        S_mean = sn * t_sys / ( gain * math.sqrt(2. * float(time_detection) * bandwidth)) *\
                 math.sqrt( w_equiv_bins / (num_bins - w_equiv_bins)) * 1000.
        #constants to make uncertainty calc easier
        S_mean_cons = t_sys / ( math.sqrt(2. * float(time_detection) * bandwidth)) *\
                 math.sqrt( w_equiv_bins / (num_bins - w_equiv_bins)) * 1000. 
        u_S_mean = math.sqrt( math.pow(S_mean_cons * u_sn / gain , 2)  +\
                              math.pow(sn * S_mean_cons * u_gain / math.pow(gain,2) , 2) )  

        logger.info('Smean {0:.2f} +/- {1:.2f} mJy'.format(S_mean, u_S_mean))
    else:
        logger.info("Profile is scattered. Flux cannot be estimated")
        S_mean = None
        u_S_mean = None

    logger.debug("T_sys {0} K".format(t_sys))
    logger.debug("Gain {0} K/Jy".format(gain))

    #calc obstype
    if (maxfreq - minfreq) == 23:
        obstype = 1
    else:
        obstype = 2
        
    #calc scattering 
    scat_height = max(profile) / 2.71828
    scat_bins = 0
    for p in profile:
        if p > scat_height:
            scat_bins = scat_bins + 1
    scattering = float(scat_bins + 1) * float(period) /1000. #in s
    u_scattering = 1. * float(period) /1000. # assumes the uncertainty is one bin

    #calc sub-bands
    subbands = 1
    for b in range(len(channels)):
        if b == 0:
            continue
        if not (channels[b] - channels[b-1]) == 1:
            subbands = subbands + 1
    
    #get cal id
    if coh:
        cal_list = client.calibrator_list(web_address, auth)
        cal_already_created = False
        for c in cal_list:
            if ( c[u'observationid'] == int(args.cal_id) ) and ( c[u'caltype'] == calibrator_type ):
                cal_already_created = True
                cal_db_id = c[u'id']
        if not cal_already_created:
            cal_db_id = int(client.calibrator_create(web_address, auth,
                                                  observationid = str(args.cal_id),
                                                  caltype = calibrator_type)[u'id'])
    else:
        cal_db_id = None
    
    if S_mean is not None:
        #prevent TypeError caused by trying to format Nones given to fluxes for 
        #highly scattered pulsars
        S_mean = float("{0:.2f}".format(S_mean))
        u_S_mean = float("{0:.2f}".format(u_S_mean))
    
    try:
        client.detection_create(web_address, auth, 
                               observationid = int(obsid),
                               pulsar = str(pulsar), 
                               subband = int(subbands), 
                               coherent = coh,
                               observation_type = int(obstype),
                               calibrator = cal_db_id,
                               startcchan = int(minfreq), stopcchan = int(maxfreq), 
                               flux = S_mean,
                               flux_error = u_S_mean,
                               width = float("{0:.2f}".format(w_equiv_ms)),
                               width_error = float("{0:.2f}".format(u_w_equiv_ms)),
                               scattering = float("{0:.5f}".format(scattering)), 
                               scattering_error = float("{0:.5f}".format(u_scattering)),
                               dm = float(dm))
    except:
        logger.info("Detection already on database so updating the values")
        client.detection_update(web_address, auth, 
                               observationid = int(obsid),
                               pulsar = str(pulsar), 
                               subband = str(subbands), 
                               coherent = coh,
                               observation_type = int(obstype),
                               calibrator = cal_db_id,
                               startcchan = int(minfreq), stopcchan = int(maxfreq), 
                               flux = S_mean,
                               flux_error = u_S_mean,
                               width = float("{0:.2f}".format(w_equiv_ms)),
                               width_error = float("{0:.2f}".format(u_w_equiv_ms)),
                               scattering = float("{0:.5f}".format(scattering)), 
                               scattering_error = float("{0:.5f}".format(u_scattering)),
                               dm = float(dm))
                           
    logger.info("Observation submitted to database")
                              
    return subbands

"""
Testi set:
Run these form the vcstools directory

Scattered detection (crab):
python submit_to_database.py -o 1127939368 --cal_id 1127939368 -p J0534+2200 --bestprof tests/test_files/1127939368_J05342200.bestprof -L DEBUG 
Expected flux: 7500

Weak detection (it's a 10 min detection):
python submit_to_database.py -o 1222697776 --cal_id 1222695592 -p J0034-0721 --bestprof /group/mwaops/vcs/1222697776/pointings/00:34:08.87_-07:21:53.40/1222697776_PSR_J0034-0721.pfd.bestprof
Expected flux: 640

Medium detection:
python submit_to_database.py -o 1222697776 --cal_id 1222695592 -p J2330-2005 --bestprof /group/mwaops/xuemy/pol_census/1222697776/pfold/1222697776_PSR_J2330-2005.pfd.bestprof
Expected flux: ~180

Strong detection:
python submit_to_database.py -o 1226062160 --cal_id 1226054696 -p J2330-2005 --bestprof tests/test_files/1226062160_J2330-2005.bestprof -L DEBUG
S/N: 51.51 +/- 1.16
Flux: 156.04 +/- 34.13 mJy
"""

if __name__ == "__main__":
    # Dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING)
    
    #Argument parsing
    parser = argparse.ArgumentParser(formatter_class=LineWrapRawTextHelpFormatter, description=_textwrap.dedent("""
    This code is used to submit pulsar detections and observation calibrations to the MWA Pulsar Database. It can calculate all need values (flux density, width, scattering) for each detection and submit them to the pulsar database without having to manually input them. It can also submit diagnostic files such as pulse profiles and calibrations. This code assumes that the observation is coherent and correlated using the RTS so please use --incoh and --andre if your observation is incoherent or correlated using Andre's Tools respectively.

    A common use case is to simply upload a calibration file (this can be done before a detection is uploaded)
    > submit_to_database.py --cal_dir_to_tar <calibration directory> --srclist <srclist> -o <obs ID> -O <calibrator obs ID>
    Another is to record a pulsar detection and its paramters, such as flux density and pulse width, to the database. To do this the bestprof file is needed 
    > submit_to_database.py -b <bestprof file location>
    Diagnostic files (see Upload Options) such as PRESTO plots can be uploaded using:
    > submit_to_database.py  -o <obs ID> -O <calibrator obs ID> -p <pulsar> --ppps <PRESTO prepfold output post script file location>
    The diagnostic files can be created using the Dspsr Calculation Options (not robustly tested).
    """))
    parser.add_argument('-o','--obsid',type=str,help='The observation ID (eg. 1221399680).')
    parser.add_argument('-O', '--cal_id',type=str,help='The observation ID of the calibrator.')
    parser.add_argument('-p','--pulsar',type=str,help='The pulsar J name.')
    parser.add_argument('--incoh',action='store_true',help='Used for incoherent detections to accurately calculate gain. The default is coherent.')
    parser.add_argument('--andre',action='store_true',help="Used for calibrations done using Andre Offrina's tools. Default is RTS.")
    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO", 
                        choices=loglevels.keys(), default="INFO")
    parser.add_argument("-V", "--version", action='store_true', help="Print version and quit")
        
    calcargs = parser.add_argument_group('Detection Calculation Options', 'All the values of a pulsar detection (such as flux density, width, scattering) can be calculated by this code using either a .besprof file or a soon to be implemented DSPSR equivalent and automatically uploaded to the MWA Pulsar Database. Analysis files can still be uploaded before this step but this will leave the pulsar values as nulls.')
    calcargs.add_argument('-b','--bestprof',type=str,help='The location of the .bestprof file. Using this option will cause the code to calculate the needed parameters to be uploaded to the database (such as flux density, width and scattering). Using this option can be used instead of inputting the observation ID and pulsar name.')
    calcargs.add_argument('--ascii',type=str,help='The location of the ascii file (pulsar profile output of DSPSR). Using this option will cause the code to calculate the needed parameters to be uploaded to the database (such as flux density, width and scattering).')
    calcargs.add_argument('--start', type=int,
            help="The start time of the detection in GPS format.")
    calcargs.add_argument('--stop', type=int,
            help="The stop time of the detection in GPS format.")
    calcargs.add_argument('--trcvr',type=str, default = "/group/mwaops/PULSAR/MWA_Trcvr_tile_56.csv", help='File location of the receiver temperatures to be used. Only required if you do not want to use the default values located in %(default)s.')

    uploadargs = parser.add_argument_group('Upload Options', 'The different options for each file type that can be uploaded to the pulsar database. Will cause an error if the wrong file type is being uploaded.')
    uploadargs.add_argument('--cal_dir_to_tar',type=str,help='The calibration directory of a calibration solution that you would like to tar and upload to the database (eg. /group/mwaops/vcs/1221832280/cal/1221831856/rts). Must be used with --srclist so the correct source list is uploaded. If the calibration files are in the default positions then they will be tared and uploaded.')
    uploadargs.add_argument('--srclist',type=str,help='Used with --cal_dir to indicate the source list file location. eg /group/mwaops/vcs/1221832280/cal/1221831856/vis/srclist_pumav3_EoR0aegean_EoR1pietro+ForA_1221831856_patch1000.txt.')
    uploadargs.add_argument('-c','--calibration',type=str,help='The calibration solution file location to be uploaded to the database. Expects a single file so please zip or tar up the bandpass calibrations, the DI Jones matrices, the flagged_channels.txt file, the flagged_tiles.txt file, the rts.in file and the source file.')
    uploadargs.add_argument('-a','--archive',type=str,help="The DSPSR archive file location to be uploaded to the database. Expects a single file that is the output of DSPSR using the pulsar's ephemeris.")
    uploadargs.add_argument('--single_pulse_series',type=str,help='The single pulse series file location to be uploaded to the database. Expects a single file that is the output of DSPSR in single pulse mode (the -s option).')
    uploadargs.add_argument('--ppps',type=str,help="The Presto Prepfold PostScript file location to be uploaded to the database. Expects a single file that is the output of PRESTO's prepfold script.")
    uploadargs.add_argument('-i','--ippd',type=str,help="The Intergrated Pulse Profile (sometimes called a pulse profile) file location. Expects a single file that is the output of DSPSR's pav script.")
    uploadargs.add_argument('-w','--waterfall',type=str,help="The file location of a waterfall plot of pulse phase vs frequency. Expects a single file that is the output of DSPSR's psrplot.")
        
    dspsrargs = parser.add_argument_group('DSPSR Calculation Options', "Requires the --fits_files option. These options are all boolean flags that when used will send off DSPSR jobs to process the needed files that can be uploaded to the database. The files will be uploaded automatically when the DSPS scripts are tested more") #TODO remove when I'm confident with dspsr
    dspsrargs.add_argument('-f','--fits_files',type=str,help='The fits files location to be used in any Detection Calculation processing. Recommended to end in *.fits and surrounded by quotation marks.', default = "/group/mwaops/vcs/${obsid}/fits/*fits")
    dspsrargs.add_argument('--u_archive', action='store_true',help='Used to create an archive file using DSPSR.')
    dspsrargs.add_argument('--u_single_pulse_series', action='store_true',help='Used to create a single pulse archive using DSPSR.')
    dspsrargs.add_argument('--u_ppps', action='store_true', help="Used to create a Presto Prepfold PostScript file using PRESTO")
    dspsrargs.add_argument('--u_ippd', action='store_true', help="Used to create an Intergrated Pulse Profile (sometimes called a pulse profile) using DSPSR")
    dspsrargs.add_argument('--u_waterfall', action='store_true', help="Used to create a waterfall plot of pulse phase vs frequency using DSPSR.")
    args=parser.parse_args()


    if args.version:
        try:
            import version
            logger.info(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            logger.error("Couldn't import version.py - have you installed vcstools?")
            logger.error("ImportError: {0}".format(ie))
            sys.exit(0)

    # set up the logger for stand-alone execution
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.propagate = False

    #Checks for MWA database usernames and passwords
    if 'MWA_PULSAR_DB_USER' in os.environ and 'MWA_PULSAR_DB_PASS' in os.environ:
        auth = (os.environ['MWA_PULSAR_DB_USER'],os.environ['MWA_PULSAR_DB_PASS'])
    else:
        auth = ('mwapulsar','veovys9OUTY=')
        logging.warning("No MWA Pulsar Database username and password found so using the defaults.")
        logging.warning('Please add the following to your .bashrc: ')
        logging.warning('export MWA_PULSAR_DB_USER="<username>"')
        logging.warning('export MWA_PULSAR_DB_PASS="<password>"')
        logging.warning('replacing <username> <password> with your MWA Pulsar Database username and password.')

    #defaults for incoh and calibrator type
    if args.incoh:
        coh = False
        calibrator_type = None
    else:
        coh = True
        if not args.cal_id:
            logger.error("Please include --cal_id for coherent observations")
            quit()
        if args.andre:
            calibrator_type = 1
        else:
            calibrator_type = 2

    if args.fits_files:
        fits_files_loc = args.fits_files
    else:
        fits_files_loc = '/scratch2/mwaops/vcs/'+str(obsid)+'/fits/*.fits'

    #get info from .bestprof file
    if args.bestprof:
        bestprof_data = get_from_bestprof(args.bestprof)
        obsid, pulsar, dm, period, period_uncer, obsstart, time_detection,\
               profile, num_bins = bestprof_data
        if obsid is None and args.obsid:
            obsid = args.obsid
        elif obsid is None:
            logger.error("Please use --obsid. Exiting")
            sys.exit(1)
    elif args.ascii:
        profile, num_bins = get_from_ascii(args.ascii)
        if args.obsid:
            obsid = args.obsid
        if args.pulsar:
            pulsar = args.pulsar
        dm, period = get_pulsar_dm_p(pulsar)
        period_uncer = float(period) / 10.**10 #assuming it's miniscule
        if args.stop and args.start:
            time_detection = args.stop - args.start
        else:
            logger.error("Please use --start and --stop for ascii files. Exiting")
            sys.exit(1)
        bestprof_data = [obsid, pulsar, dm, period, period_uncer, args.start,
                         time_detection, profile, num_bins]
    elif args.obsid and args.pulsar:
        num_bins = 128 #used in dspsr calculations
        obsid = args.obsid
        pulsar = args.pulsar
    elif args.obsid and (args.calibration or args.cal_dir_to_tar) and not\
         (args.archive or args.single_pulse_series or args.ppps or args.ippd\
          or args.waterfall or args.u_archive or args.u_single_pulse_series\
          or args.u_ppps or args.u_ippd  or args.u_waterfall):
        #uploaded calibrator
        obsid = args.obsid
    else:
        logger.error("Please us either (--obsid and --pulsar) or --bestprof")
        quit()
        
    if args.pulsar or args.bestprof or args.ascii:
        #Checks to see if the pulsar is already on the database
        pul_list_dict = client.pulsar_list(web_address, auth)
        pul_list_str = ''
        for p in pul_list_dict:
            pul_list_str = pul_list_str + p[u'name']
        if pulsar in pul_list_str:
            logger.info('This pulsar is already on the database')
            
            #gets Ra and DEC from PSRCAT
            pulsar_ra_dec = fpio.get_psrcat_ra_dec(pulsar_list=[pulsar])
            pulsar_name, pul_ra, pul_dec = pulsar_ra_dec[0]
        else:
            logger.info('Congratulations you have detected ' + pulsar + ' for the first time with the MWA')
            #gets Ra and DEC from PSRCAT
            pulsar_ra_dec = fpio.get_psrcat_ra_dec(pulsar_list=[pulsar])
            pulsar_name, pul_ra, pul_dec = pulsar_ra_dec[0]

            #then adds it to the database
            client.pulsar_create(web_address, auth, name = pulsar, ra = pul_ra, dec = pul_dec)  
        


        

    #get meta data from obsid
    metadata = get_common_obs_metadata(obsid)
    obsid,ra_obs,dec_obs,time_obs,delays,centrefreq,channels = metadata
    minfreq = float(min(channels))
    maxfreq = float(max(channels))

    bandwidth = 30720000. #In Hz
    #bandwidth = 30720000. / 12. #In Hz

    #calc obstype
    if (maxfreq - minfreq) == 23:
        obstype = 1
    else:
        obstype = 2


    if args.bestprof or args.ascii:
        #Does the flux calculation and submits the results to the MWA pulsar database
        subbands = flux_cal_and_submit(time_detection, time_obs, metadata, bestprof_data,
                            pul_ra, pul_dec, coh, auth,
                            trcvr=args.trcvr)

    if args.cal_dir_to_tar:
        if not args.srclist:
            logger.error("You must use --srclist to define the srclist file location. Exiting")
            quit()
        args.calibration = zip_calibration_files(args.cal_dir_to_tar, args.cal_id, args.srclist)


    if args.pulsar and not (args.bestprof or args.ascii):  
        #calc sub-bands
        subbands = 1
        for b in range(len(channels)):
            if b == 0:
                continue
            if not (channels[b] - channels[b-1]) == 1:
                subbands = subbands + 1
                
        if coh:
            cal_list = client.calibrator_list(web_address, auth)
            cal_already_created = False
            for c in cal_list:
                if ( c[u'observationid'] == int(args.cal_id) ) and ( c[u'caltype'] == calibrator_type ):
                    cal_already_created = True
                    cal_db_id = c[u'id']
            if not cal_already_created:
                cal_db_id = client.calibrator_create(web_address, auth,
                                                      observationid = str(args.cal_id),
                                                      caltype = calibrator_type)[u'id']
        elif not args.cal_id:
            cal_db_id = None
            
        #uploads files to database if there's the no calc option
        #checks if the observation is on the database
        try:
            temp_dict = client.detection_get(web_address, auth, observationid = str(obsid))
        except:
            client.detection_create(web_address, auth, 
                                    observationid = str(obsid),
                                    pulsar = str(pulsar),
                                    calibrator = int(cal_db_id),
                                    subband = int(subbands),
                                    coherent = coh,
                                    startcchan = int(minfreq), stopcchan = int(maxfreq), 
                                    observation_type = int(obstype))  
            temp_dict = client.detection_get(web_address, auth, observationid =str(obsid))  
        
        if not temp_dict:
            #no obsid so creats a blank one and assumes the subbands are continuous
            client.detection_create(web_address, auth, 
                                    observationid = str(obsid),
                                    pulsar = str(pulsar),
                                    calibrator = int(cal_db_id),
                                    subband = int(subbands),
                                    coherent = coh,
                                    startcchan = int(minfreq), stopcchan = int(maxfreq), 
                                    observation_type = int(obstype))  
            temp_dict = client.detection_get(web_address, auth, observationid = str(obsid)) 
        
        pulsar_dict_check = False
        for t in range(len(temp_dict)):
            if pulsar == temp_dict[t][u'pulsar']:
                subbands = temp_dict[t][u'subband']
                pulsar_dict_check = True
        
        if not pulsar_dict_check:
            client.detection_create(web_address, auth, 
                                    observationid = str(obsid),
                                    pulsar = str(pulsar),
                                    calibrator = int(cal_db_id),
                                    subband = int(subbands),
                                    coherent = coh,
                                    startcchan = int(minfreq), stopcchan = int(maxfreq), 
                                    observation_type = int(obstype))  
            temp_dict = client.detection_get(web_address, auth, observationid = str(obsid))  

    #Upload analysis files to the database
    if args.bestprof:
        logger.info("Uploading bestprof file to database")
        cp(str(args.bestprof) ,str(obsid) + "_" + str(pulsar) + ".bestprof")
        d_file_loc = str(obsid) + "_" + str(pulsar) + ".bestprof"
        client.detection_file_upload(web_address, auth, 
                            observationid = str(obsid),
                            pulsar = str(pulsar), 
                            subband = int(subbands),
                            coherent = coh,
                            filetype = 5,
                            filepath = str(d_file_loc))
        os.system("rm " + d_file_loc)

    #Archive files
    if args.archive:
        logger.info("Uploading archive file to database")
        client.detection_file_upload(web_address, auth, 
                                    observationid = str(obsid),
                                    pulsar = str(pulsar), 
                                    subband = int(subbands),
                                    coherent = coh,
                                    filetype = 1,
                                    filepath = str(args.archive))

    if args.single_pulse_series:
        logger.info("Uploading single_pulse_series file to database")
        client.detection_file_upload(web_address, auth,
                                    observationid = str(obsid),
                                    pulsar = str(pulsar), 
                                    subband = int(subbands),
                                    coherent = coh,
                                    filetype = 2,
                                    filepath = str(args.single_pulse_series))

    if args.ppps:
        cp(str(args.ppps) ,str(obsid) + "_" + str(pulsar) + ".prepfold.ps")
        d_file_loc = str(obsid) + "_" + str(pulsar) + ".prepfold.ps"
        logger.info("Uploading Presto Prepfold PostScript file to database")
        client.detection_file_upload(web_address, auth, 
                            observationid = str(obsid),
                            pulsar = str(pulsar), 
                            subband = int(subbands),
                            coherent = coh,
                            filetype = 3,
                            filepath = str(d_file_loc))
        os.system("rm " + d_file_loc)
        
    if args.ippd:
        cp(str(args.ippd),str(obsid) + "_" + str(pulsar) + ".prof.ps")
        d_file_loc = str(obsid) + "_" + str(pulsar) + ".prof.ps"
        logger.info("Uploading Intergrates Pulse Profile file to database")
        client.detection_file_upload(web_address, auth, 
                            observationid = str(obsid),
                            pulsar = str(pulsar), 
                            subband = int(subbands),
                            coherent = coh,
                            filetype = 3,
                            filepath = str(d_file_loc))
        os.system("rm " + d_file_loc)
        
    if args.waterfall:
        cp(str(args.waterfall),str(obsid) + "_" + str(pulsar) + ".freq.vs.phase.ps")
        d_file_loc = str(obsid) + "_" + str(pulsar) + ".freq.vs.phase.ps"
        logger.info("Uploading waterfall file to database")
        client.detection_file_upload(web_address, auth, 
                            observationid = str(obsid),
                            pulsar = str(pulsar), 
                            subband = int(subbands),
                            coherent = coh,
                            filetype = 3,
                            filepath = str(d_file_loc))
        os.system("rm " + d_file_loc)
        
    if args.calibration:
        cal_list = client.calibrator_list(web_address, auth)
        cal_already_created = False
        for c in cal_list:
            if ( c[u'observationid'] == int(args.cal_id) ) and ( c[u'caltype'] == calibrator_type ):
                cal_already_created = True
                cal_db_id = c[u'id']
        if not cal_already_created:
            cal_db_id = client.calibrator_create(web_address, auth,
                                                  observationid = str(args.cal_id),
                                                  caltype = calibrator_type)#[u'id']
        
        if args.andre:
            cp(str(args.calibration),str(args.cal_id) + "_andre_calibrator.bin")
            cal_file_loc = str(args.cal_id) + "_andre_calibrator.bin"
        else:
            cp(str(args.calibration),str(args.cal_id) + "_rts_calibrator.tar")
            cal_file_loc = str(args.cal_id) + "_rts_calibrator.tar"
        
        logger.info("Uploading calibration solution to database")
        client.calibrator_file_upload(web_address, auth, 
                                       observationid = str(args.cal_id),
                                       caltype = calibrator_type, 
                                       filepath = str(cal_file_loc))
        os.system("rm " + cal_file_loc)
        
        #result = client.detection_find_calibrator(web_address, auth,detection_obsid = 35)
        
        
                                    
        


    if args.u_archive or args.u_single_pulse_series or args.u_ppps or args.u_ippd  or args.u_waterfall:
        if not glob.glob(fits_files_loc):
            logger.error("No fits files in given location. Use -f to provide file location.")
            quit()
        
        #runs all needed jobs to create all files
        from job_submit import submit_slurm
        dspsr_batch = "{0}_{1}".format(obsid, pulsar)
        commands = []
        n_omp_threads = 20
        commands.append("ncpus={0}".format(n_omp_threads))
        commands.append("export OMP_NUM_THREADS={0}".format(n_omp_threads))
        commands.append("psrcat -e {0} > {0}.par".format(pulsar))
        commands.append("fits=({0})".format(fits_files_loc))
        if args.archive:
            commands.append("ar_loc=" + str(args.archive)[:-3])
        elif args.single_pulse_series:
            commands.append("ar_loc=" + str(args.single_pulse_series)[:-3])
        elif args.u_ippd  or args.u_waterfall or args.u_archive:
            commands.append("srun -n 1 -c $ncpus dspsr -U 600 -E {0}.par -b {1} -A -cont -O {2}_{0} {3}".format(pulsar, num_bins, obsid, fits_files_loc))
            commands.append("ar_loc={0}_{1}".format(obsid,pulsar))
        if args.u_ippd:
            commands.append("pav -CDFTp -N1,1 -g {0}_{1}.prof.ps/cps ".format(obsid,pulsar) +\
                                "${ar_loc}.ar")
        if args.u_waterfall:
            commands.append('psrplot -pG -jCDTp -j "B {0}" -D {1}_{2}.freq.vs.phase.ps/cps '.\
                                format(num_bins,obsid,pulsar)+ "${ar_loc}.ar") 
        if args.u_ppps:
            commands.append("psrcat -e {0} > {0}.eph".format(pulsar))
            commands.append("srun -n 1 -c $ncpus prepfold -ncpus $ncpus -o {0} -topo -runavg -noclip -par {1}.eph -nsub 256 ".format(obsid, pulsar) + "${fits}")
            commands.append("rm {0}.eph".format(pulsar))
        if args.u_single_pulse_series:
            commands.append("srun -n 1 dspsr -U 600 -E {0}.par -b {1} -cont -s -K ".\
                                format(pulsar,num_bins) + "${fits}")
            commands.append('psraddstring="psradd -o {0}_{1}.ts.ar "'.format(obsid,pulsar))
            commands.append("ts=(pulse*.ar)")
            commands.append('for ((i=0;i<${#ts[@]};i++)); do psraddstring=${psraddstring}" "${ts[i]} ; done')
            commands.append("srun -n 1 -c $ncpus $psraddstring")
            commands.append("rm pulse*.ar")
        commands.append("rm {0}.par".format(pulsar))
        job_id = submit_slurm(dspsr_batch, commands,
                              batch_dir="./",
                              slurm_kwargs={"time": "6:50:00", "partition": "workq"},
                              submit=True)
    
