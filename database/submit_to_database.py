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
import glob
import textwrap as _textwrap

#MWA software imports
from vcstools import data_load

from mwa_pulsar_client import client
from mwa_metadb_utils import get_common_obs_metadata
import find_pulsar_in_obs as fpio
import sn_flux_est as snfe
import prof_utils


import logging
logger = logging.getLogger(__name__)

class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _split_lines(self, text, width):
        text = _textwrap.dedent(self._whitespace_matcher.sub(' ', text).strip())
        return _textwrap.wrap(text, width)
class NoAuthError(Exception):
    """Raise when pulsar database authentication is not found"""
    pass

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

def from_power_to_gain(powers,cfreq,n,coh=True):
    from astropy.constants import c,k_B
    from math import sqrt

    obswl = c.value/cfreq
    #for coherent
    if coh:
        coeff = obswl**2*16*n/(4*np.pi*k_B.value)
    else:
        coeff = obswl**2*16*sqrt(n)/(4*np.pi*k_B.value)
    logger.debug("Wavelength {} m".format(obswl))
    logger.debug("Gain coefficient: {}".format(coeff))
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

def get_subbands(metadata):
    """
    Figures out the number of subbands in an observation

    Parameters
    ----------
    metadata: list
        The output of mwa_metadb_utils.get_common_obs_metadata()

    Returns
    -------
    subbands: int
        The number of subbands in the observation
    """
    channels = metadata[-1]
    subbands = 1
    for b, _ in enumerate(channels):
        if b == 0:
            continue
        if not (channels[b] - channels[b-1]) == 1:
            subbands = subbands + 1
    return subbands

def get_db_auth_addr():
    """
    Checks for MWA database usernames and passwords

    Returns:
    --------
    auth: tuple
        The username and password for the pulsar databse
    web_address: string
        The web address of the pulsar database
    """
    web_address = 'https://pulsar-cat.icrar.uwa.edu.au/'
    if 'MWA_PULSAR_DB_USER' in os.environ and 'MWA_PULSAR_DB_PASS' in os.environ:
        auth = (os.environ['MWA_PULSAR_DB_USER'],os.environ['MWA_PULSAR_DB_PASS'])
    else:
        auth = None
        raise NoAuthError(  """
                            No MWA Pulsar Database username/password found
                            Please add the following to your .bashrc:
                            'export MWA_PULSAR_DB_USER="<username>"'
                            'export MWA_PULSAR_DB_PASS="<password>"'
                            'replacing <username> <password> with your MWA Pulsar Database username and password
                            """ )

    return web_address, auth

def check_db_and_create_det(pulsar):
    """
    Checks to see if a pulsar is already on the database. If not, will create a new entry

    Parameters:
    -----------
    puslar: str
        The name of the pulsar

    Returns:
    --------
    new_pulsar: boolean
        If True, this is a new pulsar detection
    """
    web_address, auth = get_db_auth_addr()
    pul_list_dict = client.pulsar_list(web_address, auth)
    pul_list_str = ''
    for p in pul_list_dict:
        pul_list_str = pul_list_str + p[u'name']
    if pulsar in pul_list_str:
        logger.info('This pulsar is already on the database')
        #gets Ra and DEC from PSRCAT
        pulsar_ra_dec = fpio.get_psrcat_ra_dec(pulsar_list=[pulsar])
        pulsar_name, pul_ra, pul_dec = pulsar_ra_dec[0]
        new_pulsar = False
    else:
        logger.info('Congratulations you have detected ' + pulsar + ' for the first time with the MWA')
        #gets Ra and DEC from PSRCAT
        pulsar_ra_dec = fpio.get_psrcat_ra_dec(pulsar_list=[pulsar])
        pulsar_name, pul_ra, pul_dec = pulsar_ra_dec[0]
        #then adds it to the database
        client.pulsar_create(web_address, auth, name = pulsar, ra = pul_ra, dec = pul_dec)
        new_pulsar = True
    return new_pulsar

def get_filetypes_from_db(obsid, pulsar, filetype):
    """
    Searches the pulsar database and returns the given obsid/pulsar/filetype files

    Parameters:
    -----------
    obsid: int
        Observation ID
    pulsar: string
        The name of the puslar
    filetype: int
        The type of file to search for. Options are:
        1: Archive, 2: Timeseries, 3: Diagnostics, 4: Calibration Solution, 5: Bestprof

    Returns:
    --------
    myfiles: list
        A list of filenames fitting the given parameters
    """
    def find_obj(search_key, arrdict, search_for): #for searching the various dictionaries in client
        for mydict in arrdict:
            if mydict[search_key] == search_for:
                return mydict

    #Check if any files of the filetype exist
    web_address, auth = get_db_auth_addr()
    obs_dets = client.detection_get(web_address, auth, observationid=obsid)
    if not obs_dets:
        logger.warn("No detections for this obs ID")
        return []

    detection       = find_obj("pulsar", obs_dets, pulsar)
    if not detection:
        logger.warn("No detections for this pulsar")
        return []

    allfiles_list   = detection["detection_files"]
    is_filetype     = bool(find_obj("filetype", allfiles_list, filetype))
    if not is_filetype:
        logger.warn("No files of the specified type available for this pulsar & obsid")
        return []

    #find the files with the filetype
    myfiles = []
    for detfile in allfiles_list:
        if detfile["filetype"] == filetype:
            myfiles.append(detfile["filename"])

    return myfiles

def upload_cal_files(obsid, cal_id, cal_dir_to_tar, srclist, caltype=2):
    """
    Uploads the calibrator solutions to the MWA pulsar database

    Parameters:
    -----------
    obsid: int
        The observation ID
    cal_id: int
        The calibrator ID
    cal_dir_to_tar: string
        The location of the calibrator RTS directory
    srclist: string
        The path to the sourcelist
    caltype: int
        The type of calibrator. Default: 2
    """
    auth, web_address = get_db_auth_addr()
    client_files_dict = client.calibration_file_by_observation_id(web_address, auth, obsid=cal_id)
    if client_files_dict:
        logger.info("This calibrator already has solutions on the database. Not uploading")
    else:
        zip_loc = zip_calibration_files(cal_dir_to_tar, cal_id, srclist)

        # Check if calibrator has to be made
        cal_list = client.calibrator_list(web_address, auth)
        if not cal_id in [c['observationid'] for c in cal_list]:
            client.calibrator_create(web_address, auth, observationid=str(cal_id))

        # Upload Calibration file
        try:
            client.calibrator_file_upload(web_address, auth, observationid=str(cal_id), filepath=str(zip_loc), caltype=2)
            logger.info("Uploaded calibrator solutions from {} to the database".format(cal_id))
        except:
            logger.warn("Failed to upload calibration files")
        os.remove(zip_loc)

def filename_prefix(obsid, pulsar, bins=None, cal=None):
    """
    Creates a filename prefix depending on inputs

    Parameters:
    -----------
    obsid: int
        The observation ID
    pulsar: str
        The name of the puslar
    bins: int
        The number of bins in the profile
    cal: int
        The calibrator ID

    Returns:
    --------
    pref: str
        The filename prefix
    """
    pref = "{0}_{1}".format(obsid, pulsar)
    if cal:
        pref += "_c{}".format(cal)
    if bins:
        pref += "_b{}".format(bins)
    return pref

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

def upload_file_to_db(obsid, pulsar, filepath, filetype, metadata=None, coh=True):
    """
    Uploads a file to the pulsar database

    Parameters:
    -----------
    obsid: str
        The observation ID
    pulsar: str
        The name of the pulsar
    filepath: str
        The file path of the file to upload
    filetype: int
        The type of file to upload: 1: Archive, 2: Timeseries, 3: Diagnostics, 4: Calibration Solution, 5: Bestprof
    metadata: list
        OPTINOAL - The output of mwa_metadb_utils.get_common_obs_metadata(obsid). If None, will generate it.
    coh: boolean
        OPTINOAL - Whether this is a coherent detection or not. Default: True
    """
    if not metadata:
        metadata = get_common_obs_metadata(obsid)
    subbands = get_subbands(metadata)
    web_address, auth = get_db_auth_addr()
    client.detection_file_upload(web_address, auth,
                                    observationid = str(obsid),
                                    pulsar = pulsar,
                                    filetype = int(filetype),
                                    coherent = coh,
                                    subband = subbands,
                                    filepath = filepath)

def multi_upload_files(obsid, pulsar, files_dict, metadata=None, coh=True):
    """
    Uploads any number of files to the same obsid/pulsar on the database

    Parameters:
    -----------
    obsid: str
        The observation ID
    pulsar: str
        The name of the pulsar
    files_dict: dict
        A dictionary with keys 1-5 representing filetpes. Each key holds a list containing the files to be uploaded.
    metadata: list
        OPTINOAL - The output of mwa_metadb_utils.get_common_obs_metadata(obsid). If None, will generate it.
    coh: boolean
        OPTINOAL - Whether this is a coherent detection or not. Default: True
    """
    for filetype in files_dict.keys():
        for filename in files_dict[filetype]:
            logger.info("Uploading file to databse: {}".format(filename))
            upload_file_to_db(obsid, pulsar, filename, int(filetype), metadata=metadata, coh=coh)

def flux_cal_and_submit(time_obs, metadata, bestprof_data,
                        pul_ra, pul_dec, coh, auth,
                        pulsar=None, trcvr=data_load.TRCVR_FILE):
    """
    time_obs: the time in seconds of the dectection from the metadata
    metadata: list from the function get_obs_metadata
    bestprof_data: list from the function get_from_bestprof
    trcvr: the file location of antena temperatures
    """
    #unpack the bestprof_data
    #[obsid, pulsar, dm, period, period_uncer, obsstart, obslength, profile, bin_num]
    obsid, prof_psr, dm, period, _, beg, t_int, profile, num_bins = bestprof_data
    if not pulsar:
        pulsar = prof_psr
    period=float(period)
    num_bins=int(num_bins)

    #get r_sys and gain
    t_sys, _, gain, u_gain = snfe.find_t_sys_gain(pulsar, obsid, obs_metadata=metadata,\
                                    beg=beg, end=(t_int + beg - 1))

    #estimate S/N
    try:
        prof_dict = prof_utils.auto_gfit(profile,\
                    period = period, plot_name="{0}_{1}_{2}_bins_gaussian_fit.png".format(obsid, pulsar, num_bins))
    except (prof_utils.ProfileLengthError, prof_utils.NoFitError) as _:
        prof_dict=None

    if not prof_dict:
        logger.info("Profile couldn't be fit. Using old style of profile analysis")
        prof_dict = prof_utils.auto_analyse_pulse_prof(profile, period)
        if prof_dict:
            sn = prof_dict["sn"]
            u_sn = prof_dict["sn_e"]
            w_equiv_bins = prof_dict["w_equiv_bins"]
            u_w_equiv_bins =  prof_dict["w_equiv_bins_e"]
            w_equiv_ms = period/num_bins * w_equiv_bins
            u_w_equiv_ms = period/num_bins * u_w_equiv_bins
            scattering = prof_dict["scattering"]*period/num_bins/1000 #convert to seconds
            u_scattering = prof_dict["scattering_e"]*period/num_bins/1000
            scattered = prof_dict["scattered"]
        else:
            logger.warn("Profile could not be analysed using any methods")
            sn = None
            u_sn = None
            w_equiv_bins = None
            u_w_equiv_bins =  None
            w_equiv_ms = None
            u_w_equiv_ms = None
            scattering = None
            u_scattering = None
            scattered = None
            S_mean = None
            u_S_mean = None
    else:
        sn = prof_dict["sn"]
        u_sn = prof_dict["sn_e"]
        w_equiv_bins = prof_dict["Weq"]
        u_w_equiv_bins =  prof_dict["Weq_e"]
        w_equiv_ms = period/num_bins * w_equiv_bins
        u_w_equiv_ms = period/num_bins * u_w_equiv_bins
        scattering = prof_dict["Wscat"]*period/num_bins/1000 #convert to seconds
        u_scattering = prof_dict["Wscat_e"]*period/num_bins/1000
        scattered = prof_dict["scattered"]

    if prof_dict:
        logger.info("Profile scattered? {0}".format(scattered))
        logger.info("S/N: {0} +/- {1}".format(sn, u_sn))
        logger.debug("Gain {0} K/Jy".format(gain))
        logger.debug("Equivalent width in bins: {0}".format(w_equiv_bins))
        logger.debug("T_sys: {0} K".format(t_sys))
        logger.debug("Bandwidth: {0}".format(bandwidth))
        logger.debug("Detection time: {0}".format(t_int))
        logger.debug("NUmber of bins: {0}".format(num_bins))

        if scattered == False:
            #final calc of the mean fluxdesnity in mJy
            S_mean = sn * t_sys / ( gain * math.sqrt(2. * float(t_int) * bandwidth)) *\
                    math.sqrt( w_equiv_bins / (num_bins - w_equiv_bins)) * 1000.
            #constants to make uncertainty calc easier
            S_mean_cons = t_sys / ( math.sqrt(2. * float(t_int) * bandwidth)) *\
                    math.sqrt( w_equiv_bins / (num_bins - w_equiv_bins)) * 1000.
            u_S_mean = math.sqrt( math.pow(S_mean_cons * u_sn / gain , 2)  +\
                                math.pow(sn * S_mean_cons * u_gain / math.pow(gain,2) , 2) )

            logger.info('Smean {0:.2f} +/- {1:.2f} mJy'.format(S_mean, u_S_mean))
        else:
            logger.info("Profile is scattered. Flux cannot be estimated")
            S_mean = None
            u_S_mean = None

    #calc obstype
    if (maxfreq - minfreq) == 23:
        obstype = 1
    else:
        obstype = 2

    subbands = get_subbands(metadata)
    web_address, auth = get_db_auth_addr()

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
        #prevent TypeError caused by trying to format Nones given to fluxes for highly scattered pulsars
        S_mean = float("{0:.2f}".format(S_mean))
        u_S_mean = float("{0:.2f}".format(u_S_mean))

    #format data for uploading
    if w_equiv_ms:
        w_equiv_ms = float("{0:.2f}".format(w_equiv_ms))
    if u_w_equiv_ms:
        u_w_equiv_ms = float("{0:.2f}".format(u_w_equiv_ms))
    if scattering:
        scattering = float("{0:.5f}".format(scattering))
    if u_scattering:
        u_scattering = float("{0:.5f}".format(u_scattering))

    det_kwargs = {}
    det_kwargs["observationid"]     = int(obsid)
    det_kwargs["pulsar"]            = str(pulsar)
    det_kwargs["subband"]           = int(subbands)
    det_kwargs["coherent"]          = coh,
    det_kwargs["observation_type"]  = int(obstype)
    det_kwargs["calibrator"]        = cal_db_id
    det_kwargs["startcchan"]        = int(minfreq)
    det_kwargs["stopcchan"]         = int(maxfreq)
    det_kwargs["flux"]              = S_mean
    det_kwargs["flux_error"]        = u_S_mean
    det_kwargs["width"]             = w_equiv_ms
    det_kwargs["width_error"]       = u_w_equiv_ms
    det_kwargs["scattering"]        = scattering
    det_kwargs["scattering_error"]  = u_scattering
    det_kwargs["dm"]                = float(dm)
    try:
        client.detection_create(web_address, auth, **det_kwargs)
    except:
        logger.info("Detection already on database so updating the values")
        client.detection_update(web_address, auth, **det_kwargs)

    logger.info("Observation submitted to database")
    return subbands

"""
Test set:
Run these from the vcstools/database directory

Scattered detection (crab):
python submit_to_database.py -o 1127939368 --cal_id 1127939368 -p J0534+2200 --bestprof tests/test_files/1127939368_J05342200.bestprof -L DEBUG
Expected flux: 7500

Weak detection (it's a 10 min detection):
python submit_to_database.py -o 1222697776 --cal_id 1222695592 -p J0034-0721 --bestprof ../tests/test_files/1222697776_PSR_J0034-0721.pfd.bestprof
Expected flux: 640

Medium detection:
python submit_to_database.py -o 1222697776 --cal_id 1222695592 -p J2330-2005 --bestprof tests/test_files/1222697776_PSR_J2330-2005.pfd.bestprof
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
    calcargs.add_argument('--trcvr',type=str, default=data_load.TRCVR_FILE, help='File location of the receiver temperatures to be used. Only required if you do not want to use the default values located in %(default)s.')

    uploadargs = parser.add_argument_group('Upload Options', 'The different options for each file type that can be uploaded to the pulsar database. Will cause an error if the wrong file type is being uploaded.')
    uploadargs.add_argument('--cal_dir_to_tar',type=str,help='The calibration directory of a calibration solution that you would like to tar and upload to the database (eg. /group/mwavcs/vcs/1221832280/cal/1221831856/rts). Must be used with --srclist so the correct source list is uploaded. If the calibration files are in the default positions then they will be tared and uploaded.')
    uploadargs.add_argument('--srclist',type=str,help='Used with --cal_dir to indicate the source list file location. eg /group/mwavcs/vcs/1221832280/cal/1221831856/vis/srclist_pumav3_EoR0aegean_EoR1pietro+ForA_1221831856_patch1000.txt.')
    uploadargs.add_argument('-c','--calibration',type=str,help='The calibration solution file location to be uploaded to the database. Expects a single file so please zip or tar up the bandpass calibrations, the DI Jones matrices, the flagged_channels.txt file, the flagged_tiles.txt file, the rts.in file and the source file.')
    uploadargs.add_argument('-a','--archive',type=str,help="The DSPSR archive file location to be uploaded to the database. Expects a single file that is the output of DSPSR using the pulsar's ephemeris.")
    uploadargs.add_argument('--single_pulse_series',type=str,help='The single pulse series file location to be uploaded to the database. Expects a single file that is the output of DSPSR in single pulse mode (the -s option).')
    uploadargs.add_argument('--ppps',type=str,help="The Presto Prepfold PostScript file location to be uploaded to the database. Expects a single file that is the output of PRESTO's prepfold script.")
    uploadargs.add_argument('-i','--ippd',type=str,help="The Intergrated Pulse Profile (sometimes called a pulse profile) file location. Expects a single file that is the output of DSPSR's pav script.")
    uploadargs.add_argument('-w','--waterfall',type=str,help="The file location of a waterfall plot of pulse phase vs frequency. Expects a single file that is the output of DSPSR's psrplot.")
    args = parser.parse_args()
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

    #defaults for incoh and calibrator type
    if args.incoh:
        coh = False
        calibrator_type = None
    else:
        coh = True
        if not args.cal_id:
            logger.error("Please include --cal_id for coherent observations")
            sys.exit(1)
        if args.andre:
            calibrator_type = 1
        else:
            calibrator_type = 2

    if args.incoh:
        fits_files_loc = '/group/mwavcs/vcs/{0}/incoh/*.fits'.format(args.obsid)
    else:
        #TODO add a pointing for this
        fits_files_loc = '/group/mwavcs/vcs/{0}/pointings/*.fits'.format(args.obsid)

    if not args.pulsar or not args.obsid:
        logger.error("Please supply both --obsid and --pulsar")

    #get db web address and authentication
    web_address, auth = get_db_auth_addr()

    #check if pulsar already on db, if not, create it
    if args.pulsar or args.bestprof or args.ascii:
        check_db_and_create_det(args.pulsar)

    #get meta data from obsid
    metadata = get_common_obs_metadata(args.obsid)
    _, ra_obs, _, time_obs, delays, centrefreq, channels = metadata
    minfreq = float(min(channels))
    maxfreq = float(max(channels))
    bandwidth = 30720000. #In Hz

    #calc obstype
    if (maxfreq - minfreq) == 23:
        obstype = 1
    else:
        obstype = 2

    if args.bestprof:
        #Does the flux calculation and submits the results to the MWA pulsar database
        bestprof_data = prof_utils.get_from_bestprof(args.bestprof)
    elif args.ascii:
        if not args.stop or not args.start:
            logger.error("Please supply both start and stop times of the detection for ascii files")
            sys.exit(1)
        profile, num_bins = prof_utils.get_from_ascii(args.ascii)
        dm, period = get_pulsar_dm_p(args.pulsar)
        time_detection = args.stop - args.start
        bestprof_data = [args.obsid, args.pulsar, dm, period, None, args.start,
                         time_detection, profile, num_bins]

    if args.bestprof or args.ascii:
        _, pul_ra, pul_dec = fpio.get_psrcat_ra_dec(pulsar_list=[args.pulsar])[0]
        subbands = flux_cal_and_submit(time_obs, metadata, bestprof_data,
                            pul_ra, pul_dec, coh, auth, pulsar=args.pulsar, trcvr=args.trcvr)

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
        det_kwargs = {}
        det_kwargs["observationid"]         = str(args.obsid)
        det_kwargs["pulsar"]                = args.pulsar
        det_kwargs["calibrator"]            = int(cal_db_id)
        det_kwargs["subband"]               = int(subbands)
        det_kwargs["coherent"]              = coh
        det_kwargs["startcchan"]            = int(minfreq)
        det_kwargs["stopcchan"]             = int(maxfreq)
        det_kwargs["observation_type"]      = int(obstype)
        try:
            temp_dict = client.detection_get(web_address, auth, observationid = str(args.obsid))
        except:
            client.detection_create(web_address, auth, **det_kwargs)
            temp_dict = client.detection_get(web_address, auth, observationid =str(args.obsid))

        if not temp_dict:
            #no obsid so creats a blank one and assumes the subbands are continuous
            client.detection_create(web_address, auth, **det_kwargs)
            temp_dict = client.detection_get(web_address, auth, observationid = str(args.obsid))

        pulsar_dict_check = False
        for t in range(len(temp_dict)):
            if args.pulsar == temp_dict[t][u'pulsar']:
                subbands = temp_dict[t][u'subband']
                pulsar_dict_check = True

        if not pulsar_dict_check:
            client.detection_create(web_address, auth, **det_kwargs)
            temp_dict = client.detection_get(web_address, auth, observationid = str(args.obsid))

    #Upload analysis files to the database
    #create filname prefix
    bins=None
    if args.bestprof:
        bins = prof_utils.get_from_bestprof(args.bestprof)[-1]
    fname_pref = filename_prefix(args.obsid, args.pulsar, bins=bins, cal=args.cal_id)
    upfiles_dict={"1":[], "2":[], "3":[], "4":[], "5":[]}
    remove_list=[]

    if args.bestprof:
        cp(str(args.bestprof) ,"{}.bestprof".format(fname_pref))
        d_file_loc = "{}.bestprof".format(fname_pref)
        upfiles_dict["5"].append(d_file_loc)
        remove_list.append(d_file_loc)

    if args.archive:
        cp(str(args.archive) ,"{}.ar".format(fname_pref))
        d_file_loc = "{}.ar".format(fname_pref)
        upfiles_dict["1"].append(d_file_loc)
        remove_list.append(d_file_loc)

    if args.single_pulse_series:
        logger.info("Uploading single_pulse_series file to database")
        upfiles_dict["2"].append(args.single_pulse_series)

    if args.ppps:
        cp(str(args.ppps) ,"{}.prepfold.ps".format(fname_pref))
        d_file_loc = "{}.prepfold.ps".format(fname_pref)
        upfiles_dict["3"].append(d_file_loc)
        remove_list.append(d_file_loc)

    if args.ippd:
        cp(str(args.ippd), "{}.prof.ps".format(fname_pref))
        d_file_loc = "{}.prof.ps".format(fname_pref)
        upfiles_dict["3"].append(d_file_loc)
        remove_list.append(d_file_loc)

    if args.waterfall:
        cp(str(args.waterfall), "{}.freq.vs.phase.ps".format(fname_pref))
        d_file_loc = "{}.freq.vs.phase.ps".format(fname_pref)
        upfiles_dict["3"].append(d_file_loc)
        remove_list.append(d_file_loc)

    #upload all files and remove copies
    multi_upload_files(str(args.obsid), args.pulsar, upfiles_dict, metadata=metadata, coh=coh)
    for filename in remove_list:
        os.remove(filename)

    #cal upload
    if args.calibration:
        logger.info("Uploading calibration solution to database")
        if args.andre:
            cp(str(args.calibration),str(args.cal_id) + "_andre_calibrator.bin")
            cal_file_loc = str(args.cal_id) + "_andre_calibrator.bin"
        else:
            cp(str(args.calibration),str(args.cal_id) + "_rts_calibrator.tar")
            cal_file_loc = str(args.cal_id) + "_rts_calibrator.tar"

        upload_cal_files(str(args.obsid), str(cal_id), cal_file_loc, args.srclist, caltype=calibrator_type)
        os.remove(cal_file_loc)