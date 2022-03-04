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
import subprocess
import sys
from shutil import copyfile as cp
import glob
import textwrap as _textwrap
import numpy as np
from astropy.table import Table
import csv

import logging
logger = logging.getLogger(__name__)

#MWA software imports
from mwa_pulsar_client import client

from vcstools import data_load
from vcstools.metadb_utils import get_common_obs_metadata
from vcstools.catalogue_utils import get_psrcat_ra_dec, get_psrcat_dm_period
from vcstools import prof_utils
from vcstools.general_utils import setup_logger
from vcstools.radiometer_equation import analyise_and_flux_cal



class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _split_lines(self, text, width):
        text = _textwrap.dedent(self._whitespace_matcher.sub(' ', text).strip())
        return _textwrap.wrap(text, width)

class NoAuthError(Exception):

    """Raise when pulsar database authentication is not found"""
    pass


def get_subbands(common_metadata):
    """
    Figures out the number of frequency subbands in an observation

    Parameters
    ----------
    common_metadata: list
        The output of mwa_metadb_utils.get_common_obs_metadata()

    Returns
    -------
    subbands : `int`
        The number of subbands in the observation
    """
    channels = common_metadata[-1]
    # Count how many subbands the frequency is split into
    subbands = 1
    for ci in range(1, len(channels)):
        if (channels[ci] - channels[ci]) != 1:
            # Previous channel is not contiguous
            subbands = subbands + 1
    return subbands

def get_db_auth_addr():
    """
    Checks for MWA database usernames and passwords

    Returns
    -------
    auth: tuple
        The username and password for the pulsar databse
    web_address : `str`
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

def check_db_and_create_pulsar(pulsar):
    """
    Checks to see if a pulsar is already on the database. If not, will create a new entry

    Parameters
    ----------
    puslar : `str`
        The name of the pulsar

    Returns
    -------
    new_pulsar : `boolean`
        If True, this is a new pulsar detection
    """
    web_address, auth = get_db_auth_addr()
    pul_list_dict = client.pulsar_list(web_address, auth)
    logger.debug("pul_list_dict: {}".format(pul_list_dict))
    pul_list_str = ''
    for p in pul_list_dict:
        pul_list_str = pul_list_str + p[u'name']
    logger.debug("pul_list_str: {}".format(pul_list_str))
    if pulsar in pul_list_str:
        logger.info('This pulsar is already on the database')
        #gets Ra and DEC from PSRCAT
        pulsar_ra_dec = get_psrcat_ra_dec(pulsar_list=[pulsar])
        _, pul_ra, pul_dec = pulsar_ra_dec[0]
        new_pulsar = False
    else:
        logger.info('Congratulations you have detected ' + pulsar + ' for the first time with the MWA')
        #gets Ra and DEC from PSRCAT
        pulsar_ra_dec = get_psrcat_ra_dec(pulsar_list=[pulsar])
        _, pul_ra, pul_dec = pulsar_ra_dec[0]
        #then adds it to the database
        client.pulsar_create(web_address, auth, name = pulsar, ra = pul_ra, dec = pul_dec)
        new_pulsar = True
    return new_pulsar

def get_filetypes_from_db(obsid, pulsar, filetype):
    """
    Searches the pulsar database and returns the given obsid/pulsar/filetype files

    Parameters
    ----------
    obsid : `int`
        Observation ID
    pulsar : `str`
        The name of the puslar
    filetype : `int`
        The type of file to search for. Options are:
        1: Archive, 2: Timeseries, 3: Diagnostics, 4: Calibration Solution, 5: Bestprof

    Returns
    -------
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

def upload_cal_files(obsid, calid, cal_dir_to_tar, srclist, caltype=2):
    """
    Uploads the calibrator solutions to the MWA pulsar database

    Parameters
    ----------
    obsid : `int`
        The observation ID
    calid : `int`
        The calibrator ID
    cal_dir_to_tar : `str`
        The location of the calibrator RTS directory
    srclist : `str`
        The path to the sourcelist
    caltype : `int`
        The type of calibrator. Default: 2
    """
    auth, web_address = get_db_auth_addr()
    client_files_dict = client.calibration_file_by_observation_id(web_address, auth, obsid=calid)
    if client_files_dict:
        logger.info("This calibrator already has solutions on the database. Not uploading")
    else:
        zip_loc = zip_calibration_files(cal_dir_to_tar, calid, srclist)

        # Check if calibrator has to be made
        cal_list = client.calibrator_list(web_address, auth)
        if not calid in [c['observationid'] for c in cal_list]:
            client.calibrator_create(web_address, auth, observationid=str(calid))

        # Upload Calibration file
        try:
            client.calibrator_file_upload(web_address, auth, observationid=str(calid), filepath=str(zip_loc), caltype=2)
            logger.info("Uploaded calibrator solutions from {} to the database".format(calid))
        except:
            logger.warn("Failed to upload calibration files")
        os.remove(zip_loc)

def filename_prefix(obsid, pulsar, bins=None, cal=None):
    """
    Creates a filename prefix depending on inputs

    Parameters
    ----------
    obsid : `int`
        The observation ID
    pulsar : `str`
        The name of the puslar
    bins : `int`
        The number of bins in the profile
    cal : `int`
        The calibrator ID

    Returns
    -------
    pref : `str`
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

def upload_file_to_db(obsid, pulsar, filepath, filetype, subbands=None, coh=True):
    """
    Uploads a file to the pulsar database

    Parameters
    ----------
    obsid : `str`
        The observation ID
    pulsar : `str`
        The name of the pulsar
    filepath : `str`
        The file path of the file to upload
    filetype : `int`
        The type of file to upload: 1: Archive, 2: Timeseries, 3: Diagnostics, 4: Calibration Solution, 5: Bestprof
    subbands : `int`
        OPTINOAL - The number of frequency subbands in an observation. If None, will calculate it.
    coh : `boolean`
        OPTINOAL - Whether this is a coherent detection or not. Default: True
    """
    if not subbands:
        subbands = get_subbands(get_common_obs_metadata(obsid))
    web_address, auth = get_db_auth_addr()
    client.detection_file_upload(web_address, auth,
                                    observationid = str(obsid),
                                    pulsar = pulsar,
                                    filetype = int(filetype),
                                    coherent = coh,
                                    subband = subbands,
                                    filepath = filepath)

def multi_upload_files(obsid, pulsar, files_dict, subbands=None, coh=True):
    """
    Uploads any number of files to the same obsid/pulsar on the database

    Parameters
    ----------
    obsid : `str`
        The observation ID
    pulsar : `str`
        The name of the pulsar
    files_dict: dict
        A dictionary with keys 1-5 representing filetpes. Each key holds a list containing the files to be uploaded.
    subbands : `int`
        OPTINOAL - The number of frequency subbands in an observation. If None, will calculate it.
    coh : `boolean`
        OPTINOAL - Whether this is a coherent detection or not. Default: True
    """
    for filetype in files_dict.keys():
        for filename in files_dict[filetype]:
            logger.info("Uploading file to databse: {}".format(filename))
            upload_file_to_db(obsid, pulsar, filename, int(filetype), subbands=subbands, coh=coh)


"""
Test set:
Run these from the vcstools/database directory

Scattered detection (crab):
python submit_to_database.py -o 1127939368 --calid 1127939368 -p J0534+2200 --bestprof tests/test_files/1127939368_J05342200.bestprof -L DEBUG
Expected flux: 7500

Weak detection (it's a 10 min detection):
python submit_to_database.py -o 1222697776 --calid 1222695592 -p J0034-0721 --bestprof ../tests/test_files/1222697776_PSR_J0034-0721.pfd.bestprof
Expected flux: 640

Medium detection:
python submit_to_database.py -o 1222697776 --calid 1222695592 -p J2330-2005 --bestprof tests/test_files/1222697776_PSR_J2330-2005.pfd.bestprof
Expected flux: ~180

Strong detection:
python submit_to_database.py -o 1226062160 --calid 1226054696 -p J2330-2005 --bestprof tests/test_files/1226062160_J2330-2005.bestprof -L DEBUG
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
    parser.add_argument('-o', '--obsid', type=str,
            help='The observation ID (eg. 1221399680).')
    parser.add_argument('-O', '--calid', type=str,
            help='The observation ID of the calibrator.')
    parser.add_argument('-p','--pulsar', type=str,
            help='The pulsar J name.')
    parser.add_argument('--incoh', action='store_true',
            help='Used for incoherent detections to accurately calculate gain. The default is coherent.')
    parser.add_argument('--andre', action='store_true',
            help="Used for calibrations done using Andre Offrina's tools. Default is RTS.")
    parser.add_argument('--sefd_file', type=str,
            help='The output file of the pabeam.py code to be used for accurate flux density calculations.')
    parser.add_argument('--simple_sefd', action='store_true',
            help="Just use the tile beam to estimate the SEFD (T_sys and gain) instead of submiting a job to do a full tied-array beam simulation. Default False")
    parser.add_argument("--vcstools_version", type=str, default="master",
            help="VCSTools version to load in jobs (i.e. on the queues) ")
    parser.add_argument("-L", "--loglvl", type=str, choices=loglevels.keys(), default="INFO",
            help="Logger verbosity level. Default: INFO")
    parser.add_argument("-V", "--version", action='store_true', help="Print version and quit")

    calcargs = parser.add_argument_group('Detection Calculation Options', 'All the values of a pulsar detection (such as flux density, width, scattering) can be calculated by this code using either a .besprof file or a soon to be implemented DSPSR equivalent and automatically uploaded to the MWA Pulsar Database. Analysis files can still be uploaded before this step but this will leave the pulsar values as nulls.')
    calcargs.add_argument('-b', '--bestprof', type=str,
            help='The location of the .bestprof file. Using this option will cause the code to calculate the needed parameters to be uploaded to the database (such as flux density, width and scattering). Using this option can be used instead of inputting the observation ID and pulsar name.')
    calcargs.add_argument('--ascii', type=str,
            help='The location of the ascii file (pulsar profile output of DSPSR). Using this option will cause the code to calculate the needed parameters to be uploaded to the database (such as flux density, width and scattering).')
    calcargs.add_argument('--dont_upload', action='store_true',
            help='Will not upload the results or files to the database. Instead will just output the results of the analysis and flux calibration to a file.')
    calcargs.add_argument('--pointing', type=str,
            help="The pointing of the detection with the RA and Dec seperated by _ in the format HH:MM:SS_+DD:MM:SS, e.g. \"19:23:48.53_-20:31:52.95 19:23:40.00_-20:31:50.00\". This is only required if the bestprof has a non standard input fits file.")
    calcargs.add_argument('--start', type=int,
            help="The start time of the detection in GPS format.")
    calcargs.add_argument('--stop', type=int,
            help="The stop time of the detection in GPS format.")
    calcargs.add_argument('--trcvr', type=str, default=data_load.TRCVR_FILE,
            help='File location of the receiver temperatures to be used. Only required if you do not want to use the default values located in %(default)s.')

    uploadargs = parser.add_argument_group('Upload Options', 'The different options for each file type that can be uploaded to the pulsar database. Will cause an error if the wrong file type is being uploaded.')
    uploadargs.add_argument('--cal_dir_to_tar', type=str,
            help='The calibration directory of a calibration solution that you would like to tar and upload to the database (eg. /group/mwavcs/vcs/1221832280/cal/1221831856/rts). Must be used with --srclist so the correct source list is uploaded. If the calibration files are in the default positions then they will be tared and uploaded.')
    uploadargs.add_argument('--srclist', type=str,
            help='Used with --cal_dir to indicate the source list file location. eg /group/mwavcs/vcs/1221832280/cal/1221831856/vis/srclist_pumav3_EoR0aegean_EoR1pietro+ForA_1221831856_patch1000.txt.')
    uploadargs.add_argument('-c', '--calibration', type=str,
            help='The calibration solution file location to be uploaded to the database. Expects a single file so please zip or tar up the bandpass calibrations, the DI Jones matrices, the flagged_channels.txt file, the flagged_tiles.txt file, the rts.in file and the source file.')
    uploadargs.add_argument('-a', '--archive', type=str,
            help="The DSPSR archive file location to be uploaded to the database. Expects a single file that is the output of DSPSR using the pulsar's ephemeris.")
    uploadargs.add_argument('--single_pulse_series', type=str,
            help='The single pulse series file location to be uploaded to the database. Expects a single file that is the output of DSPSR in single pulse mode (the -s option).')
    uploadargs.add_argument('--ppps', type=str,
            help="The Presto Prepfold PostScript file location to be uploaded to the database. Expects a single file that is the output of PRESTO's prepfold script.")
    uploadargs.add_argument('-i', '--ippd', type=str,
            help="The Intergrated Pulse Profile (sometimes called a pulse profile) file location. Expects a single file that is the output of DSPSR's pav script.")
    uploadargs.add_argument('-w', '--waterfall', type=str,
            help="The file location of a waterfall plot of pulse phase vs frequency. Expects a single file that is the output of DSPSR's psrplot.")
    args = parser.parse_args()

    # Set up the logger for stand-alone execution
    logger = setup_logger(logger, log_level=loglevels[args.loglvl])

    if args.version:
        try:
            import version
            logger.info(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            logger.error("Couldn't import version.py - have you installed vcstools?")
            logger.error("ImportError: {0}".format(ie))
            sys.exit(0)

    if not args.pulsar or not args.obsid:
        logger.error("Please supply both --obsid and --pulsar")

    # Get meta data from obsid
    common_metadata = get_common_obs_metadata(args.obsid,)
    _, ra_obs, _, time_obs, delays, centrefreq, channels = common_metadata
    # Perform frequency calculations
    minfreq = float(min(channels))
    maxfreq = float(max(channels))
    bandwidth = 30720000. #In Hz
    # Check if the data is contiguous in frequency
    if (maxfreq - minfreq) == 23:
        obstype = 1 # contiguous
    else:
        obstype = 2 # picket fence
    subbands = get_subbands(common_metadata)


    # Read in the pulsar profile -----------------------------------------------
    if args.bestprof:
        #Does the flux calculation and submits the results to the MWA pulsar database
        bestprof_data = prof_utils.get_from_bestprof(args.bestprof, pointing_input=args.pointing)
    elif args.ascii:
        if not args.stop or not args.start:
            logger.error("Please supply both start and stop times of the detection for ascii files")
            sys.exit(1)
        profile, num_bins = prof_utils.get_from_ascii(args.ascii)
        _, dm, period = get_psrcat_dm_period(args.pulsar)[0]
        time_detection = args.stop - args.start
        bestprof_data = [args.obsid, args.pulsar, dm, period, None, args.start,
                         time_detection, profile, num_bins]
        # TODO work out the best way to add pointing to ascii files


    # Perform the profile analysis and flux calibration if profile available ------------------------
    if args.bestprof or args.ascii:
        logger.info("Performing profile analaysis and flux density calculation")
        det_kwargs, sn, u_sn = analyise_and_flux_cal(args.pulsar,
                bestprof_data, args.calid,
                common_metadata=common_metadata,
                trcvr=args.trcvr,
                simple_sefd=args.simple_sefd, sefd_file=args.sefd_file,
                vcstools_version=args.vcstools_version,
                args=args)
        det_kwargs["dm"] = float(bestprof_data[2]) # DM from profile
    else:
        # No profile provided so just create pulsar table without the pulsar properties
        det_kwargs = {}
    det_kwargs["observationid"]         = str(args.obsid)
    det_kwargs["pulsar"]                = args.pulsar
    det_kwargs["subband"]               = int(subbands)
    det_kwargs["startcchan"]            = int(minfreq)
    det_kwargs["stopcchan"]             = int(maxfreq)
    det_kwargs["observation_type"]      = int(obstype)

    if args.dont_upload:
        # Just output a results file instead of uploading to database
        file_name = "{}_{}_{}_flux_results.csv".format(args.pulsar, args.obsid, args.calid)
        w = csv.writer(open(file_name, "w"))
        # loop over dictionary keys and values
        for key, val in det_kwargs.items():
            w.writerow([key, val])
        # also output the calculated SN
        w.writerow(["sn", sn])
        w.writerow(["u_sn", u_sn])
        logger.info("Outut file {} written. Exiting".format(file_name))
        sys.exit(0)


    # Create all the necessary database tables and upload the data/files ---------------------------------------
    web_address, auth = get_db_auth_addr()

    # Create pulsar table if required
    check_db_and_create_pulsar(args.pulsar)

    # Defaults for incoh and calibrator type
    if args.incoh:
        coh = False
        calibrator_type = None
        cal_db_id = None
    else:
        coh = True
        if not args.calid:
            logger.error("Please include --calid for coherent observations")
            sys.exit(1)
        if args.andre:
            calibrator_type = 1
        else:
            calibrator_type = 2
    # Create calibration table if required
    if coh:
        cal_list = client.calibrator_list(web_address, auth)
        cal_already_created = False
        for c in cal_list:
            if ( c[u'observationid'] == int(args.calid) ) and ( c[u'caltype'] == calibrator_type ):
                cal_already_created = True
                cal_db_id = c[u'id']
        if not cal_already_created:
            cal_db_id = client.calibrator_create(web_address, auth,
                                                 observationid = str(args.calid),
                                                 caltype = calibrator_type)[u'id']
    det_kwargs["coherent"] = coh
    det_kwargs["calibrator"] = cal_db_id

    # Create detection table
    try:
        # Try to create a new detection
        client.detection_create(web_address, auth, **det_kwargs)
        logger.info("Detection uploaded to database")
    except:
        client.detection_update(web_address, auth, **det_kwargs)
        logger.info("Detection already on database so updating the values")


    # Upload analysis files to the database --------------------------------------------------------

    if args.bestprof or args.archive or args.single_pulse_series or \
       args.ppps or args.ippd or args.waterfall:
        # Create filname prefix
        bins = bestprof_data[8]
        fname_pref = filename_prefix(args.obsid, args.pulsar, bins=bins, cal=args.calid)
        upfiles_dict={"1":[], "2":[], "3":[], "4":[], "5":[]}
        remove_list=[]

    # Change file names to standard format
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

    if args.bestprof or args.archive or args.single_pulse_series or \
       args.ppps or args.ippd or args.waterfall:
        # Upload all files and remove copies
        multi_upload_files(str(args.obsid), args.pulsar, upfiles_dict, subbands=subbands, coh=coh)
        for filename in remove_list:
            os.remove(filename)

    # Upload calibration solution
    if args.cal_dir_to_tar:
        if not args.srclist:
            logger.error("You must use --srclist to define the srclist file location. Exiting")
            sys.exit(0)
        # Zip up the calibration files and return the directory of the zipped file
        args.calibration = zip_calibration_files(args.cal_dir_to_tar, args.calid, args.srclist)
    if args.calibration:
        logger.info("Uploading calibration solution to database")
        if args.andre:
            cp(str(args.calibration),str(args.calid) + "_andre_calibrator.bin")
            cal_file_loc = str(args.calid) + "_andre_calibrator.bin"
        else:
            cp(str(args.calibration),str(args.calid) + "_rts_calibrator.tar")
            cal_file_loc = str(args.calid) + "_rts_calibrator.tar"

        upload_cal_files(str(args.obsid), str(args.calid), cal_file_loc, args.srclist, caltype=calibrator_type)
        os.remove(cal_file_loc)