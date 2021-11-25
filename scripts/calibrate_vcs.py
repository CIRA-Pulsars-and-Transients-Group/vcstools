#!/usr/bin/env python3

"""Script for RTS calibration jobs in the VCS Pulsar pipeline.

The class defined in this script attempts to abstract away any of the underlying
difficulties in creating a RTS configuration file. It can handle both types of VCS observing:
    "contiguous bandwidth" observations (i.e. all 24 coarse channels are adjacent), and
    "picket-fence" observations (i.e. where the 24 coarse channels are arbitrarily distributed

Author: Bradley Meyers
Date: 18-December-2017 (v0.9)
"""

import os
import sys
import numpy as np
import re
import argparse
from itertools import groupby
from operator import itemgetter
import glob
import logging
from time import strptime, strftime
import socket

from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.io import fits

from vcstools.job_submit import submit_slurm
from vcstools.general_utils import mdir, setup_logger
from vcstools.metadb_utils import getmeta
from vcstools.config import load_config_file

logger = logging.getLogger(__name__)


class CalibrationError(Exception):
    pass


class BaseRTSconfig(object):
    """A class to hold the base information for a given RTS configuration file.

    Parameters
    ----------
    obsid  : `int`
        The target observation ID.
    cal_obsid  : `int`
        The calibrator observation ID.
    metafits  : `str`
        Path to the metafits file for the calibrator observation.
    srclist  : `str`
        Path to the source list created for this calibrator observation (using srclist_by_beam.py).
    datadir  : `str`
        Path to visibility data to be used to produce a calibration solution.
    outdir  : `str`, optional
        Path to write RTS configuration and relevant flagging information.
    offline  : `str`, optional
        Whether the visibility data were produced with the offline correlator. Default is False.
    beam_model : `str`, optional
        Which beam model to use for calibrating. Can be either 'FEE2016' (default) or 'ANALYTIC'
    vcstools_version : `str`, optional
        Which version of vcstools to load - only used when correlating offline. Default is 'master'


    Attributes
    ----------
    obsid  : `int`
        The target observation ID
    cal_obsid  : `int`
        The calibrator observation ID
    offline : bool
        Whether the calibrator data was produced by the offline correlator.
    utctime  : `str`
        The start UTC time (as YYYY-MM-DDThh:mm:ss.ss)
    nfine_chan  : `int`
        The number of fine channels per coarse channel
    channels : list of ints
        The list of coarse channels recorded within a subband. Can range from length of 1-24.
    fine_cbw  : `float`
        The fine channel bandwidth in MHz
    max_frequnecy  : `float`
        The maximum frequency in the provided data, used by the RTS to calcualte the decorrelation across the band
    corr_dump_time  : `float`
        The time scale on which the correlator dumped the visibilities (i.e. the integration time).
    n_corr_dumps_to_average  : `int`
        The number of correlator dumps to use.
        Must be such that `corr_dump_time * n_corr_dumps_to_average` is <= than the total amount of data
        available for the calibrator observation.
    n_integration_bins  : `int`
        The number of visbility groups to construct (and then integrate over). 5 seems to be ok - but this will change with configuration.
    PB_HA  : `float`
        The primary beam Hour Angle (in degrees)
    PB_DEC  : `float`
        The primary beam Declination (in degrees)
    freq_base :  float
        The starting frequency from which to count the subband IDs (coarse channels)
    JD  : `float`
        Julian day conversion of `utctime`
    metafits_RTSform  : `str`
        A modified string of the user-define metafits file location. Truncates "_metafits_pdds.fits".
    ArrayPositionLat  : `float`
        The MWA's latitude (in degrees)
    ArrayPositionLong  : `float`
        The MWA's longitude (in degrees)
    base_str  : `str`
        The basic RTS configuration skeleton script with some relevant information (built in BaseRTSconfig)
    data_dir  : `str`
        Path to look for visibilities to use for calibration.
    output_dir  : `str`
        Path to write RTS configuration and relevant flagging information.
    batch_dir  : `str`
        Path to where the SLURM scripts and their output are to be written.
    metafits  : `str`
        Path to the original metafits file for the calibrator observation.
    source_list : `str`
        Path to the source list created for this calibrator observation (using srclist_by_beam.py).
    useCorrInput  : `int`
        Option value for RTS configuration regarding interpreting correlator streamed data.
        For offline-correlated data, `useCorrInput=1` and `readDirect=0`.
    readDirect  : `int`
        Option value for RTS configuration regarding reading data files from disk.
        For online-correlated data, `readDirect=1` and `useCorrInput=0`.

    Raises
    ------
    CalibrationError
        When there is a problem with some of the observation information and/or its manipulation.
    """

    def __init__(self, obsid, cal_obsid, metafits, srclist, n_int_bins=6, datadir=None, outdir=None, offline=False, beam_model="FEE2016", vcstools_version="master"):
        self.obsid = obsid  # target observation ID
        self.cal_obsid = cal_obsid  # calibrator observation ID
        self.offline = offline  # switch to decide if offline correlated data or not
        self.utctime = None  # start UTC time
        self.nfine_chan = None  # number of fine channels
        self.channels = None  # actual channel numbers
        self.fine_cbw = None  # fine channel bandwidth
        self.max_frequency = None # the maximum frequency used by the RTS to calculate decorrelation
        self.corr_dump_time = None  # correlator dump times (i.e. integration time)
        self.n_dumps_to_average = None  # number of integration times to use for calibration
        self.PB_HA = None  # primary beam HA
        self.PB_DEC = None  # primary beam Dec
        self.freq_base = None  # frequency base for RTS
        self.JD = None  # time base for RTS
        self.metafits_RTSform = None  # modified metafits file name for RTS
        self.ArrayPositionLat = -26.70331940  # MWA latitude
        self.ArrayPositionLong = 116.6708152  # MWA longitude
        self.n_integration_bins = n_int_bins # number of visibility integration groups for RTS
        self.base_str = None  # base string to be written to file, will be editted by RTScal
        self.beam_model = beam_model # The beam model to use for the calibration solutions. Either 'ANALYTIC' or 'FEE2016'
        self.beam_model_bool = None
        self.vcstools_version = vcstools_version

        comp_config = load_config_file()

        # Check to make sure paths and files exist:
        # First, check that the actual data directory exists
        if datadir is None:
            # use the default data path
            self.data_dir = os.path.join(comp_config['base_data_dir'], str(obsid), "cal", str(cal_obsid), "vis")
            logger.info("Using default calibrator data path: {0}".format(self.data_dir))
            if os.path.exists(os.path.realpath(self.data_dir)) is False:
                errmsg = "Default data directory ({0}) does not exist. Aborting.".format(self.data_dir)
                logger.error(errmsg)
                raise CalibrationError(errmsg)
        elif os.path.isdir(datadir):
            self.data_dir = os.path.realpath(datadir)
            logger.info("Using the user specified data directory: {0}".format(datadir))
        else:
            errmsg = "Data directory ({0}) does not exist. Aborting.".format(datadir)
            logger.error(errmsg)
            raise CalibrationError(errmsg)

        # Then check if the specified output and batch directories exists
        if outdir is None:
            # this is the default
            logger.info("Assuming default directory structure...")
            self.output_dir = os.path.join(comp_config['base_data_dir'], str(self.obsid), "cal", str(self.cal_obsid), "rts")
            self.batch_dir =os.path.join(comp_config['base_data_dir'], str(self.obsid), "batch")
            logger.debug("RTS output directory is {0}".format(self.output_dir))
            logger.debug("Batch directory is {0}".format(self.batch_dir))
            mdir(self.output_dir, "RTS", gid=comp_config['gid'])
            mdir(self.batch_dir, "Batch", gid=comp_config['gid'])
        else:
            # mdir handles if the directory already exists
            self.output_dir = os.path.realpath(outdir + "/rts")
            self.batch_dir = os.path.realpath(outdir + "/batch")
            logger.warning("Non-standard RTS output path: {0}".format(self.output_dir))
            logger.warning("Non-standard batch directory path: {0}".format(self.batch_dir))
            mdir(self.output_dir, "RTS", gid=comp_config['gid'])
            mdir(self.batch_dir, "Batch", gid=comp_config['gid'])

        # Then check that the metafits file exists
        if os.path.isfile(metafits) is False:
            # file doesn't exist
            errmsg = "Given metafits file ({0}) does not exist.".format(metafits)
            logger.error(errmsg)
            raise CalibrationError(errmsg)
        elif "_ppds" not in metafits:
            # file doesn't have the correct naming convention
            errmsg = "Looks like you have an old-style metafits. " \
                     "You'll need to download the new version, which is named like: " \
                     "{0}_metafits_ppds.fits.".format(obsid)
            logger.error(errmsg)
            raise CalibrationError(errmsg)
        else:
            logger.info("Metafits file exists and is named correctly.")
            self.metafits = os.path.realpath(metafits)
            logger.debug("    {0}".format(self.metafits))

        # the check that the source list exists
        if os.path.isfile(srclist) is False:
            # file doesn't exist
            errmsg = "Given source list file ({0}) does not exist.".format(srclist)
            logger.error(errmsg)
            raise CalibrationError(errmsg)
        else:
            logger.info("Checking source list file exists... Ok")
            self.source_list = os.path.realpath(srclist)

        # Check the 'beam_model' is one of the correct choices
        choices = ("FEE2016", "ANALYTIC")
        if self.beam_model not in choices:
            errmsg = "Given beam model: {0} not an available choice: {1}".format(self.beam_model, choices)
            logger.error(errmsg)
            raise CalibrationError(errmsg)
        else:
            logger.info("Using {0} beam model for calibration solution".format(self.beam_model))
            self.beam_model_bool=int(bool(self.beam_model == "ANALYTIC")) #produces 1 for ANALYTIC, 0 for FEE2016

        # set some RTS flags based on if we have offline correlated data or not
        logger.info("Setting RTS data input flags...")
        if self.offline:
            self.useCorrInput = 1
            self.readDirect = 0
            logger.debug("Offline correlation")
        else:
            self.useCorrInput = 0
            self.readDirect = 1
            logger.debug("Online correlation")


    def power_of_2_less_than(self, n):
        """Return the largest power of 2 that is less than or equal to n"""
        # using the magic of Python integers actually being objects
        return 2**(int(n).bit_length() - 1)


    def get_info_from_data_header(self):
        """Read information from the FITS file header to figure out calibration configuration.

        Raises
        ------
        CalibrationError
            When there is a problem with some of the observation information and/or its manipulation.
        """
        # first determine the UTC time from the file name
        logger.info("Gathering information from data headers...")
        file_glob = "{0}/*_gpubox*.fits".format(self.data_dir)
        files = sorted(glob.glob(file_glob))
        len_files = len(files)
        if len_files == 0:
            errmsg = "No *_gpubox*.fits files found in {0}.".format(self.data_dir)
            logger.error(errmsg)
            raise CalibrationError(errmsg)
        elif len_files % 24 != 0:
            errmsg = "Number of *_gpubox*.fits files is not divisible by 24!?"
            logger.critical(errmsg)
            raise CalibrationError(errmsg)

        first_file = files[0]
        self.utctime = os.path.splitext(os.path.basename(first_file))[0].split("_")[1]

        if self.offline is False:
            # now figure out how much data we have in total by counting the number of data HDUs
            # then open the header to access the frequency spacing and inetgration times
            hdulist = fits.open(first_file)
            nfine_chan = int(hdulist[1].header['NAXIS2'])
            inttime = float(hdulist[1].header['INTTIME'])

            # coarse channel BW is 1.28 MHz
            self.fine_cbw = 1.28 / nfine_chan
            self.nfine_chan = nfine_chan

            # the correlator dump time is whatever is in the header
            self.corr_dump_time = inttime

            ngroups = len_files // 24
            fgrouped = np.array(np.array_split(files, ngroups))
            ndumps = 0
            for item in fgrouped[:,0]:
                # count how many units are present, subtract one (primary HDU)
                ndumps += len(fits.open(item)) - 1

            # use all of the data to calibrate
            # TODO: the RTS currently (14 Nov 2017) dies when using more than 1 set of visibilities
            # Currently require that the number of dumps to average is a power of 2 <= calculated ndumps here (27 Feb 2018)
            self.n_dumps_to_average = self.power_of_2_less_than(ndumps)

        else:
            # we have to figure out how much data has been correlated and what the frequency/time resolution is
            # offline correlated data don't have a primary HDU, so we get the info from the first header
            hdulist = fits.open(first_file)
            nfine_chan = int(hdulist[0].header['NAXIS2'])
            inttime = float(hdulist[0].header['INTTIME'])

            self.fine_cbw = 1.28 / nfine_chan
            self.nfine_chan = nfine_chan
            self.corr_dump_time = inttime

            # for offline correlation, each file is one integration time - there has been no concatenation
            # TODO: this assumes that the offline correlator ALWAYS produces 1 second FITS files
            ndumps = len(fits.open(first_file)) * (len_files // 24)
            self.n_dumps_to_average = self.power_of_2_less_than(ndumps)

        logger.info("Number of fine channels: {0}".format(self.nfine_chan))
        logger.info("Fine channel bandwidth (MHz): {0}".format(self.fine_cbw))
        logger.info("Integration time (s): {0}".format(self.corr_dump_time))
        logger.info("Number of correlator dumps to average: {0}".format(self.n_dumps_to_average))


    def construct_base_string(self):
        """Construct the basic string to be written to the RTS config file.
        This string will then be edit with regexp to update the relevant details.

        Returns
        -------
        file_str  : `str`
            The base string to be written to an RTS configuration srcipt after manipulation.

        Raises
        ------
        CalibrationError
            When there is a problem with some of the observation information and/or its manipulation.
        """
        # get calibrator observation information from database
        logger.info("Querying metadata database for obsevation information...")
        obsinfo = getmeta(service='obs', params={'obs_id': str(self.cal_obsid)})

        # quick check to make sure what's returned is actually real data
        if obsinfo[u'metadata'] is None:
            errmsg = "Metadata database error (metadata empty). Maybe an invalid obs ID?"
            logger.error(errmsg)
            raise CalibrationError(errmsg)

        # get the RA/Dec pointing for the primary beam
        ra_pointing_degs = obsinfo['metadata']['ra_pointing']
        dec_pointing_degs = obsinfo['metadata']['dec_pointing']

        # now get the absolute channels
        self.channels = obsinfo[u'rfstreams'][u"0"][u'frequencies']
        # and figure out the MaxFrequency parameter
        self.max_frequency = 1.28 * (max(self.channels) + 1) # this way we ensure that MaxFrequency is applicable for ALL subbands

        # use the same operations as in timeconvert.py for our specific need
        logger.info("Converting times with astropy")
        #mwa_loc = EarthLocation.of_site('Murchison Widefield Array')
        mwa_loc = EarthLocation.from_geodetic(lon="116:40:14.93", lat="-26:42:11.95", height=377.8)
        #Astropy formating
        utctime = strptime(self.utctime, '%Y%m%d%H%M%S')
        a_time = strftime('%Y-%m-%dT%H:%M:%S', utctime)
        obstime = Time(a_time, format='fits', scale='utc', location=mwa_loc)
        lst_in_hours = obstime.sidereal_time('apparent').hourangle
        jd = obstime.jd
        logger.info("   LST: {0}".format(lst_in_hours))
        logger.info("   JD : {0}".format(jd))

        # set the HA of the image centre to the primary beam HA
        logger.debug("Determining HA and DEC for primary beam")
        self.JD = jd
        pb_ha_hours = (ra_pointing_degs / 360.0) * 24.0
        self.PB_HA = lst_in_hours - pb_ha_hours
        self.PB_DEC = dec_pointing_degs

        logger.debug("Primary beam: HA = {0} hrs, Dec = {1} deg".format(self.PB_HA, self.PB_DEC))
        logger.debug("JD = {0}".format(self.JD))

        # get the lowest frequency channel
        freqs = obsinfo['rfstreams']['0']['frequencies']
        start_channel = freqs[0]

        self.freq_base = start_channel * 1.28e6 - 0.64e6 + 15e3  # frequency base in Hz (based on Steve Ord's logic)
        self.freq_base = self.freq_base / 1.0e6  # convert to MHz

        logger.debug("Frequency lower edge = {0} MHz".format(self.freq_base))

        # make metafits file formatted for RTS
        self.metafits_RTSform = self.metafits.split("_metafits_ppds.fits")[0]
        logger.debug("RTS form metafits location: {0}".format(self.metafits_RTSform))

        # create the final file string, expanding symlinks to real paths
        logger.info("Constructing base RTS configuration script content")
        file_str = """
ReadAllFromSingleFile=
BaseFilename={base}/*_gpubox
ReadGpuboxDirect={read_direct}
UseCorrelatorInput={use_corr_input}

ReadMetafitsFile=1
MetafitsFilename={metafits_file}

DoCalibration=
doMWArxCorrections=1
doRawDataCorrections=1
doRFIflagging=0
useFastPrimaryBeamModels={beam_model_bool}
generateDIjones=1
applyDIcalibration=1
UsePacketInput=0
UseThreadedVI=1

MaxFrequency={max_freq}
ObservationFrequencyBase={base_freq}
ObservationTimeBase={base_time}
ObservationPointCentreHA={obs_ha}
ObservationPointCentreDec={obs_dec}
ChannelBandwidth={fcbw}
NumberOfChannels={nchan}

CorrDumpsPerCadence={corr_dumps_per_cadence}
CorrDumpTime={corr_dump_time}
NumberOfIntegrationBins={n_int_bins}
NumberOfIterations=1

StartProcessingAt=0

ArrayPositionLat={array_lat}
ArrayPositionLong={array_lon}
ArrayNumberOfStations=128

ArrayFile=

SourceCatalogueFile={source_list}
NumberOfCalibrators=1
NumberOfSourcesToPeel=0
calBaselineMin=20.0
calShortBaselineTaper=40.0
FieldOfViewDegrees=1""".format(base=self.data_dir,
                               read_direct=self.readDirect,
                               use_corr_input=self.useCorrInput,
                               metafits_file=self.metafits_RTSform,
                               beam_model_bool=self.beam_model_bool,
                               max_freq=self.max_frequency,
                               base_freq=self.freq_base,
                               base_time=self.JD,
                               obs_ha=self.PB_HA,
                               obs_dec=self.PB_DEC,
                               fcbw=self.fine_cbw,
                               nchan=self.nfine_chan,
                               corr_dumps_per_cadence=self.n_dumps_to_average,
                               corr_dump_time=self.corr_dump_time,
                               n_int_bins=self.n_integration_bins,
                               array_lat=self.ArrayPositionLat,
                               array_lon=self.ArrayPositionLong,
                               source_list=self.source_list)

        return file_str


    def write_flag_files(self):
        """Given the output directory, write initial flagging files based on bad tiles in metafits and number
        of fine channels.
        """
        metafits = fits.open(self.metafits)  # read metafits file
        bad_tiles = metafits[1].data['Flag'][::2]  # both polarisation recorded, so we just want every second value
        bad_tiles = np.where(bad_tiles == 1)[0]
        flagged_tiles = "{0}/flagged_tiles.txt".format(self.output_dir)

        if os.path.isfile(flagged_tiles):
            logger.warning("{0} already exists. Not overwriting.".format(flagged_tiles))
        else:
            with open(flagged_tiles, 'w') as fid:
                for b in bad_tiles:
                    fid.write("{0}\n".format(b))

        # figure out how many edge channels to flag based on the fact that with 128, we flag the edge 8
        ntoflag = int(8 * self.nfine_chan / 128.0)
        logger.debug("Flagging {0} edge channels".format(ntoflag))
        chans = np.arange(self.nfine_chan)
        start_chans = chans[:ntoflag]
        end_chans = chans[-ntoflag:]
        center_chan = [self.nfine_chan // 2]
        bad_chans = np.hstack((start_chans, center_chan, end_chans))
        flagged_channels = "{0}/flagged_channels.txt".format(self.output_dir)
        if os.path.isfile(flagged_channels):
            logger.warning("{0} already exists. Not overwriting.".format(flagged_channels))
        else:
            with open(flagged_channels, 'w') as fid:
                for b in bad_chans:
                    fid.write("{0}\n".format(b))


    def run(self):
        """Run through the pipeline to produce the RTS file string and write the channel/tile flags to disk."""
        self.get_info_from_data_header()
        self.base_str = self.construct_base_string()
        self.write_flag_files()


class RTScal(object):
    """Class to contain calibration information, create and submit RTS calibration jobs.
    It is able to distinguish picket-fence and "normal" observations, and write the
    appropriate RTS configuration files and batch scripts for each circumstance.

    The "run()" method ensures that information is gathered in order and that the
    correct methods are called at the appropriate time. It should be the only method
    call after initialising a RTScal object.


    Attributes
    ----------
    obsid  : `int`
        The target observation ID
    cal_obsid  : `int`
        The calibrator observation ID
    channels : list of ints
        The list of coarse channels recorded within a subband. Can range from length of 1-24.
    hichans : list of ints
        The coarse channel numbers >= 129, where the PFB begins to use frequencies in reverse order (i.e. channel 129 is
        higher in real-frequency than channel 130)"
    lochans : list of ints
        The coarse channel numbers < 129
    picketfence : bool
        True is the observation is picket fence (i.e. one or more subbands are not adjacent).
    rts_out_dir  : `str`
        Path to where to write the relevant files to disk (i.e. RTS configuration files and flagging information)
    batch_dir  : `str`
        Path to where the SLURM scripts and their output are to be written.
    base_str  : `str`
        The basic RTS configuration skeleton script with some relevant information (built in BaseRTSconfig)
    script_body : list of strs
        Each element is a line to be written into the SBATCH script which will run the job on the compute nodes.
    submit : bool
        Whether to submit the created scripts to the SLURM queue. True = submit, False = don't submit.
    """

    def __init__(self, rts_base, submit=True):
        """Initiliase the class attributes from the rts_base information

        Paramters
        ---------
        rts_base : BaseRTSconfig
            An instance of BaseRTSconfig which has been initialised to contain the basic observation information.
        submit  : `boolean`
            Whether to submit the created scripts to the SLURM queue. True = submit, False = don't submit.
        """
        # initilaise class attributes with user-specified options
        self.obsid = rts_base.obsid
        self.cal_obsid = rts_base.cal_obsid

        # initilaise attributes to None and then assign values later
        self.channels = None
        self.hichans = None
        self.lochans = None
        self.picket_fence = None

        # Set up the product directory according to options passed
        self.rts_out_dir = rts_base.output_dir
        self.batch_dir = rts_base.batch_dir

        # use the channels from rts_base rather than re-query the database
        self.channels = rts_base.channels

        # use the base string from RTS base as a template
        self.base_str = rts_base.base_str

        # setup the basic body for all SLURM scripts
        self.script_body = []
        self.script_body.append("cd {0}".format(self.rts_out_dir))  # change to output directory and run RTS there

        # boolean to control batch job submission. True = submit to queue, False = just write the files
        self.submit = submit

        # vcstools version to use for offline correlation
        self.offline = rts_base.offline
        self.vcstools_version = rts_base.vcstools_version


    def summary(self):
        """Print a nice looking summary of relevant attributes."""
        logger.debug("Summary of Calibration object contents:")
        logger.debug("Observation ID:             {0}".format(self.obsid))
        logger.debug("Calibrator Observation ID:  {0}".format(self.cal_obsid))
        logger.debug("RTS output directory:       {0}".format(self.rts_out_dir))
        logger.debug("Batch directory:            {0}".format(self.batch_dir))
        logger.debug("Observed absolute channels: {0}".format(self.channels))
        logger.debug("Is picket fence?            {0}".format(self.picket_fence))
        logger.debug("Submit jobs?                {0}".format(self.submit))


    def __summary(self):
        # debugging use only
        print(self.__dict__)


    def submit_true(self):
        """Set the submit flag to True, allowing sbatch queue submission."""
        logger.info("Setting submit attribute to True: allows job submission")
        self.submit = True


    def submit_false(self):
        """Set the submit flag to False, not allowing sbatch queue submission."""
        logger.info("Setting submit attribute to False: disallows job submission")
        self.submit = False


    def is_picket_fence(self):
        """Check whether the observed channels imply picket-fence or not.
        Set boolean attribute in either case.
        """
        ch_offset = self.channels[-1] - self.channels[0]
        if ch_offset == len(self.channels) - 1:
            logger.info("This observation's channels are consecutive: this is a normal observation")
            self.picket_fence = False
        else:
            logger.info("This observation's channels are NOT consecutive: this is a picket-fence observation")
            self.picket_fence = True


    def sort_obs_channels(self):
        """Just sort the channels and split them into "high" and "low" channel lists."""
        self.hichans = [c for c in self.channels if c > 128]
        self.lochans = [c for c in self.channels if c <= 128]
        logger.debug("High channels: {0}".format(self.hichans))
        logger.debug("Low channels:  {0}".format(self.lochans))


    def construct_subbands(self):
        """Group the channels into consecutive subbands, being aware of the high channel ordering reversal.
        If a subband is split across the channel 129 boundary, it is split into a "low" and "high" sub-subband.

        Returns
        -------
        hichan_groups, lochan_groups : list of lists of ints (combined size = 24)
            Each item will be a list containing 1 or more coarse channel numbers (>1 if the channels are consecutive).
        """
        logger.info("Grouping individual channels into consecutive subbands (if possible)")

        # check if the 128,129 channel pair exists and note to user that the subband will be divided on that boundary
        if any([128, 129] == self.channels[i:i + 2] for i in range(len(self.channels) - 1)):
            logger.info("A subband crosses the channel 129 boundary. This subband will be split on that boundary.")

        # create a list of lists with consecutive channels grouped within a list
        hichan_groups = [list(map(itemgetter(1), g))[::-1] for k, g in groupby(enumerate(self.hichans), lambda ix: ix[0] - ix[1])][::-1]  # reversed order (both of internal lists and globally)
        lochan_groups = [list(map(itemgetter(1), g)) for k, g in groupby(enumerate(self.lochans), lambda ix: ix[0] - ix[1])]
        logger.debug("High channels (grouped): {0}".format(hichan_groups))
        logger.debug("Low channels  (grouped): {0}".format(lochan_groups))

        return hichan_groups, lochan_groups


    def write_cont_scripts(self):
        """Function to write RTS submission script in the case of a "standard"
        contiguous bandwidth observation. This is relatively simple and just
        requires us to use the standard format RTS configuration file.

        Returns
        -------
        jobids : list of ints
            A list of all the job IDs submitted to the system compute queue.
        """
        logger.info("Writing RTS configuration script for contiguous bandwidth observation")
        fname = "{0}/rts_{1}.in".format(self.rts_out_dir, self.cal_obsid)  # output file name to write
        with open(fname, 'w') as f:
            f.write(self.base_str)

        jobids = []
        nnodes = 25  # number of required GPU nodes - 1 per coarse channels + 1 master node
        rts_batch = "RTS_{0}".format(self.cal_obsid)
        slurm_kwargs = {"chdir": "{0}".format(self.rts_out_dir),
                        "time": "2:00:00",
                        "nodes": "{0}".format(nnodes),
                        "cpus-per-gpu": "1"}
        module_list = ["RTS/master"]
        commands = list(self.script_body)  # make a copy of body to then extend
        commands.append("export UCX_MEMTYPE_CACHE=n")
        commands.append("srun --export=all -N {0} -n {0} rts_gpu {1}".format(nnodes, fname))
        hostname = socket.gethostname()
        if hostname.startswith("galaxy"):
            mem = 1024
        else:
            mem = 10240
        jobid = submit_slurm(rts_batch, commands,
                             slurm_kwargs=slurm_kwargs,
                             module_list=module_list,
                             batch_dir=self.batch_dir,
                             submit=self.submit,
                             queue='gpuq',
                             export="NONE",
                             mem=mem,
                             load_vcstools=self.offline, #load if offline
                             vcstools_version = self.vcstools_version)
        jobids.append(jobid)

        return jobids


    def get_subband_config(self, chan_groups, basepath, chan_type, count):
        """Function to make the appropriate changes to a base copy of the RTS configuration scripts.
        This will interrogate the channels and decide on the "subbandIDs" and
        "ObservationBaseFrequency" that are appropriate for each subband.

        Parameters
        ----------
        chan_groups : list of lists of ints, or list of ints
            The list of coarse channels to be used when creating the RTS configuration information.
        basepath  : `str`
            Path to where the script will eventually be written.
        chan_type  : `str`
            "low" is coarse channels are <129, "high" if coarse channels are >=129
        count  : `int`
            Counter variable that keeps track of re-ordering of channels.

        Returns
        -------
        chan_file_dict : dict
            Dictionary containing the number of nodes required for each subband and the appropriate filename.
            Key = filename, value = number of nodes required.

        count  : `int`
            Counter variable that keeps track of re-ordering of channels.

        Raises
        ------
        CalibrationError
            When there is a problem with some of the observation information and/or its manipulation.
        """
        chan_file_dict = {}  # used to keep track of which file needs how many nodes

        cc = 0
        logger.info("Looping over {0} channel groups...".format(chan_type))
        for c in chan_groups:
            if len(c) == 1:
                # single channel group, write its own RTS config file
                subid = str(count + 1)

                if chan_type == "low":
                    # offset is just how many channels along we are in the list
                    offset = cc
                elif chan_type == "high":
                    # offset is now how many subbandsIDs away from 24 we are
                    offset = 24 - int(subid)
                else:
                    errmsg = "Invalid channel group type: must be 'low' or 'high'"
                    logger.error(errmsg)
                    raise CalibrationError(errmsg)

                # the base frequency is then
                basefreq = 1.28 * (c[0] - offset) - 0.625

                # and the coarse channel centre frequency is
                cfreq = 1.28 * c[0] - 0.625

                string = self.base_str
                # find the pattern and select the entire line for replacement
                string = re.sub("ObservationFrequencyBase=.*\n", "ObservationFrequencyBase={0}\n".format(basefreq),
                                string)
                # include the SubBandIDs tag
                string = re.sub("StartProcessingAt=0\n", "StartProcessingAt=0\nSubBandIDs={0}\n\n".format(subid),
                                string)

                fname = "{0}/rts_{1}_chan{2}.in".format(basepath, self.cal_obsid, c[0])
                chan_file_dict[fname] = 1  # this particular subband consists of only 1 channel

                with open(fname, 'wb') as f:
                    f.write(string)

                logger.debug(
                    "Single channel :: (subband id, abs. chan, abs. freq) = ({0}, {1}, {2})".format(subid, c[0], cfreq))

                count += 1
                cc += 1

            elif len(c) > 1:
                # multiple consecutive channels
                subids = [str(count + i + 1) for i in range(len(c))]

                # as with single channels, compute the offset
                if chan_type == "low":
                    offset = cc
                elif chan_type == "high":
                    # for high channels, there is now multiple offsets
                    offset = np.array([24 - int(x) for x in subids])
                else:
                    errmsg = "Invalid channel group type: must be 'low' or 'high'"
                    logger.error(errmsg)
                    raise CalibrationError(errmsg)

                # basefreq is then the minmium of the calculated quantities
                basefreq = min(1.28 * (np.array(c) - offset) - 0.625)

                # multiple coarse channel centre frequencies
                cfreqs = 1.28 * np.array(c) - 0.625

                string = self.base_str
                # find the pattern and select the entire line for replacement
                string = re.sub("ObservationFrequencyBase=.*\n", "ObservationFrequencyBase={0}\n".format(basefreq),
                                string)
                # include the SubBandIDs tag
                string = re.sub("StartProcessingAt=0\n",
                                "StartProcessingAt=0\nSubBandIDs={0}\n\n".format(",".join(subids)), string)

                logger.debug("Multiple channels :: (subband ids, abs. chans, abs. freqs) = {0}".format(
                    ", ".join("({0}, {1}, {2})".format(i, j, k) for i, j, k in zip(subids, c, cfreqs))))

                if chan_type == "low":
                    fname = "{0}/rts_{1}_chan{2}-{3}.in".format(basepath, self.cal_obsid, min(c), max(c))
                elif chan_type == "high":
                    fname = "{0}/rts_{1}_chan{2}-{3}.in".format(basepath, self.cal_obsid, max(c), min(c))

                chan_file_dict[fname] = len(c)
                with open(fname, 'w') as f:
                    f.write(string)

                count += len(c)
                cc += len(c)

            else:
                errmsg = "Reached a channel group with no entries!?"
                logger.error(errmsg)
                raise CalibrationError(errmsg)

        return chan_file_dict, count


    def write_picket_fence_scripts(self):
        """Function to write RTS submission scripts in the case of a picket-fence
        observation. A significant amount of extra information needs to be
        determined and placed in the RTS configuration files. There will be:

            1 RTS configuration file per set of adjacent single channels (subband)

        The only exception to that rule is where the 129 boundary is crossed,
        in which case that subband will be split in two.

        Returns
        -------
        jobids : list of ints
            A list of all the job IDs submitted to the system compute queue.
        """
        logger.info("Sorting picket-fence channels and determining subband info...")
        self.sort_obs_channels()
        hichan_groups, lochan_groups = self.construct_subbands()

        # write out the RTS config files and keep track of the number of nodes required for each
        count = 0
        lodict, count = self.get_subband_config(lochan_groups, self.rts_out_dir, "low", count)
        hidict, count = self.get_subband_config(hichan_groups, self.rts_out_dir, "high", count)
        chan_file_dict = lodict.copy()
        chan_file_dict.update(hidict)

        # Now submit the RTS jobs
        logger.info("Writing individual subband RTS configuration scripts")

        hostname = socket.gethostname()
        if hostname.startswith("galaxy"):
            mem = 1024
        else:
            mem = 10240

        jobids = []
        for k, v in chan_file_dict.items():
            nnodes = v + 1
            chans = k.split('_')[-1].split(".")[0]
            rts_batch = "RTS_{0}_{1}".format(self.cal_obsid, chans)
            slurm_kwargs = {"chdir": "{0}".format(self.rts_out_dir),
                            "time": "2:00:00",
                            "nodes": "{0}".format(nnodes),
                            "cpus-per-gpu": "1"}
            module_list= ["RTS/master"]
            commands = list(self.script_body)  # make a copy of body to then extend
            commands.append("export UCX_MEMTYPE_CACHE=n")
            commands.append("srun --export=all -N {0} -n {0} rts_gpu {1}".format(nnodes, k))
            jobid = submit_slurm(rts_batch, commands,
                                 slurm_kwargs=slurm_kwargs,
                                 module_list=module_list,
                                 batch_dir=self.batch_dir,
                                 submit=self.submit,
                                 queue='gpuq',
                                 export="NONE",
                                 mem=mem,
                                 load_vcstools=False)
            jobids.append(jobid)

        return jobids


    def run(self):
        """Only function that needs to be called after creating the RTScal object.
        Ensures functions are called in correct order to update/evaluate attributes and
        produce the RTS channel-specific files if required.

        Returns
        -------
        jobids : list of ints
            A list of all the job IDs submitted to the system compute queue.
        """
        self.is_picket_fence()  # determine if picket-fence obs or not
        self.summary()

        if self.picket_fence:
            jobids = self.write_picket_fence_scripts()
        else:
            jobids = self.write_cont_scripts()

        logger.info("Done!")
        return jobids




if __name__ == '__main__':

    # Dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING)

    # Option parsing
    parser = argparse.ArgumentParser(
        description="Script for creating RTS configuration files and submitting relevant jobs in the "
                    "VCS pulsar pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-o", "--obsid", type=int, help="Observation ID of target", default=None, required=True)

    parser.add_argument("-O", "--cal_obsid", type=int,
                        help="Calibrator observation ID for the target observation. "
                             "Can be the same as --obs if using in-beam visibilities", default=None, required=True)

    parser.add_argument("-m", "--metafits", type=str,
                        help="Path to the relevant calibration observation metafits file", default=None, required=True)

    parser.add_argument("-s", "--srclist", type=str, help="Path to the desired source list",
                        default=None, required=True)

    parser.add_argument("--gpubox_dir", type=str,
                        help="Where the *_gpubox*.fits files are located "
                             "(ONLY USE IF YOU WANT THE NON-STANDARD LOCATIONS.)", default=None)

    parser.add_argument("--rts_output_dir", type=str,
                        help="Parent directory where you want the /rts directory and /batch directory to be created "
                             "(ONLY USE IF YOU WANT THE NON-STANDARD LOCATIONS.)", default=None)

    parser.add_argument("--beam_model", type=str, choices=["ANALYTIC", "FEE2016"], default="FEE2016",
                        help="Which beam model to use to create the calibration solution")

    parser.add_argument("--n_vis_grp", type=int,
                        help="The number of visbility groups for the RTS to construct. "
                        "This changes with the MWA configuration "
                        "(i.e. long baselines require more visibility groups to combat decorrelation)",
                        default=6)

    parser.add_argument("--offline", action='store_true',
                        help="Tell the RTS to read calibrator data in the offline correlated data format")

    parser.add_argument("--vcstools_version", type=str, default="master",
                        help="The version of vcstools to use. Only relevant for offline correlation")

    parser.add_argument("--nosubmit", action='store_false', help="Write jobs scripts but DO NOT submit to the queue")

    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO",
                        choices=loglevels.keys(), default="INFO")

    parser.add_argument("-V", "--version", action='store_true', help="Print version and quit")

    args = parser.parse_args()

    if args.version:
        try:
            import version
            print(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            print("Couldn't import version.py - have you installed vcstools?")
            print("ImportError: {0}".format(ie))
            sys.exit(0)

    # set up the logger for stand-alone execution
    logger = setup_logger(logger, log_level=loglevels[args.loglvl])

    logger.info("Creating BaseRTSconfig instance - compiling basic RTS configuration script")
    try:
        baseRTSconfig = BaseRTSconfig(args.obsid, args.cal_obsid, args.metafits, args.srclist,
                                        datadir=args.gpubox_dir, outdir=args.rts_output_dir,
                                        offline=args.offline, n_int_bins=args.n_vis_grp,
                                        beam_model=args.beam_model, vcstools_version=args.vcstools_version)
        baseRTSconfig.run()
    except CalibrationError as e:
        logger.critical("Caught CalibrationError: {0}".format(e))
        sys.exit(1)

    logger.info("Creating RTScal instance - determining specific configuration for this observation")
    try:
        calobj = RTScal(baseRTSconfig, args.nosubmit)
        jobs = calobj.run()
        logger.info("Job IDs: {0}".format(jobs))
    except CalibrationError as e:
        logger.critical("Caught CalibrationError: {0}".format(e))
        sys.exit(1)
