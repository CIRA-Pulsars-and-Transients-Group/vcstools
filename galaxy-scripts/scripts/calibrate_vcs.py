#!/usr/bin/env python

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
from astropy.io import fits
from process_vcs import submit_slurm  # need to get this moved out of process_vcs.py
from mdir import mdir
from mwa_metadb_utils import getmeta
from itertools import groupby
from operator import itemgetter
import distutils.spawn
import subprocess
import glob
import logging

logger = logging.getLogger(__name__)


class CalibrationError(Exception):
    pass


class RTScal(object):
    """Class to contain calibration information, create and submit RTS calibration jobs.
    It is able to distinguish picket-fence and "normal" observations, and write the 
    appropriate RTS configuration files and batch scripts for each circumstance.

    The "run()" method ensures that information is gathered in order and that the 
    correct methods are called at the appropriate time. It should be the only method 
    call after initialising a RTScal object.


    Attributes
    ----------
    obsid : int
        The target observation ID
    cal_obsid : int
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
    rts_out_dir : str
        Path to where to write the relevant files to disk (i.e. RTS configuration files and flagging information)
    batch_dir : str
        Path to where the SLURM scripts and their output are to be written.
    base_str : str
        The basic RTS configuration skeleton script with some relevant information (built in BaseRTSconfig)
    script_body : list of strs
        Each element is a line to be written into the SBATCH script which will run the job on the compute nodes.
    submit : bool
        Whether to submit the created scripts to the SLURM queue. True = submit, False = don't submit.
    """

    def __init__(self, rts_base, submit=False):
        """Initiliase the class attributes from the rts_base information

        Paramters
        ---------
        rts_base : BaseRTSconfig
            An instance of BaseRTSconfig which has been initialised to contain the basic observation information.
        submit : boolean
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
        self.script_body.append("module load cudatoolkit")  # need CUDA for RTS
        self.script_body.append("module load cfitsio")  # need cfitsio for reading fits files
        self.script_body.append("cd {0}".format(self.rts_out_dir))  # change to output directory and run RTS there

        # boolean to control batch job submission. True = submit to queue, False = just write the files
        self.submit = submit

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
        print self.__dict__

    def submit_true(self):
        """Set the submit flag to True, allowing sbatch queue submission."""
        logger.info("Setting submit attribute to True")
        self.submit = True

    def submit_false(self):
        """Set the submit flag to False, not allowing sbatch queue submission."""
        logger.info("Setting submit attribute to False")
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
        logger.info("Grouping individual channels into consecutive chunks (if possible)")

        # check if the 128,129 channel pair exists and note to user that the subband will be divided on that boundary
        if any([128, 129] == self.channels[i:i + 2] for i in xrange(len(self.channels) - 1)):
            logger.info("A subband crosses the channel 129 boundary. This subband will be split on that boundary.")

        # create a list of lists with consecutive channels grouped within a list
        hichan_groups = [map(itemgetter(1), g)[::-1] for k, g in
                         groupby(enumerate(self.hichans), lambda (i, x): i - x)][
                        ::-1]  # reversed order (both of internal lists and globally)
        lochan_groups = [map(itemgetter(1), g) for k, g in groupby(enumerate(self.lochans), lambda (i, x): i - x)]
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
        logger.info("Writing RTS configuration script for contiguous bandwith observation")
        fname = "{0}/rts_{1}.in".format(self.rts_out_dir, self.cal_obsid)  # output file name to write
        with open(fname, 'wb') as f:
            f.write(self.base_str)

        jobids = []
        nnodes = 25  # number of required GPU nodes - 1 per coarse channels + 1 master node
        rts_batch = "RTS_{0}".format(self.cal_obsid)
        slurm_kwargs = {"partition": "gpuq", "workdir": "{0}".format(self.rts_out_dir), "time": "00:20:00",
                        "nodes": "{0}".format(nnodes)}
        commands = list(self.script_body)  # make a copy of body to then extend
        commands.append("aprun -n {0} -N 1  rts_gpu {1}".format(nnodes, fname))
        jobid = submit_slurm(rts_batch, commands, slurm_kwargs=slurm_kwargs, batch_dir=self.batch_dir,
                             submit=self.submit)
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
        basepath : str
            Path to where the script will eventually be written.
        chan_type : str
            "low" is coarse channels are <129, "high" if coarse channels are >=129
        count : int
            Counter variable that keeps track of re-ordering of channels.

        Returns
        -------
        chan_file_dict : dict
            Dictionary containing the number of nodes required for each subband and the appropriate filename. 
            Key = filename, value = number of nodes required.
        
        count : int
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
                    errmsg = "Invalid channel group type: must be \"low\" or \"high\"."
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
                    errmsg = "Invalid channel group type: must be \"low\" or \"high\"."
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
                with open(fname, 'wb') as f:
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

        # lolengths = set([len(l) for l in lochan_groups])  # figure out the unique lengths
        # hilengths = set([len(h) for h in hichan_groups])
        # lengths = lolengths.union(hilengths)  # combine the sets

        # Now submit the RTS jobs
        logger.info("Writing individual subband RTS config scripts")

        jobids = []
        for k, v in chan_file_dict.iteritems():
            nnodes = v + 1
            chans = k.split('_')[-1].split(".")[0]
            rts_batch = "RTS_{0}_{1}".format(self.cal_obsid, chans)
            slurm_kwargs = {"partition": "gpuq", "workdir": "{0}".format(self.rts_out_dir), "time": "00:45:00",
                            "nodes": "{0}".format(nnodes)}
            commands = list(self.script_body)  # make a copy of body to then extend
            commands.append("aprun -n {0} -N 1  rts_gpu {1}".format(nnodes, k))
            jobid = submit_slurm(rts_batch, commands, slurm_kwargs=slurm_kwargs, batch_dir=self.batch_dir,
                                 submit=self.submit)
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
        # self.get_obs_channels() # fetch channels
        self.is_picket_fence()  # determine if picket-fence obs or not
        self.summary()

        if self.picket_fence:
            jobids = self.write_picket_fence_scripts()
        else:
            jobids = self.write_cont_scripts()

        logger.info("Done!")
        return jobids


class BaseRTSconfig(object):
    """A class to hold the base information for a given RTS configuration file.
    
    Parameters
    ----------
    obsid : int
        The target observation ID.
    cal_obsid : int
        The calibrator observation ID.
    metafits : str
        Path to the metafits file for the calibrator observation.
    srclist : str
        Path to the source list created for this calibrator observation (using srclist_by_beam.py).
    datadir : str
        Path to visibility data to be used to produce a calibration solution.
    outdir : str, optional
        Path to write RTS configuration and relevant flagging information.
    offline : str, optional
        Whether the visibility data were produced with the offline correlator. Default is False.


    Attributes
    ----------
    obsid : int
        The target observation ID
    cal_obsid : int
        The calibrator observation ID
    offline : bool
        Whether the calibrator data was produced by the offline correlator.
    utctime : str
        The start UTC time (as YYYY-MM-DDThh:mm:ss.ss)
    nfine_chan : int
        The number of fine channels per coarse channel
    channels : list of ints
        The list of coarse channels recorded within a subband. Can range from length of 1-24.
    fine_cbw : float
        The fine channel bandwidth in MHz
    corr_dump_time : float
        The time scale on which the correlator dumped the visibilities (i.e. the integration time).
    n_corr_dumps_to_average : int
        The number of correlator dumps to use. 
        Must be such that `corr_dump_time * n_corr_dumps_to_average` is <= than the total amount of data 
        available for the calibrator observation.
    PB_HA : float
        The primary beam Hour Angle (in degrees)
    PB_DEC : float
        The primary beam Declination (in degrees)
    freq_base :  float
        The starting frequency from which to count the subband IDs (coarse channels)
    JD : float
        Julian day conversion of `utctime`
    metafits_RTSform : str
        A modified string of the user-define metafits file location. Truncates "_metafits_pdds.fits".
    ArrayPositionLat : float
        The MWA's latitude (in degrees)
    ArrayPositionLong : float
        The MWA's longitude (in degrees)
    base_str : str
        The basic RTS configuration skeleton script with some relevant information (built in BaseRTSconfig)
    data_dir : str
        Path to look for visibilities to use for calibration.
    output_dir : str
        Path to write RTS configuration and relevant flagging information.
    batch_dir : str
        Path to where the SLURM scripts and their output are to be written.    
    metafits : str
        Path to the original metafits file for the calibrator observation.
    source_list: str 
        Path to the source list created for this calibrator observation (using srclist_by_beam.py).
    useCorrInput : int
        Option value for RTS configuration regarding interpreting correlator streamed data.
        For offline-correlated data, `useCorrInput=1` and `readDirect=0`.
    readDirect : int
        Option value for RTS configuration regarding reading data files from disk. 
        For online-correlated data, `readDirect=1` and `useCorrInput=0`.

    Raises
    ------
    CalibrationError
        When there is a problem with some of the observation information and/or its manipulation.
    """

    def __init__(self, obsid, cal_obsid, metafits, srclist, datadir, outdir=None, offline=False):
        self.obsid = obsid  # target observation ID
        self.cal_obsid = cal_obsid  # calibrator observation ID
        self.offline = offline  # switch to decide if offline correlated data or not
        self.utctime = None  # start UTC time
        self.nfine_chan = None  # number of fine channels
        self.channels = None  # actual channel numbers
        self.fine_cbw = None  # fine channel bandwidth
        self.corr_dump_time = None  # correlator dump times (i.e. integration time)
        self.n_dumps_to_average = None  # number of integration times to use for calibration
        self.PB_HA = None  # primary beam HA
        self.PB_DEC = None  # primary beam Dec
        self.freq_base = None  # frequency base for RTS
        self.JD = None  # time base for RTS
        self.metafits_RTSform = None  # modified metafits file name for RTS
        self.ArrayPositionLat = -26.70331940  # MWA latitude
        self.ArrayPositionLong = 116.6708152  # MWA longitude
        self.base_str = None  # base string to be written to file, will be editted by RTScal

        # Check to make sure paths and files exist:
        # First, check that the actual data directory exists
        if os.path.isdir(datadir):
            self.data_dir = os.path.abspath(datadir)
            logger.info("Checking data directory exists... Ok")
        else:
            errmsg = "Data directory ({0}) does not exists. Aborting.".format(datadir)
            logger.error(errmsg)
            raise CalibrationError(errmsg)

        # Then check if the specified output and batch directories exists
        if outdir is None:
            # this is the default
            logger.info("Assuming default directory structure...")
            self.output_dir = "/group/mwaops/vcs/{0}/cal/{1}/rts".format(self.obsid, self.cal_obsid)
            self.batch_dir = "/group/mwaops/vcs/{0}/batch".format(self.obsid)
            logger.debug("RTS output directory is {0}".format(self.output_dir))
            logger.debug("Batch directory is {0}".format(self.batch_dir))
            mdir(self.output_dir, "RTS")
            mdir(self.batch_dir, "Batch")
        else:
            # mdir handles if the directory already exists
            self.output_dir = os.path.abspath(outdir + "/rts")
            self.batch_dir = os.path.abspath(outdir + "/batch")
            logger.warning("Non-standard RTS output path: {0}".format(self.output_dir))
            logger.warning("Non-standard batch directory path: {0}".format(self.batch_dir))
            mdir(self.output_dir, "RTS")
            mdir(self.batch_dir, "Batch")

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
            logger.info("Checking metafits file exists and is named correctly... Ok")
            self.metafits = os.path.abspath(metafits)

        # the check that the source list exists
        if os.path.isfile(srclist) is False:
            # file doesn't exist
            errmsg = "Given source list file ({0}) does not exist.".format(srclist)
            logger.error(errmsg)
            raise CalibrationError(errmsg)
        else:
            logger.info("Checking source list file exists... Ok")
            self.source_list = os.path.abspath(srclist)

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

            ngroups = len_files / 24
            fgrouped = np.array(np.array_split(files, ngroups))
            ndumps = 0
            for item in fgrouped[:, 0]:
                # count how many units are present, subtract one (primary HDU)
                ndumps += len(fits.open(item)) - 1

            # use all of the data to calibrate
            # TODO: the RTS currently (14 Nov 2017) dies when using more than 1 set of visibilities
            self.n_dumps_to_average = ndumps

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
            ndumps = len(fits.open(first_file)) * len_files / 24
            self.n_dumps_to_average = ndumps

        logger.info("Number of fine channels: {0}".format(self.nfine_chan))
        logger.info("Fine channel bandwidth (MHz): {0}".format(self.fine_cbw))
        logger.info("Integration time (s): {0}".format(self.corr_dump_time))
        logger.info("Number of correlator dumps to average: {0}".format(self.n_dumps_to_average))

    def construct_base_string(self):
        """Construct the basic string to be written to the RTS config file. 
        This string will then be edit with regexp to update the relevant details.

        Returns
        -------
        file_str : str
            The base string to be written to an RTS configuration srcipt after manipulation.

        Raises
        ------
        CalibrationError
            When there is a problem with some of the observation information and/or its manipulation.
        """
        # get calibrator observation information from database
        # TODO: need to make this write a metafile so that we don't have to keep querying the database on every run
        logger.info("Querying metadata database for obsevation information...")
        obsinfo = getmeta(service='obs', params={'obs_id': str(self.cal_obsid)})

        # quick check to make sure what's returned is actually real data
        if len(obsinfo[u'logs']) == 0:
            errmsg = "Metadata database error (logs empty). Maybe an invalid obs ID?"
            logger.error(errmsg)
            raise CalibrationError(errmsg)

        # get the RA/Dec pointing for the primary beam
        ra_pointing_degs = obsinfo['metadata']['ra_pointing']
        dec_pointing_degs = obsinfo['metadata']['dec_pointing']

        # now get teh absolute channels
        self.channels = obsinfo[u'rfstreams'][u"0"][u'frequencies']

        # convert times using our timeconvert and get LST and JD 
        # TODO: need to make this not have to call an external script
        try:
            timeconvert = distutils.spawn.find_executable("timeconvert.py")
        except Exception:
            errmsg = "Unable to access or find the executebale \"timeconvert.py\""
            logger.error(errmsg)
            raise CalibrationError(errmsg)

        cmd = "{0} --datetime={1}".format(timeconvert, self.utctime)
        logger.debug(cmd)
        time_cmd = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

        hh = 0
        mm = 0
        ss = 0
        jd = 2451545.25  # 2001-01-01 18:00:00 UTC

        for line in time_cmd.stdout:
            if "LST" in line:
                lststring, lstflag = line.split()
                hh, mm, ss = lststring.split(":")

            if "JD" in line:
                jdflag, jd = line.split()

        lst_in_hours = float(hh) + float(mm) / 60.0 + float(ss) / 60.0 ** 2

        # set the HA of the image centre to the primary beam HA
        logger.debug("Determining HA and DEC for primary beam")
        self.JD = jd
        pb_ha_hours = (ra_pointing_degs / 360.0) * 24.0
        self.PB_HA = lst_in_hours - pb_ha_hours
        self.PB_DEC = dec_pointing_degs

        logger.debug("Primary beam: HA = {0} hrs, Dec = {1} deg".format(self.PB_HA, self.PB_DEC))
        logger.debug("JD = {0}".format(self.JD))
        logger.debug("LST = {0}:{1}:{2} = {3} hrs".format(hh, mm, ss, lst_in_hours))

        # get the lowest frequency channel
        freqs = obsinfo['rfstreams']['0']['frequencies']
        start_channel = freqs[0]

        self.freq_base = start_channel * 1.28e6 - 0.64e6 + 15e3  # frequency base in Hz
        self.freq_base /= 1e6  # convert to MHz

        logger.debug("Frequency lower edge = {0} MHz".format(self.freq_base))

        # make metafits file formatted for RTS
        self.metafits_RTSform = self.metafits.split("_metafits_ppds.fits")[0]

        logger.info("Constructing base RTS configuration script content")
        file_str = """
ReadAllFromSingleFile=
BaseFilename={0}/*_gpubox
ReadGpuboxDirect={1}
UseCorrelatorInput={2}

ReadMetafitsFile=1
MetafitsFilename={3}

DoCalibration=
doMWArxCorrections=1
doRawDataCorrections=1
doRFIflagging=0
useFastPrimaryBeamModels=1
generateDIjones=1
applyDIcalibration=1
UsePacketInput=0
UseThreadedVI=1

ObservationFrequencyBase={4}
ObservationTimeBase={5}
ObservationPointCentreHA={6}
ObservationPointCentreDec={7}
ChannelBandwidth={8}
NumberOfChannels={9}

CorrDumpsPerCadence={10}
CorrDumpTime={11}
NumberOfIntegrationBins=3
NumberOfIterations=1

StartProcessingAt=0

ArrayPositionLat={12}
ArrayPositionLong={13}
ArrayNumberOfStations=128

ArrayFile=

SourceCatalogueFile={14}
NumberOfCalibrators=1
NumberOfSourcesToPeel=0
calBaselineMin=20.0
calShortBaselineTaper=40.0
FieldOfViewDegrees=1""".format(os.path.realpath(self.data_dir),
                               self.readDirect,
                               self.useCorrInput,
                               self.metafits_RTSform,
                               self.freq_base,
                               self.JD,
                               self.PB_HA,
                               self.PB_DEC,
                               self.fine_cbw,
                               self.nfine_chan,
                               self.n_dumps_to_average,
                               self.corr_dump_time,
                               self.ArrayPositionLat,
                               self.ArrayPositionLong,
                               self.source_list)

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
        ntoflag = int(8 * self.nfine_chan / 128.)
        logger.debug("Flagging {0} edge channels".format(ntoflag))
        chans = np.arange(self.nfine_chan)
        start_chans = chans[:ntoflag]
        end_chans = chans[-ntoflag:]
        center_chan = [self.nfine_chan / 2]
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


if __name__ == '__main__':
    # Option parsing
    parser = argparse.ArgumentParser(
        description="Calibration script for creating RTS configuration files in the VCS pulsar pipeline")

    parser.add_argument("-o", "--obsid", type=int, help="Observation ID of target.", required=True)

    parser.add_argument("-O", "--cal_obsid", type=int,
                        help="Calibrator observation ID for the target observation. "
                             "Can be the same as --obs if using in-beam visibilities.", required=True)

    parser.add_argument("-m", "--metafits", type=str,
                        help="Path to the relevant calibration observation metafits file.", required=True)

    parser.add_argument("-s", "--srclist", type=str, help="Path to the desired source list.", required=True)

    parser.add_argument("--gpubox_dir", type=str, help="Where the *_gpubox*.fits files are located")

    parser.add_argument("--rts_output_dir", type=str,
                        help="Parent directory where you want the /rts directory and /batch directory to be created."
                             " ONLY USE IF YOU WANT THE NON-STANDARD LOCATIONS.", default=None)

    parser.add_argument("--offline", action='store_true',
                        help="Tell the RTS to read calibrator data in the offline correlated data format.")

    parser.add_argument("--submit", action='store_true', help="Switch to allow SLURM job submission")

    args = parser.parse_args()

    # set up the logger for stand-alone execution
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s %(thread)d  %(name)s  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logger.info("Creating BaseRTSconfig instance - setting up basic information for RTS configuration scripts")
    try:
        baseRTSconfig = BaseRTSconfig(args.obsid, args.cal_obsid, args.metafits, args.srclist, args.gpubox_dir,
                                      args.rts_output_dir, args.offline)
        baseRTSconfig.run()
    except CalibrationError as e:
        logger.critical("Caught CalibrationError: {0}".format(e))
        sys.exit(1)

    logger.info("Creating RTScal instance - determining specific config requirements for this observation")
    try:
        calobj = RTScal(baseRTSconfig, args.submit)
        jobs = calobj.run()
        logger.info("Job IDs: {0}".format(jobs))
    except CalibrationError as e:
        logger.critical("Caught CalibrationError: {0}".format(e))
        sys.exit(1)
