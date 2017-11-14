#!/usr/bin/env python

"""
Script for RTS calibration jobs in the VCS Pulsar pipeline.

The class defined in this script attempts to abstract away any of the underlying 
intricicais in creating a RTS configuration file. It can handle both:
    "contiguous bandwidth" observations (i.e. all 24 coarse channels are adjacent), and
    "picket-fence" observations (i.e. where the 24 coarse channels are arbitrarily distributed

Author: Bradley Meyers
Date: 17-October-2017 (v0.1)
"""

import os
import sys
import numpy as np
import re
import argparse
from astropy.io import fits
from process_vcs import submit_slurm
from mdir import mdir
from mwa_metadb_utils import getmeta
from itertools import groupby
from operator import itemgetter
import distutils.spawn
import subprocess
import glob

class RTScal(object):
    """ 
    Class to contain calibration information, create and submit RTS calibration jobs.
    It is able to distinguish picket-fence and "normal" observations, and write the 
    appropriate RTS configuration files and batch scripts for each circumstance.

    The "run()" method ensures that information is gathered in order and that the 
    correct methods are called at the appropriate time. It should be the only method 
    call after initialising a RTScal object.
    """

    def __init__(self, obsid=None, cal_obsid=None, rts_in_file=None, rts_out_dir=None, submit=True):
        # initilaise class attributes with user-specified options
        self.obsid     = obsid
        self.cal_obsid = cal_obsid

        # check that provided RTS config script exists, else abort
        if os.path.isfile(rts_in_file):
            self.rts_in_file = rts_in_file
        else:
            print "!!ERROR!! The base RTS configuration script you provided does not exist! Aborting."
            sys.exit(1)
            
        # initilaise attributes to None and then assign values later
        self.channels     = None
        self.hichans      = None
        self.lochans      = None
        self.picket_fence = None
        
        # Set up the product directory according to options passed
        if rts_out_dir is not None:
            # non-standard output path
            print "!!WARNING!! You are not using the standard directory structure. Be sure you know what you're doing..."
            self.rts_out_dir = rts_out_dir
            self.batch_dir   = "{0}/batch".format(self.rts_out_dir)
            mdir(self.batch_dir, "Batch")
        else:
            # use the default locations:
            self.rts_out_dir = "/group/mwaops/vcs/{0}/cal/{1}".format(self.obsid, self.cal_obsid)
            self.batch_dir   = "/group/mwaops/vcs/{0}/batch".format(self.obsid)
            mdir(self.batch_dir, "Batch") # should have already been created, but just in case...

        # re-define the rts_out_dir attribute to be the directory where the solutions are written
        self.rts_out_dir = "{0}/rts".format(self.rts_out_dir)
        mdir(self.rts_out_dir, "RTS") 

        # setup the basic body for all SLURM scripts
        self.script_body = []
        self.script_body.append("module load cudatoolkit") # need CUDA for RTS
        self.script_body.append("module load cfitsio")     # need cfitsio for reading fits files
        self.script_body.append("cd {0}".format(self.rts_out_dir)) # change to output directory and run RTS there
    
        # boolean to control batch job submission. True = submit to queue, False = just write the files
        self.submit = submit


    def summary(self):
        """
        Print a nice looking summary of relevant attributes
        """
        print
        print "Summary of Calibration object contents:"
        print "\tObservation ID:             ",self.obsid
        print "\tCalibrator Observation ID:  ",self.cal_obsid
        print "\tBase RTS configuration file:",self.rts_in_file
        print "\tRTS output directory:       ",self.rts_out_dir
        print "\tBatch directory:            ",self.batch_dir
        print "\tObserved absolute channels: ",self.channels
        print "\tIs picket fence?            ",self.picket_fence
        print "\tSubmit jobs?                ",self.submit
        print

    def __summary(self):
        # debugging use only
        print self.__dict__


    def submit_True(self):
        # set the submit attribute to True (allow sbatch submission)
        self.submit = True


    def is_picket_fence(self):
        """
        Check whether the observed channels imply picket-fence or not. 
        Set boolean attribute in either case.
        """
        ch_offset = self.channels[-1] - self.channels[0] 
        if ch_offset == len(self.channels)-1:
            print "The channels are consecutive: this is a normal observation"
            self.picket_fence = False
        else:
            print "The channels are NOT consecutive: this is a picket-fence observation"
            self.picket_fence = True


    def get_obs_channels(self):
        """
        For this observation ID, fetch the channels from the database or the metafile.
        If the metafile does not exist in the output directory, then a database query 
        is sent and a metafile is written for easier/faster access in future.
        """
        metafile = "{0}/{1}.meta".format(self.rts_out_dir, self.cal_obsid)
        metafile_exists = False

        channels = None
        if os.path.isfile(metafile):
            print "Found observation metafile: {0}".format(metafile)
            with open(metafile,'rb') as m:
                for line in m.readlines():
                    if line.startswith("channels"):
                        channels = line.strip().split(",")[1:]

        if channels == None:
            print "Channels keyword not found in metafile or metafile doesn't exits."
        else:
            metafile_exists = True
            channels = [int(c) for c in channels]

        if metafile_exists == False:
            print "Querying the database for calibrator obs ID {0}...".format(self.cal_obsid)
            obs_info = getmeta(service='obs', params={'obs_id':str(self.cal_obsid)})
            channels = obs_info[u'rfstreams'][u"0"][u'frequencies']
        
        with open(metafile,"wb") as m:
            m.write("#Metadata for obs ID {0} required to determine if: normal or picket-fence\n".format(self.cal_obsid))
            m.write("channels,{0}".format(",".join([str(c) for c in channels])))

        self.channels = channels

    
    def sort_obs_channels(self):
        """
        Just sort the channels and split them into "high" and "low" channel lists
        """
        self.hichans = [c for c in self.channels if c>128]
        self.lochans = [c for c in self.channels if c<=128]
        print "High channels:",self.hichans
        print "Low channels:",self.lochans


    def construct_subbands(self):
        """
        Group the channels into consecutive subbands, being aware of the high channel ordering reversal.
        If a subband is split across the channel 129 boundary, it is split into a "low" and "high" sub-subband.
        """
        print "Grouping individual channels into consecutive chunks (if possible)"

        # check if the 128,129 channel pair exists and note to user that the subband will be divided on that boundary
        if any([128,129] == self.channels[i:i+2] for i in xrange(len(self.channels)-1)):
            print "\tNOTE: a subband crosses the channel 129 boundary. This subband will be split on that boundary."

        # create a list of lists with consecutive channels grouped within a list
        hichan_groups = [map(itemgetter(1),g)[::-1] for k,g in groupby(enumerate(self.hichans), lambda (i, x): i-x)][::-1] # reversed order (both of internal lists and globally)
        lochan_groups = [map(itemgetter(1),g) for k,g in groupby(enumerate(self.lochans), lambda (i, x): i-x)]
        print "High channels (grouped):",hichan_groups
        print "Low channels (grouped):",lochan_groups

        return hichan_groups,lochan_groups

    
    def write_cont_scripts(self):
        """
        Function to write RTS submission script in the case of a "standard" 
        contiguous bandwidth observation. This is relatively simple and just
        requires us to use the standard format RTS configuration file.
        """
        nnodes = 25 # number of required GPU nodes - 1 per coarse channels + 1 master node
        rts_batch = "RTS_{0}".format(self.cal_obsid)
        slurm_kwargs = {"partition":"gpuq", "workdir":"{0}".format(self.rts_out_dir), "time":"00:20:00", "nodes":"{0}".format(nnodes)}
        commands = list(self.body) # make a copy of body to then extend
        commands.append("aprun -n {0} -N 1  rts_gpu {1}".format(nnodes, self.rts_in_file))
        submit_slurm(rts_batch, commands, slurm_kwargs=slurm_kwargs, batch_dir=self.batch_dir, submit=self.submit)
 

    def get_subband_config(self, chan_groups, basepath, chan_type, count):
        """
        Function to make the approriate changes to a base copy of the RTS configuration scripts.
        This will interrogate the channels and decide on the "subbandIDs" and 
        "ObservationBaseFrequency" that are appropraite for each subband.
        """
        chan_file_dict = {} # used to keep track of which file needs how many nodes

        cc = 0
        for c in chan_groups:
            if len(c) == 1:
                # single channel group, write its own RTS config file
                subid = str(count+1)

                if chan_type == "low":
                    # offset is just how many channels along we are in the list
                    offset = cc
                elif chan_type == "high":
                    # offset is now how many subbandsIDs away from 24 we are
                    offset = 24-int(subid)
                else:
                    print "Invalid channel group type: must be \"low\" or \"high\". Aborting!"
                    sys.exit(1)
                
                # the base frequency is then
                basefreq = 1.28 *(c[0] - offset) - 0.625
    
                # and the coarse channel centre frequency is
                cfreq = 1.28*c[0]-0.625

                with open(self.rts_in_file,'rb') as f:
                    string = f.read()
                    # find the pattern and select the entire line for replacement
                    string = re.sub("ObservationFrequencyBase=.*\n","ObservationFrequencyBase={0}\n".format(basefreq), string)
                    # include the SubBandIDs tag 
                    string = re.sub("StartProcessingAt=0\n","StartProcessingAt=0\nSubBandIDs={0}\n\n".format(subid), string)

                fname = "{0}/rts_{1}_chan{2}.in".format(basepath, self.cal_obsid, c[0])
                chan_file_dict[fname] = 1 # this particular subband consists of only 1 channel
            
                with open(fname,'wb') as f:
                    f.write(string)
        
                print "Single channel :: (subband id, abs. chan, abs. freq) = ({0}, {1}, {2})".format(subid, c[0], cfreq)

                count += 1
                cc += 1

            elif len(c) > 1:
                # multiple consecutive channels
                subids = [str(count+i+1) for i in range(len(c))]
        
                # as with single channels, compute the offset
                if chan_type == "low":
                    offset = cc
                elif chan_type == "high":
                    # for high channels, there is now multiple offsets
                    offset = np.array([24-int(x) for x in subids])
                else:
                    print "Invalid channel group type: must be \"low\" or \"high\". Aborting!"
                    sys.exit(1)

                # basefreq is then the minmium of the calculated quantities
                basefreq = min(1.28 *(np.array(c) - offset) - 0.625)

                # multiple coarse channel centre frequencies
                cfreqs = 1.28*(np.array(c))-0.625

                with open(self.rts_in_file,'rb') as f:
                    string = f.read()
                    # find the pattern and select the entire line for replacement
                    string = re.sub("ObservationFrequencyBase=.*\n","ObservationFrequencyBase={0}\n".format(basefreq), string)
                    # include the SubBandIDs tag
                    string = re.sub("StartProcessingAt=0\n","StartProcessingAt=0\nSubBandIDs={0}\n\n".format(",".join(subids)), string)
        
                print "Multiple channels :: (subband ids, abs. chans, abs. freqs) = {0}".format(", ".join("({0}, {1}, {2})".format(i,j,k) for i,j,k in zip(subids,c,cfreqs)))

                if chan_type == "low":
                    fname = "{0}/rts_{1}_chan{2}-{3}.in".format(basepath, self.cal_obsid, min(c), max(c))
                elif chan_type == "high":
                    fname = "{0}/rts_{1}_chan{2}-{3}.in".format(basepath, self.cal_obsid, max(c), min(c))

                chan_file_dict[fname] = len(c)
                with open(fname,'wb') as f:
                    f.write(string)

                count += len(c)
                cc += len(c)

            else:
                print "Reached a channel group with no entries!? Aborting."
                sys.exit(1)


        return chan_file_dict,count

       
            
    def write_picket_fence_scripts(self):
        """
        Function to write RTS submission scripts in the case of a picket-fence
        observation. A significant amount of extra information needs to be 
        determined and placed in the RTS configuration files. There will be:
            
            1 RTS configuration file per set of adjacent single channels (subband)
        
        The only exception to that rule is where the 129 boundary is crossed, 
        in which case that subband will be split in two.
        """
        self.sort_obs_channels()
        hichan_groups, lochan_groups = self.construct_subbands()

        # write out the RTS config files and keep track of the number of nodes required for each
        count = 0
        lodict,count = self.get_subband_config(lochan_groups,self.rts_out_dir,"low",count)
        hidict,count = self.get_subband_config(hichan_groups,self.rts_out_dir,"high",count)
        chan_file_dict = lodict.copy()
        chan_file_dict.update(hidict)

        lolengths = set([len(l) for l in lochan_groups]) # figure out the unique lengths 
        hilengths = set([len(h) for h in hichan_groups])
        lengths = lolengths.union(hilengths) # combine the sets
                            
        # Now submit the RTS jobs
        print "Writing and submitting RTS jobs"
        
        for k,v in chan_file_dict.iteritems():
            nnodes = v + 1
            chans = k.split('_')[-1].split(".")[0]
            rts_batch = "RTS_{0}_{1}".format(self.cal_obsid, chans)
            slurm_kwargs = {"partition":"gpuq", "workdir":"{0}".format(self.rts_out_dir), "time":"00:45:00", "nodes":"{0}".format(nnodes)}
            commands = list(self.script_body) # make a copy of body to then extend
            commands.append("aprun -n {0} -N 1  rts_gpu {1}".format(nnodes, k))
            submit_slurm(rts_batch, commands, slurm_kwargs=slurm_kwargs, batch_dir=self.batch_dir, submit=self.submit)



    def run(self):
        """
        Only function that needs to be called after creating the RTScal object.
        Ensures functions are called in correct order to update/evalute attributes and 
        produce the RTS channel-specific files if required.
        """
        self.get_obs_channels() # fetch channels
        self.is_picket_fence()  # determine if picket-fence obs or not
        self.summary()

        if self.picket_fence:
            self.write_picket_fence_scripts()
        else:
            self.write_cont_scripts()

        print "Done!"


class BaseRTSconfig(object):
    """
    A class to hold the base information for a given RTS configuration file.
    """

    def __init__(self, obsid, cal_obsid, metafits, srclist, datadir, outdir, offline=False):
        self.obsid = obsid # target observation ID
        self.cal_obsid = cal_obsid # calibrator observation ID
        self.offline = offline # switch to decide if offline correlated data or not
        self.utctime = None # start UTC time
        self.nfine_chan = None # number of fine channels
        self.fine_cbw = None # fine channel bandwidth
        self.corr_dump_time = None # correaltor dump times (i.e. integration time)
        self.n_corr_dumps_to_average = None # number of integration times to use for calibration
        self.PB_HA = None # primary beam HA
        self.PB_DEC = None # primary beam Dec
        self.freq_base = None # frequency base for RTS
        self.JD = None # time base for RTS
        self.metafits_RTSform = None # modified metafits file name for RTS
        self.ArrayPositionLat = -26.70331940 # MWA latitude
        self.ArrayPositionLong = 116.6708152 # MWA longitude
        self.base_str = None # base string to be written to file, will be editted by RTScal

        # Check to make sure paths and files exist:
        # First, check that the actual data directory exists
        if os.path.isdir(datadir):
            self.data_dir = os.path.abspath(datadir)
        else:
            print "Given data directory doesn't exist. Aborting."
            sys.exit(1)
        
        # Then check if the specified output directory exists
        if outdir is None:
            # this is the default
            self.output_dir = "/group/mwaops/vcs/{0}/cal/{1}/rts".format(self.obsid, self.cal_obsid)
            mdir(self.output_dir, "RTS")
        else:
            # mdir handles if the directory already exists
            self.output_dir = os.path.abspath(outdir + "/rts")
            print "non standard RTS output path", self.output_dir
            mdir(self.output_dir, "RTS")
        
        # Then check that the metafits file exists
        if os.path.isfile(metafits) is False:
            # file doesn't exist
            print "Given metafits file does not exist. Aborting."
            sys.exit(1)
        elif "_ppds" not in metafits:
            # file doesn't have the correct naming convention
            print "Looks like you have an old-style metafits. You'll need to download the new version, which is named like: {0}_metafits_ppds.fits. Aborting.".format(obsid)
            sys.exit(1)
        else:
            self.metafits = os.path.abspath(metafits)

        # the check that the source list exists
        if os.path.isfile(srclist) is False:
            # file doesn't exist
            print "Given source list file does not exist. Aborting."
            sys.exit(1)
        else:
            self.source_list = os.path.abspath(srclist)



        # set some RTS flags based on if we have offline correlated data or not
        if self.offline:
            self.useCorrInput = 1
            self.readDirect = 0
        else:
            self.useCorrInput = 0
            self.readDirect = 1


    def get_info_from_data_header(self):
        # first determine the UTC time from the file name
        file_glob = "{0}/*_gpubox*.fits".format(self.data_dir)
        files = sorted(glob.glob(file_glob))
        len_files = len(files)
        if len_files == 0:
            print "No *_gpubox*.fits files found in {0}. Aborting.".format(self.data_dir)
            sys.exit(1)
        elif len_files % 24 != 0:
            print "Number of *_gpubox*.fits files is not divisible by 24!?. Aborting."
            sys.exit(1)

        first_file = files[0]
        self.utctime = os.path.splitext(os.path.basename(first_file))[0].split("_")[1]


        if self.offline is False:
            # now figure out how much data we have in total by counting the number of data HDUs
            # then open the header to access the frequency spacing and inetgration times
            hdulist = fits.open(first_file)
            nfine_chan = int(hdulist[1].header['NAXIS2'])
            inttime = float(hdulist[1].header['INTTIME'])

            # coarse channel BW is 1.28 MHz
            self.fine_cbw = 1.28/nfine_chan 
            self.nfine_chan = nfine_chan

<<<<<<< HEAD
            # the correlator dump time is whatever is in the header
            self.corr_dump_time = inttime

            ngroups = len_files/24
            fgrouped = np.array(np.array_split(files, ngroups))
            ndumps = 0
            for item in fgrouped[:,0]:
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

            self.fine_cbw = 1.28/nfine_chan
            self.nfine_chan = nfine_chan
            self.corr_dump_time = inttime
            
            # for offline correlation, each file is one integration time - there has been no concatenation
            # TODO: this assumes that the offline correlator ALWAYS produces 1 second FITS files
            ndumps = len(fits.open(first_file)) * len_files/24
            self.n_dumps_to_average = ndumps


    def construct_base_string(self):
        # get observation information from database
        obsinfo = getmeta(service='obs', params={'obs_id':str(self.obsid)})
        
        # get the RA/Dec pointing for the primary beam
        ra_pointing_degs = obsinfo['metadata']['ra_pointing']
        dec_pointing_degs = obsinfo['metadata']['dec_pointing']

        # convert times using our timeconvert and get LST and JD 
        # TODO: need to make this not have to call an external script
        timeconvert = distutils.spawn.find_executable("timeconvert.py")
        cmd = "{0} --datetime={1}".format(timeconvert, self.utctime)
        print cmd
        time_cmd = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

        for line in time_cmd.stdout:
            if "LST" in line:
                lststring,lstflag = line.split()
                hh,mm,ss = lststring.split(":")

            if "JD" in line:
                jdflag,jd = line.split()

        lst_in_hours = float(hh) + float(mm) / 60.0 + float(ss) / 60.0**2

        # set the HA of the image centre to the primary beam HA
        self.JD = jd
        PB_HA_HOURS = (ra_pointing_degs / 360.0) * 24.0
        self.PB_HA = lst_in_hours - PB_HA_HOURS
        self.PB_DEC = dec_pointing_degs

        print " * Primary beam: HA = {0} hrs, Dec = {1} deg".format(self.PB_HA, self.PB_DEC)
        print " * JD = {0}".format(self.JD)
        print " * LST = {0}:{1}:{2} = {3} hrs".format(hh, mm, ss, lst_in_hours)

        # get the lowest frequency channel
        freqs = obsinfo['rfstreams']['0']['frequencies']
        start_channel = freqs[0]

        self.freq_base = start_channel * 1.28e6 - 0.64e6 + 15e3 # frequency base in Hz
        self.freq_base /= 1e6 # convert to MHz

        print " * Frequency lower edge = {0} MHz".format(self.freq_base)

        # make metafits file formatted for RTS
        self.metafits_RTSform = self.metafits.split("_metafits_ppds.fits")[0]


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
        FieldOfViewDegrees=1""".format(os.path.realpath(self.data_dir), \
                                        self.readDirect, \
                                        self.useCorrInput, \
                                        self.metafits_RTSform, \
                                        self.freq_base, \
                                        self.JD, \
                                        self.PB_HA, \
                                        self.PB_DEC, \
                                        self.fine_cbw, \
                                        self.nfine_chan, \
                                        self.n_dumps_to_average, \
                                        self.corr_dump_time, \
                                        self.ArrayPositionLat, \
                                        self.ArrayPositionLong, \
                                        self.source_list)



        return file_str


    def run(self):
        self.get_info_from_data_header()
        self.base_str = self.construct_base_string()



=======
##################
## END OF CLASS ##
##################
>>>>>>> ed65b95009a550e892ca34f417763e42866e3a06
if __name__ == '__main__':
    ####################
    ## OPTION PARSING ##
    ####################
    parser = argparse.ArgumentParser(description="Calibration script for creating RTS configuration files in the VCS pulsar pipeline")

<<<<<<< HEAD
    parser.add_argument("-o", "--obsid", type=int, help="Observation ID of target.", required=True)
    parser.add_argument("-O", "--cal_obsid", type=int, help="Calibrator observation ID for the target observation. \
                                                                        Can be the same as -o/--obs if using in-beam visibilities.", required=True)
    parser.add_argument("-m", "--metafits", type=str, help="Path to the relevant calibration observation metafits file.", required=True)
    parser.add_argument("-s", "--srclist", type=str, help="Path to the desired source list.", required=True)
    
    #parser.add_argument("--rts_in_file", type=str, help="Path to the base RTS configuration file - created by write_rts_in_file.py", required=True)
    parser.add_argument("--gpubox_dir", type=str, help="Where the *_gpubox*.fits files are located")

    parser.add_argument("--rts_output_dir", type=str, help="Parent directory where you want the /rts directory and /batch directory to be created. \
                                                    ONLY USE IF YOU WANT THE NON-STANDARD LOCATIONS.", default=None)
    parser.add_argument("--offline", action='store_true', help="Tell the RTS to read calibrator data in the offline correlated data format.")
    parser.add_argument("--submit", action='store_true', help="Switch to allow SLURM job submission")
    args = parser.parse_args()

    baseRTSconfig = BaseRTSconfig(args.obsid, args.cal_obsid, args.metafits, args.srclist, args.gpubox_dir, args.rts_output_dir, args.offline)
    baseRTSconfig.run()

    #calobj = RTScal(args.obs, args.cal_obsid, args.rts_in_file, args.rts_output_dir, args.submit)
    #calobj.run()
=======
    parser.add_argument("-o", "--obs", type=int, metavar="OBSID", help="Observation ID of target.", required=True)
    parser.add_argument("-O", "--cal_obs", type=int, metavar="CAL_OBSID", help="Calibrator observation ID for the target observation. \
                                                                        Can be the same as -o/--obs if using in-beam vsibilities.", required=True)
    parser.add_argument("--rts_in_file", type=str, help="Path to the base RTS configuration file - created by write_rts_in_file.py", required=True)
    parser.add_argument("--rts_output_dir", type=str, help="Parent directory where you want the /rts directory and /batch directory to be created. \
                                                    ONLY USE IF YOU WANT THE NON-STANDARD LOCATIONS.", default=None)
    parser.add_argument("--submit", action='store_true', help="Switch to allow SLURM job submission")
    args = parser.parse_args()

    calobj = RTScal(args.obs, args.cal_obs, args.rts_in_file, args.rts_output_dir, submit)
    calobj.run()
>>>>>>> ed65b95009a550e892ca34f417763e42866e3a06
