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
from process_vcs import submit_slurm, getmeta, mdir
from itertools import groupby
from operator import itemgetter


#################################
## CALIBRATION CONTAINER CLASS ##
#################################
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



##################
## END OF CLASS ##
##################

####################
## OPTION PARSING ##
####################
parser = argparse.ArgumentParser(description="Calibration script for creating RTS configuration files in the VCS pulsar pipeline")

parser.add_argument("-o", "--obs", type=int, metavar="OBSID", help="Observation ID of target.", required=True)
parser.add_argument("-O", "--cal_obs", type=int, metavar="CAL_OBSID", help="Calibrator observation ID for the target observation. \
                                                                    Can be the same as -o/--obs if using in-beam vsibilities.", required=True)
parser.add_argument("--rts_in_file", type=str, help="Path to the base RTS configuration file - created by write_rts_in_file.py", required=True)
parser.add_argument("--rts_output_dir", type=str, help="Parent directory where you want the /rts directory and /batch directory to be created. \
                                                ONLY USE IF YOU WANT THE NON-STANDARD LOCATIONS.", default=None)
parser.add_argument("--no_submit", action='store_true', help="Switch to make the program only write the RTS config scripts and not submit them")
args = parser.parse_args()

# check whether to submit jobs to queue
if args.no_submit:
    submit = False
else:
    submit = True

calobj = RTScal(args.obs, args.cal_obs, args.rts_in_file, args.rts_output_dir, submit)
calobj.run()
