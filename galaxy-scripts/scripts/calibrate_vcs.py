#!/usr/bin/env python
import os
import sys
import numpy as np
import re
from reorder_chans import *
from process_vcs import submit_slurm, getmeta, mdir

#################################
## CALIBRATION CONTAINER CLASS ##
#################################
class Calibration(object):

    def __init__(self, obsid=None, cal_obsid=None, rts_in_file=None, rts_out_dir=None):
        self.obsid        = obsid
        self.cal_obsid    = cal_obsid
        self.rts_in_file  = rts_in_file
        self.channels     = None
        self.picket_fence = None
        
        # Set up the product directory according to options passed
        if rts_out_dir is not None:
            print "!!WARNING!! You are not using the standard directory structure. Be sure you know wht you're doing!"
            self.rts_out_dir = rts_out_dir
            self.product_dir = self.rts_out_dir
        else:
            # use the default: /group/mwaops/vcs/{obsid}/cal/{cal_obsid}/rts
            self.product_dir = "/group/mwaops/vcs/{0}/cal/{1}".format(self.obsid, self.cal_obsid)
            self.rts_out_dir = self.product_dir

        # re-define what the RTS output is based on the product directory
        self.rts_out_dir  = "{0}/rts".format(self.product_dir)
        mdir(self.rts_out_dir, "RTS output")
        self.batch_dir    = "{0}/batch".format(self.product_dir)
        mdir(self.batch_dir, "Batch")
        
    def summary(self):
        print "Summary of Calibration object contents:"
        print "\tObservation ID:           ",self.obsid
        print "\tCalibrator Observation ID:",self.cal_obsid
        print "\tRTS configuration file:   ",self.rts_in_file
        print "\tRTS output directory:     ",self.rts_out_dir
        print "\tBatch directory:          ",self.batch_dir
        print "\tObservation channels:     ",self.channels
        print "\tIs picket fence?          ",self.picket_fence


    def is_picket_fence(self):
        ch_offset = self.channels[-1] - self.channels[0]
        if ch_offset == len(self.channels)-1:
            print "The channels are consecutive: this is a normal observation"
            self.picket_fence = False
        else:
            print "The channels are NOT consecutive: this is a picket-fence observation"
            self.picket_fence = True


    def construct_metafile(self):
        metafile = "{0}/{1}.meta".format(self.product_dir, self.cal_obsid)
        metafile_exists = False
        print metafile

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

##################
## END OF CLASS ##
##################



def compute_basefreq(offset, channels):
    if len(channels) == 1:
        return 1.28 *(channels[0] - offset) - 0.625
    elif len(channels) > 1:
        return min(1.28 * (np.array(channels) - offset) - 0.625)
    else:
        print "Received channel list with 0 length. Aborting."
        sys.exit(1)


def write_rts_in_files(chan_groups, basepath, chan_type, calobj, count=0):
    chan_file_dict = {} # used to keep track of which file needs how many nodes

    cc = 0
    for c in chan_groups:
        if len(c) == 1:
            #single channel group, write its own rts_in file
            subid = str(count+1)

            if chan_type == "low":
                offset = cc
            elif chan_type == "high":
                offset = 24-int(subid)
            else:
                print "Invalid channel group type: must be \"low\" or \"high\". Aborting!"
                sys.exit(1)

            basefreq = compute_basefreq(offset, c)

            # use re.compile to make search expressions
            with open(calobj.rts_in_file,'rb') as f:
                string = f.read()
                # find the pattern and select the entire line for replacement
                string = re.sub("ObservationFrequencyBase=.*\n","ObservationFrequencyBase={0}\n".format(basefreq), string)
                # include the SubBandIDs tag 
                string = re.sub("StartProcessingAt=0\n","StartProcessingAt=0\nSubBandIDs={0}\n\n".format(subid), string)

            fname = "{0}/rts_{1}_chan{2}.in".format(basepath, calobj.cal_obsid, c[0])
            chan_file_dict[fname] = 1 # this particular rts_in file has only 1 channel
        
            with open(fname,'wb') as f:
                f.write(string)
    
            print "Single channel:: (subband id, abs. chan, abs. freq) = ({0}, {1}, {2}) {3}".format(subid, c[0], basefreq, cc)

            count += 1
            cc += 1

        elif len(c) > 1:
            # multiple consecutive channels
            subids = [str(count+i+1) for i in range(len(c))]
    
            if chan_type == "low":
                offset = cc
            elif chan_type == "high":
                offset = np.array([24-int(x) for x in subids])
            else:
                print "Invalid channel group type: must be \"low\" or \"high\". Aborting!"
                sys.exit(1)

            freqs = 1.28*(np.array(c)-offset)-0.625
            basefreq = compute_basefreq(offset, c)

            # use re.compile to make search expressions
            with open(calobj.rts_in_file,'rb') as f:
                string = f.read()
                # find the pattern and select the entire line for replacement
                string = re.sub("ObservationFrequencyBase=.*\n","ObservationFrequencyBase={0}\n".format(basefreq), string)
                # include the SubBandIDs tag
                string = re.sub("StartProcessingAt=0\n","StartProcessingAt=0\nSubBandIDs={0}\n\n".format(",".join(subids)), string)
    
            print "Multiple channels:: (subband ids, abs. chans, abs. freqs) = {0}".format(", ".join("({0}, {1}, {2})".format(i,j,k) for i,j,k in zip(subids,c,freqs)))

            if chan_type == "low":
                fname = "{0}/rts_{1}_chan{2}-{3}.in".format(basepath, calobj.cal_obsid, min(c), max(c))
            elif chan_type == "high":
                fname = "{0}/rts_{1}_chan{2}-{3}.in".format(basepath, calobj.cal_obsid, max(c), min(c))

            chan_file_dict[fname] = len(c)
            with open(fname,'wb') as f:
                f.write(string)

            count += len(c)
            cc += len(c)

        else:
            print "Reached a channel group with no entries!? Aborting."
            sys.exit(1)


    return chan_file_dict,count





def construct_RTS_scripts(calobj):
    
    #define some of the header/body text to go into the slurm submission script
    body = []
    body.append("module load cudatoolkit")
    body.append("module load cfitsio")
    body.append("cd {0}".format(calobj.product_dir))
   

    if calobj.picket_fence == False:
        nnodes = 25
        rts_batch = "RTS_{0}".format(calobj.cal_obsid)
        slurm_kwargs = {"partition":"gpuq", "workdir":"{0}".format(calobj.rts_out_dir), "time":"00:20:00", "nodes":"{0}".format(nnodes)}
        commands = list(body) # make a copy of body to then extend
        commands.append("aprun -n {0} -N 1  rts_gpu {1}".format(nnodes,calobj.rts_in_file))
        submit_slurm(rts_batch, commands, slurm_kwargs=slurm_kwargs, batch_dir=calobj.batch_dir,submit=False)
    elif calobj.picket_fence == True:
        # separate channels into "high bands" and "low bands"
        hichans = [c for c in calobj.channels if c>128]
        lochans = [c for c in calobj.channels if c<=128]
        print "High channels:",hichans
        print "Low channels:",lochans

        # get the consecutive chunks within each channel subselection
        print "Grouping individual channels into consecutive chunks (if possible)"
        from itertools import groupby
        from operator import itemgetter
        # create a list of lists with consecutive channels grouped wihtin a list
        hichan_groups = [map(itemgetter(1),g)[::-1] for k,g in groupby(enumerate(hichans), lambda (i, x): i-x)][::-1] # reversed order (both of internal lists and globally)
        lochan_groups = [map(itemgetter(1),g) for k,g in groupby(enumerate(lochans), lambda (i, x): i-x)]
        print "High channels (grouped):",hichan_groups
        print "Low channels (grouped):",lochan_groups

        basepath = os.path.dirname(calobj.rts_in_file)

        # write out the RTS in files and keep track of the number of nodes required for each
        count = 0
        lodict,count = write_rts_in_files(lochan_groups,basepath,"low",calobj,count)
        hidict,count = write_rts_in_files(hichan_groups,basepath,"high",calobj,count)
        chan_file_dict = lodict.copy()
        chan_file_dict.update(hidict)

        lolengths = set([len(l) for l in lochan_groups]) # figure out the unique lengths 
        hilengths = set([len(h) for h in hichan_groups])
        lengths = lolengths.union(hilengths) # combine the sets
                            
        # Now submit the RTS jobs
        print "Writing and submitting RTS jobs"
        
        for k,v in chan_file_dict.iteritems():
            nnodes = v + 1
            channels = k.split('_')[-1].split(".")[0]
            rts_batch = "RTS_{0}_{1}".format(calobj.cal_obsid, channels)
            slurm_kwargs = {"partition":"gpuq", "workdir":"{0}".format(calobj.rts_out_dir), "time":"00:45:00", "nodes":"{0}".format(nnodes)}
            commands = list(body) # make a copy of body to then extend
            commands.append("aprun -n {0} -N 1  rts_gpu {1}".format(nnodes,k))
            submit_slurm(rts_batch, commands, slurm_kwargs=slurm_kwargs, batch_dir=calobj.batch_dir,submit=False)

    else:
        print "Calibration object picket fence flag was never set! Aborting."
        sys.exit(1)





a = Calibration(1099415632,1099415632,"/astro/mwaops/bmeyers/test2/rts/rts-chan15-20.in","/astro/mwaops/bmeyers/test2")
#a = Calibration(1099414416,1099414416,"test_rts.in","/astro/mwaops/bmeyers/test2")
a.summary()
a.construct_metafile()
a.is_picket_fence()
a.summary()
construct_RTS_scripts(a)
