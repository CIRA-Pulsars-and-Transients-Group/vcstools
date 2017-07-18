#!/usr/bin/env python
import os
import sys
from reorder_chans import *
from process_vcs import submit_slurm, getmeta, mdir



class Calibration(object):

    def __init__(self, obsid=None, cal_obsid=None, rts_in_file=None, rts_out_dir=None):
        self.obsid       = obsid
        self.cal_obsid   = cal_obsid
        self.rts_in_file = rts_in_file
        # Set up the product directory according to options passed
        if rts_out_dir is not None:
            print "!!WARNING!! You are not using the standard directory structure. Be sure you know wht you're doing!"
            self.rts_out_dir = rts_out_dir
            self.product_dir = self.rts_out_dir
            mdir(self.product_dir, "RTS output")
        else:
            # use the default: /group/mwaops/vcs/{obsid}/cal/{cal_obsid}/rts
            self.product_dir = "/group/mwaops/vcs/{0}/cal/{1}/rts".format(self.obsid, self.cal_obsid)
            self.rts_out_dir = self.product_dir
            mdir(self.product_dir, "RTS output")

        self.channels = None


    def print_summary(self):
        print "Summary of Calibration object contents:"
        print "\tObservation ID:           ",self.obsid
        print "\tCalibrator Observation ID:",self.cal_obsid
        print "\tRTS configuration file:   ",self.rts_in_file
        print "\tRTS output directory:     ",self.rts_out_dir
        print "\tProduct directory:        ",self.product_dir

    def _set_obsid(self, obsid):
        self.obsid = obsid

    def _set_cal_obsid(self, cal_obsid):
        self.cal_obsid = cal_obsid

    def _set_rts_in_file(self, rts_in_file):
        self.rts_in_file = rts_in_file

    def _set_rts_out_dir(self, rts_out_dir):
        self.rts_out_dir = rts_out_dir

    def _set_product_dir(self, product_dir):
        self.product_dir = product_dir


    def construct_metafile(self):
        metafile = "{0}/{1}.meta".format(self.product_dir, self.cal_obsid)
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

    







a = Calibration(1099415632,1099415632,"test_rts.in","/astro/mwaops/bmeyers/test2",)
a.print_summary()
a.construct_metafile()
print a.channels

