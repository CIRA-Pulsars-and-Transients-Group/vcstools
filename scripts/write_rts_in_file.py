#!/usr/bin/env python

"""
Script to build the required input configuration file for RTS jobs 
"""

import argparse
import glob
import os
import sys
import subprocess
import math
from astropy.io import fits
import numpy as np

import urllib
import urllib2
import json
import distutils.spawn

from process_vcs import mdir

BASEURL = 'http://mwa-metadata01.pawsey.org.au/metadata/'

def getmeta(service='obs', params=None):
	"""
	Function to call a JSON web service and return a dictionary:
	Given a JSON web service ('obs', find, or 'con') and a set of parameters as
	a Python dictionary, return a Python dictionary containing the result.
	Taken verbatim from http://mwa-lfd.haystack.mit.edu/twiki/bin/view/Main/MetaDataWeb
	"""

	if params:
		data = urllib.urlencode(params)  # Turn the dictionary into a string with encoded 'name=value' pairs
	else:
		data = ''

	if service.strip().lower() in ['obs', 'find', 'con']:
		service = service.strip().lower()
	else:
		print "invalid service name: %s" % service
		return

	try:
		result = json.load(urllib2.urlopen(BASEURL + service + '?' + data))
	except urllib2.HTTPError as error:
		print "HTTP error from server: code=%d, response:\n %s" % (error.code, error.read())
		return
	except urllib2.URLError as error:
		print "URL or network error: %s" % error.reason
		return

	return result


def write_rts_in_file(obsid,utc_time,data_dir,metafits_file,srclist_file,rts_fname,chan_bw,nchan,dumptime,ndumps):
	"""
	Gather the respective information and write the input conifiguration file for the RTS
	"""
	# get observation information from database
	obsinfo = getmeta(service='obs', params={'obs_id':str(obsid)})
	
	# get the RA/Dec pointing for the primary beam
	ra_pointing_degs = obsinfo['metadata']['ra_pointing']
    	dec_pointing_degs = obsinfo['metadata']['dec_pointing']

	# convert times using our timeconvert and get LST and JD 
	timeconvert = distutils.spawn.find_executable("timeconvert.py")
	cmd = "{0} --datetime={1}".format(timeconvert,utctime)
	print cmd
	time_cmd = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
	
	for line in time_cmd.stdout:
		if "LST" in line:
			lststring,lstflag = line.split()
			hh,mm,ss = lststring.split(":")

		if "JD" in line:
			jdflag,jd = line.split()
	
	lst_in_hours = float(hh) + float(mm) / 60.0 + float(ss) / 60.0**2
	
	# set the HA of the image centre to the primary beam HA
	PB_HA_HOURS = (ra_pointing_degs / 360.0) * 24.0
	PB_HA = lst_in_hours - PB_HA_HOURS
	PB_DEC = dec_pointing_degs

	print " * Primary beam: HA = {0} hrs, Dec = {1} deg".format(PB_HA,PB_DEC)
	print " * JD = {0}".format(jd)
	print " * LST = {0}:{1}:{2} = {3} hrs".format(hh,mm,ss,lst_in_hours)

	# get the lowest frequncy channel
	freqs = obsinfo['rfstreams']['0']['frequencies']
	start_channel = freqs[0]

	freq_base = start_channel * 1.28e6 - 0.64e6 + 15e3 	# frequency base in Hz
	freq_base /= 1e6 					# convert to MHz

	print " * Frequency lower edge = {0} MHz".format(freq_base)

	# make metafits file formatted for RTS
	metafits_rtsform = metafits_file.split("_metafits_ppds.fits")[0]	
	
	# define array position
	ArrayPositionLat = -26.70331940
    	ArrayPositionLong = 116.6708152	

	# write the file
	try:
        	fid = open(rts_fname, 'w')
    	except IOError:
        	print "Can't open file {0} for writing.".format(rts_fname)
        	sys.exit(0)
	
	# read the gpubox files
	fid.write("\n")
	fid.write("ReadAllFromSingleFile=\n")
	fid.write("BaseFilename={0}/*_gpubox\n".format(os.path.realpath(data_dir))) # there are symlinks here, so expand those
	# reading gpubox files from offline correlator (as of 23 March 2017) actually requires
	fid.write("ReadGpuboxDirect=0\n")
        fid.write("UseCorrelatorInput=1\n") # this can handle BOTH online and offline correaltor products (though it throughs errors which are seemingly harmless)
	#fid.write("ReadGpuboxDirect=1\n") # read gpubox files directly from disk
	#fid.write("UseCorrelatorInput=0\n") # don't expect an input stream from the correlator
	fid.write("\n") 

	# read the metafits file
	fid.write("ReadMetafitsFile=1\n")
	fid.write("MetafitsFilename={0}\n".format(metafits_rtsform))
	fid.write("\n")

	# set the calibration configuration
	fid.write("DoCalibration=\n")
	fid.write("doMWArxCorrections=1\n")
	fid.write("doRawDataCorrections=1\n") 
	fid.write("doRFIflagging=0\n")
	fid.write("useFastPrimaryBeamModels=1\n")
        fid.write("generateDIjones=1\n")
        fid.write("applyDIcalibration=1\n")
	fid.write("UsePacketInput=0\n")
	fid.write("UseThreadedVI=1\n")
	fid.write("\n")
	
	# observation details
	fid.write("ObservationFrequencyBase={0}\n".format(freq_base))
	fid.write("ObservationTimeBase={0}\n".format(jd))
	fid.write("ObservationPointCentreHA={0}\n".format(PB_HA))
	fid.write("ObservationPointCentreDec={0}\n".format(PB_DEC))
	fid.write("ChannelBandwidth={0}\n".format(chan_bw))
	fid.write("NumberOfChannels={0}\n".format(nchan))
	fid.write("\n")
	
	# set the time integration and iterations
	fid.write("CorrDumpsPerCadence={0}\n".format(ndumps))
	fid.write("CorrDumpTime={0}\n".format(dumptime))
	fid.write("NumberOfIntegrationBins=3\n")
	fid.write("NumberOfIterations=1\n")
	fid.write("\n")
	fid.write("StartProcessingAt=0\n")
	fid.write("\n")

	# array information
	fid.write("ArrayPositionLat={0}\n".format(ArrayPositionLat))
	fid.write("ArrayPositionLong={0}\n".format(ArrayPositionLong))
	fid.write("ArrayNumberOfStations=128\n")
	fid.write("\n")

	# turns out the RTS still needs this option, but it can be left blank as the relevent info is retrieved from the metafits file
	fid.write("ArrayFile=\n")
	fid.write("\n")

	# source list and calibration info
	fid.write("SourceCatalogueFile={0}\n".format(srclist_file))
	fid.write("NumberOfCalibrators=1\n")
	fid.write("NumberOfSourcesToPeel=0\n")
	fid.write("calBaselineMin=20.0\n")
	fid.write("calShortBaselineTaper=40.0\n")
	fid.write("FieldOfViewDegrees=1\n")
	fid.write("\n")
	
	fid.close()

	print "\nRTS configuration setup written to: {0}".format(rts_fname)
	

def write_flag_files(odir, metafits_file, nchan):
	"""
	Given the output directory, write initial flagging files based on bad tiles in metafits and number of fine channels
	""" 
	metafits = fits.open(metafits_file) # read metafits file
	bad_tiles = metafits[1].data['Flag'][::2] # both polarisation recorded, so we just want every second value 
	bad_tiles = np.where(bad_tiles == 1)[0]
	flagged_tiles = "{0}/flagged_tiles.txt".format(odir) 
	with open(flagged_tiles,'w') as fid:
		for b in bad_tiles:
			fid.write("{0}\n".format(b))

	# figure out how many edge channels to flag based on the fact that with 128, we flag the edge 8
	ntoflag = 8 * nchan/128 
	chans = np.arange(nchan)
	start_chans = chans[:ntoflag]
	end_chans = chans[-ntoflag:]
	center_chan = [nchan/2 - 1] # zero based, so -1 from nchan/2
	bad_chans = np.hstack((start_chans,center_chan,end_chans))
	flagged_channels = "{0}/flagged_channels.txt".format(odir)
	with open(flagged_channels,'w') as fid:
		for b in bad_chans:
			fid.write("{0}\n".format(b))
	
	
	
parser = argparse.ArgumentParser(description="Gather calibration information and prepare a RTS configuration file.")
parser.add_argument("-o",metavar="obsID",type=str,help="The observation ID",required=True)
parser.add_argument("-f",metavar="metafits_file",type=str,help="Where the observation metafits file is located",required=True)
parser.add_argument("-s",metavar="srclist_file",type=str,help="Where the source list file created by srclist_by_beam.py is located",required=True)
parser.add_argument("--fine_chan_bw",type=float,help="Fine channel bandwidth for the observation in MHz (default: 0.01)",default=0.01)
parser.add_argument("--corr_dump_time",type=float,help="Correlator dump time in seconds (default: 2.0)",default=2.0)
parser.add_argument("--ndumps_to_average",type=int,help="Number of correlator dumps to average together (default: 16)",default=16)
parser.add_argument("--gpubox_dir",type=str,help="Where the *_gpubox files are located (default: `pwd`)",default='`pwd`')
parser.add_argument("--output_dir",type=str,help="Where you want the RTS configuration file to be written (default: `pwd`)",default='`pwd`')

if len(sys.argv)==1:
	# no arguments passed, print help and quit
	parser.print_help()
	sys.exit(0)

args = parser.parse_args()
print "\n"
# check for default gpubox_dir and make absolute paths
if args.gpubox_dir == "`pwd`":
	gpubox_dir = os.path.abspath(os.getcwd())
else:
	gpubox_dir = os.path.abspath(args.gpubox_dir)

# check for default output_dir and make absolute paths
if args.output_dir == "`pwd`":
	output_dir = os.path.abspath(os.getcwd())
else:
	output_dir = os.path.abspath(args.output_dir) 


# check to make sure output_dir and gpubox_dir exists
if os.path.isdir(output_dir) == False:
	print "output directory does not exist at:"
	print "\t {0}".format(output_dir)
	print "Aborting here."
	sys.exit(0)
elif os.path.isdir(gpubox_dir) == False:
	print "gpubox directory does not exist at:"
        print "\t {0}".format(gpubox_dir)
        print "Aborting here."
        sys.exit(0)

# make absolute path to metafits, check it exists and make sure it's the "new" version
metafits = os.path.abspath(args.f)
if os.path.isfile(metafits) == False:
	print "metafits file does not exist at:"
        print "\t {0}".format(metafits)
        print "Aborting here."
        sys.exit(0)
elif "_ppds" not in metafits:
	print "Looks like you have an old-style metafits. You'll need to download the new version, which is named like: {0}_metafits_ppds.fit".format(args.o)
	print "Aborting here."
	sys.exit(0)


# check that the frequency resolution makes sense
if args.fine_chan_bw < 0.01:
	print "!!!WARNING!!! :: Typically for VCS observations, the finest channel bandwith we have is 10 kHz."
	print "Continuing, but this frequency resolution will probably cause errors later down the pipeline.\n"

# calculate how many fine channels per coarse channel there are
nfine_per_coarse = int(math.ceil(1.28/args.fine_chan_bw))
print "Determined there are {0} fine channels per coarse channel\n".format(nfine_per_coarse)

# calcaulte how much data in seconds are required
req_data_len = args.corr_dump_time * args.ndumps_to_average
print "Options passed require at least {0} seconds of data".format(req_data_len)
print "This amount needs to be less than the total amount of data supplied for calibration, else the solutions will be degraded and/or the RTS will crash.\n"
# TODO: need to figure out a way that we can smartly asses how much data is available (from metafits?). This is fine for online, but not so clear cut for offline products.
if req_data_len < 10:
	print "!!!WARNING!!! :: This is a small amount of data to get a decent calibration solution from."
	print "Continuing, but the RTS may produce poor quality solutions and/or crash.\n"


# figure out the start time from the first gpubox file
print "Finding start time from gpubox files"
gpubox_file_glob = "{0}/*_gpubox*.fits".format(gpubox_dir)
print "Globbing {0}".format(gpubox_file_glob)
gpubox_files = sorted(glob.glob(gpubox_file_glob))

if len(gpubox_files) == 0:
	print "No gpubox files found in {0}".format(gpubox_dir)
	print "Aborting here."
	sys.exit(0)

first_gpubox_file = gpubox_files[0]
# get the utc time assuming gpubox files are named: [prefix]_UTCtime_gpubox*.fits
utctime = os.path.splitext(os.path.basename(first_gpubox_file))[0].split("_")[1]
print "Determined start time is {0} from {1}\n".format(utctime,first_gpubox_file)


#first, create an rts/ subdirectory in the output_dir
write_dir = "{0}/rts".format(output_dir)
if os.path.isdir(write_dir):
	print "rts subdirectory already prepared"
else:
	mdir(write_dir,"rts")

# set RTS configuration file name and location
fname = "{0}/rts_{1}.in".format(write_dir,args.o)

# check the source list file 
srclist = os.path.abspath(args.s)
if os.path.isfile(srclist) == False:
	# if it's not at the path given, check in the rts directory
	print "srclist file does not exist at:"
	print "\t {0}".format(srclist)
	print "Checking in rts/ subdirectory"
	if os.path.basename(srclist) in os.listdir(write_dir):
		print "Found in rts/ subdirectory"
	else:
		print "Could not find srclist, please create a new one"
		sys.exit(0)
else:	
	# move the srclist file into the rts directory
	os.system("mv {0} {1}".format(srclist,write_dir))
	
srclist = "{0}/rts/{1}".format(os.path.dirname(srclist),os.path.basename(srclist))
print "Found source list at {0}".format(srclist)

# now write the configuration file
write_rts_in_file(args.o,utctime,gpubox_dir,metafits,srclist,fname,args.fine_chan_bw,nfine_per_coarse,args.corr_dump_time,args.ndumps_to_average)
# now write the initial flagged_tiles.txt and flagged_channels.txt files
write_flag_files(write_dir,metafits,nfine_per_coarse)
