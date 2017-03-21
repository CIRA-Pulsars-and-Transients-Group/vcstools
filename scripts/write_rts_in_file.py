#!/usr/bin/env python

"""
Script to build the required input configuration file for RTS jobs 
"""

import argparse
import glob
import os
import sys
import subprocess

import urllib
import urllib2
import json
import distutils.spawn

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


def write_rts_in_file(obsid,utc_time,data_dir,metafits_file,srclist_file,rts_fname):
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
	fid.write("BaseFilename=%s/*_gpubox\n" % data_dir)
	fid.write("ReadGpuboxDirect=1\n")
	fid.write("UseCorrelatorInput=0\n")
	fid.write("\n") 

	# read the metafits file
	fid.write("ReadMetafits=1\n")
	fid.write("MetafitsFilename=%s\n" % metafits_file)
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
	fid.write("ObservationFrequencyBase=%f\n" % freq_base )
	fid.write("ObservationTimeBase=%s\n" % jd )
	fid.write("ObservationPointCentreHA=%f\n" % PB_HA )
	fid.write("ObservationPointCentreDec=%f\n" % PB_DEC )
	fid.write("ChannelBandwidth=0.01\n")
	fid.write("NumberOfChannels=128\n")
	fid.write("\n")
	
	# set the time integration and iterations
	fid.write("CorrDumpsPerCadence=16\n")
	fid.write("CorrDumpTime=2.0\n")
	fid.write("NumberOfIntegrationBins=3\n")
	fid.write("NumberOfIterations=1\n")
	fid.write("\n")
	fid.write("StartProcessingAt=0\n")
	fid.write("\n")

	# array information
	fid.write("ArrayPositionLat=%f\n" % ArrayPositionLat )
	fid.write("ArrayPositionLong=%f\n" % ArrayPositionLong )
	fid.write("\n")
	fid.write("ArrayFile=%s/array_file.txt\n" % data_dir)
	fid.write("ArrayNumberOfStations=128\n")

	# source list and calibration info
	fid.write("SourceCatalogueFile=%s\n" % srclist_file)
	fid.write("NumberOfCalibrators=1\n")
	fid.write("NumberOfSourcesToPeel=0\n")
	fid.write("calBaselineMin=20.0\n")
	fid.write("calShortBaselineTaper=40.0\n" )
	fid.write("FieldOfViewDegrees=1\n")
	fid.write("\n")
	fid.close()

	print "\nRTS configuration setup written to: {0}".format(rts_fname)




	
parser = argparse.ArgumentParser(description="Gather calibration information and prepare a RTS configuration file.")
parser.add_argument("-o",metavar="obsID",type=str,help="The observation ID")
parser.add_argument("-f",metavar="metafits_file",type=str,help="Full path to the metafits file for the observation")
parser.add_argument("-s",metavar="srclist_file",type=str,help="Full path to the source list file created by srclist_by_beam.py")
parser.add_argument("--gpubox_dir",type=str,help="Full path to where the *_gpubox files are located")
parser.add_argument("--output_dir",type=str,help="Full path to where you want the RTS configuration file to be written")


if len(sys.argv)==1:
	# no arguments passed, print help and quit
	parser.print_help()
	sys.exit(0)

args = parser.parse_args()

# figure out RTS config file name
fname = "{0}/rts_{1}.in".format(os.path.abspath(args.output_dir),args.o)

# figure out the start time from the first gpubox file
print "Finding start time from gpubox files"
gpubox_file_glob = "{0}/*_gpubox*.fits".format(os.path.abspath(args.gpubox_dir))
gpubox_files = sorted(glob.glob(gpubox_file_glob))
first_gpubox_file = gpubox_files[0]
utctime = os.path.splitext(os.path.basename(first_gpubox_file))[0].split("_")[1]
print "Determined start time is {0} from {1}".format(utctime,first_gpubox_file)

gpuboxes = os.path.abspath(args.gpubox_dir)
metafits = os.path.abspath(args.f)
srclist = os.path.abspath(args.s)

write_rts_in_file(args.o,utctime,gpuboxes,metafits,srclist,fname)
