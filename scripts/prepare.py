#!/usr/bin/env python

import sys
import glob
import os
import subprocess as subp
import string
import re
import getopt
import array
import time
import calendar
import astropy.io.fits as pyfits

import urllib
import urllib2
import json

import datetime

# Append the service name to this base URL, eg 'con', 'obs', etc.
BASEURL = 'http://ngas01.ivec.org/metadata/'

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

u2d = "u2d"
get_delays = "get_delays"
read_pfb = "read_pfb"
mwac_offline = "mwac_offline"
outdir = "notset"
timeconvert = "timeconvert.py"
get_info = "sql_get_observation_info.py"

def sfreq(freqs):

    if len(freqs) != 24:
        print "There are not 24 coarse chans defined for this obs. Got: %s" % freqs
        return

    freqs.sort()   # It should already be sorted, but just in case...
    lowchans = [f for f in freqs if f <= 128]
    highchans = [f for f in freqs if f > 128]
    highchans.reverse()
    freqs = lowchans + highchans
    return freqs

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

#def get_frequencies(obs_id):
#    obsinfo = getmeta(service='obs', params={'obs_id':str(obs_id)})
#    freq_array = obsinfo['rfstreams']['0']['frequencies']
#    return sfreq(freq_array)

def get_frequencies(metafits):
    hdulist = pyfits.open(metafits)
    freq_array = hdulist[0].header['CHANNELS']
    return sfreq(freq_array.split(','))

def build_rts_in_file(obs_id,utctime,data_dir,rts_filename):

    # we need to pick up the pointing HA
    # the time
    # the start frequency
    # we will generate one of these for the lowest frequency channel

    # Get the PB pointing from the sql_get_observation_info.py script

    obsinfo = getmeta(service='obs', params={'obs_id':str(obs_id)})

    ra_pointing_degs = obsinfo['metadata']['ra_pointing']
    dec_pointing_degs = obsinfo['metadata']['dec_pointing']

    # Need the LST

    # Best to use the timeconvert functionality as we are already using that

    cmd_line = "%s --datetime=%s" % (timeconvert,utctime)
    print cmd_line
    time_cmd = subp.Popen(cmd_line,shell=True,stdout=subp.PIPE)
    for line in time_cmd.stdout:
        if "LST" in line:
            (lststring,lstflag) = line.split()
            (hh,mm,ss) = lststring.split(':')
        if "JD" in line:
            (JDflag,jd) = line.split()

    lst_in_hours = float(hh) + float(mm)/60.0 + float(ss)/(60.0*60.0)

    # if the RA of the image centre is not set, set it to the centre of the PB FOV
    PB_RA_HOURS = (ra_pointing_degs/360.0) * 24.0
    PB_HA = lst_in_hours - PB_RA_HOURS
    PB_Dec = dec_pointing_degs



    print " * Primary beam: HA = %f hrs, Dec = %f deg" % (PB_HA, PB_Dec)
    print " * JD = %s" % jd
    print " * LST = %s:%s:%s = %f (hrs) " % (hh,mm,ss, lst_in_hours)


    # lowest frequency is from the start channel
    freq_array = obsinfo['rfstreams']['0']['frequencies']
    StartChannel = freq_array[0]

    freqbase = StartChannel*1.28e6 - 0.64e6 + 15e3
    freqbase = freqbase/1e6

    print " * Freq lower edge = %f MHz" % freqbase

    ArrayPositionLat = -26.70331940
    ArrayPositionLong = 116.67081524

    # open the tle file
    try:
        fid = open(rts_filename, 'w')
    except IOError:
        print 'Can\'t open file \'' + rts_filename + '\' for writing.'
        sys.exit(0)

    fid.write( "\n" )
    fid.write( "ReadAllFromSingleFile=\n")
    fid.write( "BaseFilename=%s/*_gpubox\n" % data_dir)
    fid.write( "\n" )

    fid.write( "ObservationFrequencyBase=%f\n" % freqbase )
    fid.write( "ChannelBandwidth=0.04\n")
    fid.write( "NumberOfChannels=32\n")

    fid.write( "ObservationTimeBase=%s\n" % jd )
    fid.write( "\n" )
    fid.write( "ObservationPointCentreHA=%f\n" % PB_HA )
    fid.write( "ObservationPointCentreDec=%f\n" % PB_Dec )
    fid.write( "CorrDumpsPerCadence=8\n")

    fid.write( "CorrDumpTime=1\n")
    fid.write( "NumberOfIntegrationBins=3\n")
    fid.write( "NumberOfIterations=25\n")
    fid.write( "\n" )
    fid.write( "StartProcessingAt=0\n")
    fid.write( "\n" )
    fid.write( "// -------------------------------------------------------------------------- //\n" )
    fid.write( "\n" )
    fid.write( "DoCalibration=\n")
    fid.write( "doMWArxCorrections=1\n" )
    fid.write( "doRawDataCorrections=1\n")
    fid.write( "useFastPrimaryBeamModels=1\n")
    fid.write( "generateDIjones=1\n")
    fid.write( "applyDIcalibration=1\n")
    fid.write( "\n" )
    fid.write( "UseCorrelatorInput=1\n" )
    fid.write( "UsePacketInput=0\n" )
    fid.write( "UseThreadedVI=1\n" )
    fid.write( "\n" )
    fid.write( "ArrayPositionLat=%f\n" % ArrayPositionLat )
    fid.write( "ArrayPositionLong=%f\n" % ArrayPositionLong )
    fid.write( "\n" )
    fid.write( "ArrayFile=array_file.txt\n")
    fid.write(  "ArrayNumberOfStations=128\n")
    fid.write( "SourceCatalogueFile=\n")

    fid.write( "NumberOfCalibrators=1\n")
    fid.write( "NumberOfSourcesToPeel=0\n")
    fid.write( "calBaselineMin=20.0\n")
    fid.write( "calShortBaselineTaper=40.0\n" )
    fid.write( "FieldOfViewDegrees=1\n")
    fid.write( "\n" )
    fid.close()

def usage (opts={}):

    print "prepare a for the beamformer/correlator run\n"
    print "Options"
    print "-f Metafile\n"
    print "-r <J2000 RA>\n"
    print "-d <J2000 Dec>\n"
    print "-s [outdir] Submits the correlator job - does not get delays\n"
    print "-g [outdir] get delays (needs the directory of the calibration solutions)\n"
    print "-r [0|1] build rts .in file\n"
    print "-m [0|1] mode 0==old 1==new\n"
    print "-o RTS mode only\n" 
    
    sys.exit()

if __name__ == '__main__':

    from sys import argv


    
    old_mode = 0
    new_mode = 1
    edge_num = 20;

    
    freq_array=[]
    f=[]

    the_options = {'ra': "05:34:34", 'dec': "22:00:10", 'delays': False, 'submit': False, 'rts' : 0, 'rts_only' : 0, 'extn' : "pfb", 'mode' : 1,'metafile': "null"}

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hr:d:e:f:g:m:or:s:")
    except getopt.GetoptError:
        usage(the_options)
        sys.exit()
    finally:
        if len(sys.argv) < 2:
            usage(the_options)


    for opt,arg in opts:
        
        if (opt == "-h"):
            usage(the_options)
        elif (opt == "-r"):
            the_options['ra'] = arg
        elif (opt == "-d"):
            the_options['dec'] = arg
        elif (opt == "-e"):
            the_options['extn'] = arg
        elif (opt == "-f"):
            the_options['metafile'] = arg
        elif (opt == "-g"):
            the_options['delays'] = True
            outdir = arg
        elif (opt == "-o"):
            the_options['rts_only'] = True
        elif (opt == "-r"):
            the_options['rts'] = int(arg)
        elif (opt == "-s"):
            the_options['submit'] = True
            outdir = arg
            corr_batch_file_root = outdir + "/correlator_run";
        elif (opt == "-m"):
            the_options['mode'] = int(arg)

    print "Pointing %s %s\n" % (the_options['ra'],the_options['dec'])

    if (the_options['mode'] == 1):
        old_mode = 0
        new_mode = 1
    elif( the_options['mode'] == 0):
        old_mode = 1
        new_mode = 0
    
    if (the_options['submit'] == True):
        the_options['delays'] = False
        print "Submitting Job (no delays)"
    if (the_options['rts'] == True):
        the_options['delays'] = False;
        print "Build RTS in mode (no delays)"


    if (the_options['rts_only'] == True):

        corr_files_glob = "*.fits"
        c_f = sorted(glob.glob(corr_files_glob))
        first_corr_file = c_f[0]
        (current_time,ext) = os.path.splitext(os.path.basename(first_corr_file))
        (obsid,utctime,gpubox,num) = string.split(current_time,'_')
        rts_file = "./%s_%s_rts.in" % (obsid,utctime)
        build_rts_in_file(obsid,utctime,".",rts_file)
        sys.exit(0)


    files_glob = "*.%s" % the_options['extn']
    f = sorted(glob.glob(files_glob))
    file = f[0]
    obsid=0

    print file
    (current_time,ext) = os.path.splitext(os.path.basename(file))

    
    if (old_mode == 1):
        (obsid,c,chan,x_time) = string.split(current_time,'_')


        if (int(x_time) < 1000000000):
            # filename is broken add the offset
            unix_time = int(x_time) + 1000000000;
        else:
            unix_time = int(x_time)
 
            print unix_time

            cmd_line = "%s %s" % (u2d,unix_time)

            print cmd_line
            u2d_cmd = subp.Popen(cmd_line,shell=True,stdout=subp.PIPE)
            for line in u2d_cmd.stdout:
                if "GPS" in line:
                    (unix_time,gpstime,zulutime) = line.split()
        
            utctime = zulutime[4:-1]

            freq_channel = int(chan[-2])*10 + int(chan[-1])
            freq_channel_index = freq_channel - 1

    if (new_mode == 1):

        import astropy
        from astropy.time import Time
        (obsid,gpstime,chan) = string.split(current_time,'_')
        print gpstime
        t = Time(int(gpstime.rstrip()), format='gps', scale='utc')
        utctime =  t.datetime.strftime('%Y-%m-%dT%H:%M:%S')

        freq_channel = int(chan[2])*100 + int(chan[3])*10 + int(chan[4])

    freq_array = get_frequencies("../%s.metafits" % obsid)
    print freq_array
    print freq_channel
    gpubox_label = 0;
    for i,chan in enumerate(freq_array):
        if (freq_channel == int(chan)):
            gpubox_label = (i+1)


    if (freq_channel < 25):

        freq_Hz = int(freq_array[freq_channel_index]) * 1.28e6 - 0.64e6
    else:
        freq_Hz = freq_channel * 1.28e6 - 0.64e6
        for relative,absolute in enumerate(freq_array):
            if (freq_channel == absolute):
                freq_channel = relative + 1
                break

    if (the_options['submit'] == True):
    
        corr_batch = "%s_%s_ch%d" % (corr_batch_file_root,obsid,freq_channel)
        with open(corr_batch, 'w') as batch_file:
            batch_file.write("#!/bin/bash -l\n#SBATCH --nodes=1\n#SBATCH --export=NONE\n")
            batch_file.write("module load cudatoolkit\nmodule load cfitsio\n")


        mkdir_line = "mkdir %s/%s" % (outdir,obsid)
        subp.Popen(mkdir_line,shell=True,stdout=subp.PIPE)
    
        to_corr = 0;
        for file in f:
            corr_line = ""
            corr_file = "%s/%s/%s" % (outdir,obsid,obsid)
		# run the correlator
            (current_time,ext) = os.path.splitext(os.path.basename(file))
            if (old_mode == 1):
                (obsid,c,chan,x_time) = string.split(current_time,'_')
        
                if (int(x_time) < 1000000000):
                    # filename is broken add the offset
                    unix_time = int(x_time) + 1000000000;
                else:
                    unix_time = int(x_time)

                corr_line = " aprun -n 1 -N 1 %s -o %s -s %d -r 1 -i 100 -e %d -f 128 -n 4 -c %02d -d %s/%s\n" % (mwac_offline,corr_file,unix_time,edge_num,freq_channel,os.getcwd(),file)

            else:

                (obsid,gpstime,chan) = string.split(current_time,'_')

                import astropy
                from astropy.time import Time
                (obsid,gpstime,chan) = string.split(current_time,'_')
                t = Time(int(gpstime), format='gps', scale='utc')
                time_str =  t.datetime.strftime('%Y-%m-%d %H:%M:%S')

                current_time = time.strptime(time_str, "%Y-%m-%d  %H:%M:%S")
                unix_time = calendar.timegm(current_time)

                corr_line = " aprun -n 1 -N 1 %s -o %s -s %d -r 1 -i 100 -f 128 -n 4 -c %02d -d %s/%s\n" % (mwac_offline,corr_file,unix_time,gpubox_label,os.getcwd(),file)
                    
            #print corr_line
            with open(corr_batch, 'a') as batch_file:
                batch_file.write(corr_line)
            
            to_corr = to_corr+1

            outfile = "%s.pfb" % (file)

            if (os.path.isfile(outfile)):
                print "already created %s" % outfile
            else:
                if (old_mode == 1):
                    pfb_line = "%s -i %s -a 128 -n 88  -o %s -4 " % (read_pfb,file,outfile)
                    subp.call(pfb_line,shell=True)

        secs_to_run = datetime.timedelta(seconds=10*to_corr)
        batch_submit_line = "sbatch --time=%s --partition=gpuq %s" % (str(secs_to_run),corr_batch)
        print batch_submit_line
        submit_cmd = subp.Popen(batch_submit_line,shell=True,stdout=subp.PIPE)
        jobid=""
        for line in submit_cmd.stdout:
            if "Submitted" in line:
                (word1,word2,word3,jobid) = line.split()


#now we have to wait until this job is finished before we move on
        queue_line = "squeue -j %s\n" % jobid

        print queue_line

        finished = False

        while not finished:
#assume finished
            queue_cmd = subp.Popen(queue_line,shell=True,stdout=subp.PIPE)
            finished = True
            for line in queue_cmd.stdout:
                print line
                if jobid in line:
# batch job still in the queue
                    finished = False;
                    print "Not finished"
                time.sleep(1)

# build the rts file for these obs -
    if (the_options['rts'] == True):
        corr_dir = "%s/%s" % (outdir,obsid)
        corr_files_glob = "%s/*.fits" % (corr_dir)
        c_f = sorted(glob.glob(corr_files_glob))
        first_corr_file = c_f[0]
        (current_time,ext) = os.path.splitext(os.path.basename(first_corr_file))
        (obsid,utctime,gpubox,num) = string.split(current_time,'_')
        rts_file = "%s/%s_%s_rts.in" % (corr_dir,obsid,utctime)
        if (os.path.isfile(rts_file)):
            print "already created %s" % rts_file
        else:
            build_rts_in_file(obsid,utctime,corr_dir,rts_file)

    if (the_options['delays'] == True):
        
        DI_file = "%s/%s/DI_JonesMatrices_node0%02d.dat" % (outdir,obsid,gpubox_label)
        print DI_file

        if (old_mode == 1):
            if (os.path.isfile(DI_file)):
                delays_line = "%s -a ./ -b %d -j %s -m %s -i -p -z %s -o %s -f %s -e %d -n 88 -w 10000 -r %s -d %s" % (get_delays,len(f),DI_file,the_options['metafile'],utctime,obsid,freq_Hz,edge_num,the_options['ra'],the_options['dec'])
            else:
                delays_line = "%s -a ./ -b %d -i -p -z %s -o %s -f %s -e %d -n 88 -w 10000 -r %s -d %s -m %s" % (get_delays,len(f),utctime,obsid,freq_Hz,edge_num,the_options['ra'],the_options['dec'],the_options['metafile'])
        elif (new_mode == 1):
            if (os.path.isfile(DI_file)):
                delays_line = "%s -a ./ -b %d -j %s -m %s -i -p -z %s -o %s -f %s -n 128 -w 10000 -r %s -d %s" % (get_delays,len(f),DI_file,the_options['metafile'],utctime,obsid,freq_Hz,the_options['ra'],the_options['dec'])
            else:
                delays_line = "%s -a ./ -b %d -i -p -z %s -o %s -f %s -n 128 -w 10000 -r %s -d %s -m %s" % (get_delays,len(f),utctime,obsid,freq_Hz,the_options['ra'],the_options['dec'],the_options['metafile'])

        print delays_line
        
        subp.call(delays_line,shell=True)
        rts_flags_file = "%s/%s/flagged_tiles.txt" % (outdir,obsid)
        if (os.path.isfile(rts_flags_file)) :
            rts_flags = []
            with open (rts_flags_file,"r") as rts:
                for line in rts:
                    try:
                        rts_flags.append(int(line))
                    except:
                        print "none integer in line"
            beamer_flags = []
            with open ("./flags.txt","r") as beamer:
                for line in beamer:
                    beamer_flags.append(float(line))
            
            for tile in rts_flags:
                beamer_flags[2*tile] = 0.0
                beamer_flags[2*tile+1] = 0.0

        with open("./flags.txt","w") as new_flags:
            for entry in beamer_flags:
                s = str(entry) + "\n"
                new_flags.write(s)


