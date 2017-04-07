#! /usr/bin/env python 

"""
Author: Nicholas Swainston
Creation Date: /05/2016

The MWA Pulsar Database was created by David Pallot and he wrote the database call functions.

This code is used to submit observations to the MWA Pulsar Database. The goal is to calculate all need values (flux density, width, scattering) for each observation and  submit them to the pulsar database without having to manually input them.

To load mwa_pulsar_client 
export /group/mwaops/PULSAR/src/mwa-voltage/
will be need to be included in your .bashrc 
"""
#TODO actually get mwa_pulsar client to load properly

__author__ = 'Nicholas Swainston'
__date__ = '2016-05-12'


import os
import requests
import argparse
import matplotlib.pyplot as plt #compatibiality issues or isn't loaded for python 2.7
import numpy as np
import urllib
import urllib2
import json
import math
import subprocess
from mwapy.pb import primary_beam
import ephem
from mwapy import ephem_utils,metadata
from scipy import constants
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
#from mwa_pulsar_client import client
#removes urllib3 ssl platform warnings
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()
from mwapy.pb import get_Tsys #TODO change when moved to MWA tools
from mwa_pulsar_client import client
import sys
from shutil import copyfile as cp


def psrcat(addr, auth, pulsar):
    """
    Return the pulsar details from psrcat.
    Args:
        pulsar: name of the pulsar.
    Returns:
        psrcat details as a python dict.
    Exception:
        pulsar not found or bad input.
    """
    path = 'https://{0}/{1}/'.format(addr, 'psrcat')
    payload = {'name': pulsar, 'format': 'json'}
    r = requests.post(url = path, auth = auth, data = payload)
    r.raise_for_status()
    return r.json()

    
def getmeta(service='obs', params=None):
    """
    getmeta(service='obs', params=None)
    
    Given a JSON web service ('obs', find, or 'con') and a set of parameters as
    a Python dictionary, return the RA and Dec in degrees from the Python dictionary.
    """
    BASEURL = 'http://mwa-metadata01.pawsey.org.au/metadata/'
    if params:
        data = urllib.urlencode(params)  # Turn the dictionary into a string with encoded 'name=value' pairs
    else:
        data = ''
    #Validate the service name
    if service.strip().lower() in ['obs', 'find', 'con']:
        service = service.strip().lower()
    else:
        print "invalid service name: %s" % service
        return
    #Get the data
    try:
        result = json.load(urllib2.urlopen(BASEURL + service + '?' + data))
    except urllib2.HTTPError as error:
        print "HTTP error from server: code=%d, response:\n %s" % (error.code, error.read())
        return
    except urllib2.URLError as error:
        print "URL or network error: %s" % error.reason
        return
    #Return the result dictionary
    return result
    
def sex2deg( ra, dec):
    """
    sex2deg( ra, dec)
    ra - the right ascension in HH:MM:SS
    dec - the declination in DD:MM:SS
    
    Convert sexagesimal coordinates to degrees.
    """ 
    c = SkyCoord( ra, dec, frame='icrs', unit=(u.hourangle,u.deg))
    
    # return RA and DEC in degrees in degrees
    return [c.ra.deg, c.dec.deg]
    
def get_from_bestprof(file_loc):
    with open(file_loc,"rb") as bestprof:
        lines = bestprof.readlines()
        obsid = lines[0][22:32]
        try:
            obsid = int(obsid)
        except:
            obsid = lines[0][28:38]
        pulsar = str(lines[1][26:-1])
        if not pulsar.startswith('J'):
            pulsar = 'J' + pulsar
        dm = lines[14][22:-1]
        period = lines[15][22:-1]
        period, period_uncer = period.split('  +/- ')
        obslength = float(lines[6][22:-1])*float(lines[5][22:-1])
        #get profile list
        orig_profile = []
        for l in lines[27:]:
            discard, temp = l.split()
            orig_profile.append(float(temp))
        bin_num = len(orig_profile)
        profile = np.zeros(bin_num)
        min_prof = min(orig_profile)
        for p in range(len(orig_profile)):
            profile[p] = orig_profile[p] - min_prof
        #maybe centre it around the pulse later        
    return [obsid, pulsar, dm, period, period_uncer,obslength, profile, bin_num]


def get_beam_power(obsid_data,
                   start,
                   stop,
                   sources,
                   dt=100,
                   centeronly=True,
                   verbose=False, option='a'):
    """
    obsid_data = [obsid,ra, dec, time, delays,centrefreq, channels]
    sources=[names,coord1,coord2] #astropy table coloumn names

    Calulates the power (gain at coordinate/gain at zenith) for each source and if it is above
    the min_power then it outputs it to the text file.

    """
    print "Calculating beam power"
    obsid,ra, dec, time, delays,centrefreq, channels = obsid_data
    
    starttimes=np.arange(start,stop,dt)
    stoptimes=starttimes+dt
    stoptimes[stoptimes>time]=time
    Ntimes=len(starttimes)
    midtimes=float(obsid)+0.5*(starttimes+stoptimes)

    mwa = ephem_utils.Obs[ephem_utils.obscode['MWA']]
    # determine the LST
    observer = ephem.Observer()
    # make sure no refraction is included
    observer.pressure = 0
    observer.long = mwa.long / ephem_utils.DEG_IN_RADIAN
    observer.lat = mwa.lat / ephem_utils.DEG_IN_RADIAN
    observer.elevation = mwa.elev

    if not centeronly:
        PowersX=np.zeros((len(sources),
                             Ntimes,
                             len(channels)))
        PowersY=np.zeros((len(sources),
                             Ntimes,
                             len(channels)))
        # in Hz
        frequencies=np.array(channels)*1.28e6
    else:
        PowersX=np.zeros((len(sources),
                             Ntimes,1))
        PowersY=np.zeros((len(sources),
                             Ntimes,1))
        frequencies=[centrefreq]

    RAs, Decs = sex2deg(sources[0][0],sources[0][1])
    theta=[]
    phi=[]
    for itime in xrange(Ntimes):
        obstime = Time(midtimes[itime],format='gps',scale='utc')
        observer.date = obstime.datetime.strftime('%Y/%m/%d %H:%M:%S')
        LST_hours = observer.sidereal_time() * ephem_utils.HRS_IN_RADIAN

        HAs = -RAs + LST_hours * 15
        Azs, Alts = ephem_utils.eq2horz(HAs, Decs, mwa.lat)
        # go from altitude to zenith angle
        # go from altitude to zenith angle
        if option == 'a':
            beam_string = "analytic"
            theta=np.radians(90-Alts)
            phi=np.radians(Azs)
            
            for ifreq in xrange(len(frequencies)):
                rX,rY=primary_beam.MWA_Tile_analytic(theta, phi,
                                                     freq=frequencies[ifreq], delays=delays,
                                                     zenithnorm=True,
                                                     power=True)
                PowersX[:,itime,ifreq]=rX
                PowersY[:,itime,ifreq]=rY
                
        if option == 'd':
            beam_string = "advanced"
            theta=np.array([np.radians(90-Alts)])
            phi=np.array([np.radians(Azs)])
            
            for ifreq in xrange(len(frequencies)):
                rX,rY=primary_beam.MWA_Tile_advanced(theta, phi,
                                                     freq=frequencies[ifreq], delays=delays,
                                                     zenithnorm=True,
                                                     power=True)
                PowersX[:,itime,ifreq]=rX
                PowersY[:,itime,ifreq]=rY
        
        if option == 'e':
            beam_string = "full_EE"
            #theta=np.radians(90-Alts)
            #phi=np.radians(Azs)
            theta.append(np.radians(90-Alts))
            phi.append(np.radians(Azs))
            #h5filepath='/group/mwaops/PULSAR/src/MWA_Tools/mwapy/pb/mwa_full_embedded_element_pattern.h5'
            
    if option == 'e':
        theta = np.array([theta])
        phi = np.array([phi])
        rX,rY=primary_beam.MWA_Tile_full_EE(theta, phi,
                                             freq=frequencies[0], delays=delays,
                                             zenithnorm=True,
                                             power=True)
        PowersX=rX
        PowersY=rY

    #Power [#sources, #times, #frequencies]
    Powers=0.5*(PowersX+PowersY)
    print 'The powers changing with time:'
    print Powers
    avg_power = np.mean(Powers)
    return avg_power

def get_pulsar_ra_dec(pulsar):
    #Gets the ra and dec from the output of PSRCAT
    cmd = ['psrcat', '-c', 'Raj', pulsar]
    output = subprocess.Popen(cmd,stdout=subprocess.PIPE).communicate()[0]
    temp = []
    lines = output.split('\n')
    print lines
    for l in lines[4:-1]: 
        columns = l.split()
        if len(columns) > 1:
            ra = columns[1]
    cmd = ['psrcat', '-c', 'Decj', pulsar]
    output = subprocess.Popen(cmd,stdout=subprocess.PIPE).communicate()[0]
    temp = []
    lines = output.split('\n')
    for l in lines[4:-1]: 
        columns = l.split()
        if len(columns) > 1:
            dec = columns[1]
    return [ra, dec]


class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()  
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)



web_address = 'mwa-pawsey-volt01.pawsey.ivec.org'
auth = ('mwapulsar','veovys9OUTY=')

parser = argparse.ArgumentParser(description="""
R|This code is used to submit observations to the MWA Pulsar Database. The goal is to calculate all need values (flux density, width, scattering) for each observation and  submit them to the pulsar database without having to manually input them.
Two common uses are to simply upload a calibration file (if there is no observation on the database it will create a blank one that won't be able to be updated without deleting it at this time)\n
> submit_to_database.py -c /scratch2/mwaops/vcs/1133775752/calibrate/1133775752-cal.tar.gz -o 1133775758 -p J2325-0530\n
Another is to give the .bestprof file to calculate all values of the detection\n
> submit_to_database.py -b 1120799848_PSR_0837-4135.pfd.bestprof\n
To process the needed databse files can be done by:
> submit_to_database.py  -o 1133775758 -p J2325-0530 --u_diagnostics --u_timeseries --u_archive
""", formatter_class=SmartFormatter)
parser.add_argument('-o','--obsid',type=str,help='The obsid.')
parser.add_argument('--incoh',action='store_true',help='Used for incoherent detections to accurately calculate gain. Default is coherent.')
parser.add_argument('--andre',action='store_true',help="Used for calibrations done using Andre Offrina's tools. Default is RTS.")
parser.add_argument('-p','--pulsar',type=str,help='The pulsar J name.')
parser.add_argument('-f','--fits_files',type=str,help='The fits files location to be used to for creating extra files. Recommended to end in *.fits and surrounded by quotation marks.')#TODO word this better
parser.add_argument('--cal_id',type=str,help='The obsid of the calibrator.')
calcargs = parser.add_argument_group('Detection Calculation Options', 'All the values of a pulsar detection for the MWA pulsar databse can be calculated by this code using either a .besprof file or a (Soon to be implimented dspsr equivalent). If neither are included then all values will be left as null and the observation(/detection) can still be used to upload files without the calculation being performed')
calcargs.add_argument('-b','--bp_file',type=str,help='The location of the .bestprof file. Using this option will Can be used instead of inputing the obsid and pulsar.')
calcargs.add_argument('--start',type=str,help='The start time of the detection in seconds. For example: 0')
calcargs.add_argument('--stop',type=str,help='The stop time of the detection in seconds. For example: 1000')
#parser.add_argument('-n','--no_calc',action='store_true',help="Exits the code before the calculations begin. Used to only upload files after a previous calculation has been completed.")
uploadargs = parser.add_argument_group('Upload options', 'The different options for each file type that can be uploaded to the pulsar database. Will cause an error if the wrong file type is being uploaded.')
uploadargs.add_argument('-a','--archive',type=str,help="The archive file location to be uploaded to the database. Expects a single file that is the output of dspsr using the pulsar's ephemeris.")
uploadargs.add_argument('-t','--timeseries',type=str,help='The timeseries file location to be uploaded to the database. Expects a single file that is the output of dspsr in single pulse mode (the -s option).')
uploadargs.add_argument('-d','--diagnostics',type=str,nargs='*',help="The diagnostics file location to be uploaded to the database. Expects three arguments of the locations of the presto's prepfold output, an intergrated pulse profile using dspsr's pav output and a waterfall plot of pulse phase vs frequency using dspsr's psrplot. If don't have all there right None instead, for example: --diagnostics None <waterfall_plot>.ps None")
uploadargs.add_argument('-c','--calibration',type=str,help='The calibration solution file location to be uploaded to the database. Expects a single file so please zip or tar up multiple files.')
dspsrargs = parser.add_argument_group('Dspsr Calculation Options', "These options will process the needed files for the database but will have to be uploaded using this script once checked. The files will be uploaded automatically when the dspsr scripts are tested more") #TODO remove when I'm confident with dspsr
dspsrargs.add_argument('--u_archive',action='store_true',help='The archive to be processed with dspsr and then uploaded to the database')
dspsrargs.add_argument('--u_timeseries',action='store_true',help='The timeseries to be processed with dspsr and then uploaded to the database')
dspsrargs.add_argument('--u_diagnostics',action='store_true',help='The diagnostics to be processed with dspsr and then uploaded to the database')
args=parser.parse_args()



#print client.detection_list(web_address, auth)

#defults for incoh and calibrator type
if args.incoh:
    incoh = True
    calibrator_type = None
else:
    incoh = False
    if not args.cal_id:
        print "Please include --cal_id for coherent observations"
        quit()
    if args.andre:
        calibrator_type = 1
    else:
        calibrator_type = 2

#get info from .bestprof file
if args.bp_file:
    obsid, pulsar, dm, period, period_uncer, time_detection, profile, num_bins = get_from_bestprof(args.bp_file)
    print pulsar
elif args.obsid and args.pulsar:
    num_bins = 128
    obsid = args.obsid
    pulsar = args.pulsar
elif args.obsid and args.calibration and not (args.archive or args.timeseries or args.diagnostics or args.u_archive or args.u_timeseries or args.u_diagnostics):
    obsid = args.obsid
else:
    print "Please us either --obsid and --pulsar or --bp_file"
    sys.exit(0)
    
"""
detects = client.detection_list(web_address,  auth)
for d in detects
    if d[u'observationid'] == obsid:
        print d
"""


if args.fits_files:
    fits_files_loc = args.fits_files
else:
    fits_files_loc = '/scratch2/mwaops/vcs/'+str(obsid)+'/fits/*.fits'

if args.pulsar or args.bp_file:
    #Checks to see if the pulsar is already on the database
    pul_list_dict = client.pulsar_list(web_address, auth)
    pul_list_str = ''
    for p in pul_list_dict:
        pul_list_str = pul_list_str + p[u'name']
    if pulsar in pul_list_str:
        print 'This pulsar is already on the database'
        #gets Ra and DEC from database (gets last one which is incorrect)
        #pul_ra = p[u'ra']
        #pul_dec = p[u'dec']
        
        #gets Ra and DEC from PSRCAT
        pul_ra , pul_dec = get_pulsar_ra_dec(pulsar)
    else:
        print 'Congratulations you have detected ' + pulsar + ' for the first time with the MWA'
        #gets Ra and DEC from PSRCAT
        pul_ra , pul_dec = get_pulsar_ra_dec(pulsar)
        #then adds it to the database
        client.pulsar_create(web_address, auth, name = pulsar, ra = pul_ra, dec = pul_dec)    


    

#get meta data from obsid
obsid,ra_obs,dec_obs,time_obs,delays,centrefreq,channels = get_Tsys.get_obs_metadata(obsid)
minfreq = float(min(channels))
maxfreq = float(max(channels))

bandwidth = 30720000.


if args.bp_file:
    #check if the obs time is entire obs. The meta data will round down to the nearest 200 seconds
    #(someitmmes 100 depending on the obs type) 
    entire_obs = False
    if 0. < (time_detection - time_obs) < 199.:
        entire_obs = True
        enter = 0.
        exit = time_detection

    if not entire_obs:
        #check if there is an input start and stop time
        if time_detection > (time_obs *5.):
            #some of the comissioning data has terrible metadata 
            #going to assume it's for the entire data set
            enter = 0.
            exit = float(time_detection)
        elif (not args.start) and args.stop:
            print "No start time input so assumed it's 0"
            enter = 0.
            exit = args.stop
            exit = float(exit)
        elif args.start and (not args.stop):
            print "No stop time input so assumed it's the end of the observation"
            exit = float(time_detection)
            enter = args.start
            enter = float(enter)
        elif args.start and args.stop:
            enter = args.start
            exit = args.stop
            enter = float(enter)
            exit = float(exit)
        else: #if there isn't assume that it's an autofold obs
            #check for a complete analytic beam output
            if os.path.isfile('/scratch2/mwaops/pulsar/incoh_census/'+obsid+'/'+obsid+'_analytic_beam.txt'):
                beam_list=open('/scratch2/mwaops/pulsar/incoh_census/'+obsid+'/'+obsid+\
                              '_analytic_beam.txt').readlines()
                for line in beam_list:
                    if line.startswith(str(pulsar)):
                        psrline=line.split()
                        #fraction of obs when the pulsar enters the beam
                        enter=float(psrline[2])*float(time_obs) 
                        #times by obs duration to turn into seconds
                        exit=float(psrline[3]) * float(time_obs)   
            else: #create file if there isn't one
                os.system('find_pulsar_in_obs.py -p '+pulsar+' Jfake -o '+obsid+' --all_volt --output ./')
                beam_list=open(obsid+'_analytic_beam.txt').readlines()
                for line in beam_list:
                    if line.startswith(str(pulsar)):
                        psrline=line.split()
                        #fraction of obs when the pulsar enters the beam
                        enter=float(psrline[2]) * float(time_obs) 
                        exit=float(psrline[3]) * float(time_obs) 
                        if 0. < (time_detection - exit) < 199.: #same bad metadata correction as above
                            exit = time_detection
                #os.remove(obsid+'_analytic_beam.txt')
        input_detection_time = exit - enter
        if not int(input_detection_time) == int(time_detection):
            print "Input detection time does not equal the dectetion time of the .bestprof file"
            print "Input time: " + str(input_detection_time)
            print "Bestprof time: " + str(time_detection)
            option = raw_input("Would you like to continue the program with possibliy incorrect values" +\
                               " (yes) or exit the program (no). (y/n)")
            if ('y' in option) or ('Y' in option) or (option == ''):
                print "Program will continue"
            elif ('n' in option) or ('N' in option):
                quit()
            else:
                print "Unknown input so program is exiting"
                quit()





    #Old gain calc
    """
    #once you know the time of the obs calc the average power
    avg_power = get_beam_power([obsid,ra_obs, dec_obs, time_detection, delays ,centrefreq, channels],
                                enter, exit, [[pul_ra, pul_dec]], option = 'e')
                                
    
    #calc system temp by assiming attena temperature is a linear fucntion of frequency
    t_sys = skytemp - 0.5 * centrefreq / 1000000.  + 125. #TODO check with engineers if this t_ant is accurate
    print "Old T_sys:",t_sys

    #once you know the time of the obs calc the average power
    avg_power = get_beam_power([obsid,ra_obs, dec_obs, time_detection, delays ,centrefreq, channels],
                                enter, exit, [[pul_ra, pul_dec]], option = 'e')
    print "Average power (fraction)" + str(avg_power)

    #gain calc

    #gain for incoherent beam using a forumla from oransaye 2015
    obswav = constants.c / (centrefreq)
    gain = math.pow(obswav,2) * 4. *math.sqrt(128.) * avg_power\
         / (constants.k * constants.pi * math.pow(10,26) )
    print "Old gain val:",gain
    """


    #New gain calc
    sys_temp, gain_table, avg_power = get_Tsys.main(obsid,[pul_ra, pul_dec],enter,exit,incoh)
    gain = np.mean(gain_table)
    t_sys = np.mean(sys_temp)

    
    
    #TODO move the calc into a def
    #ra_obs dec_obs obsid profile num_bins pulsar
    
    #gain uncertainty through beam position estimates
    mwa = ephem_utils.Obs[ephem_utils.obscode['MWA']]    
    RAs, Decs = sex2deg(ra_obs,dec_obs)
    obstime = Time(float(obsid),format='gps',scale='utc')
    observer = ephem.Observer()
    observer.date = obstime.datetime.strftime('%Y/%m/%d %H:%M:%S')
    LST_hours = observer.sidereal_time() * ephem_utils.HRS_IN_RADIAN

    HAs = -RAs + LST_hours * 15.
    Azs, Alts = ephem_utils.eq2horz(HAs, Decs, mwa.lat)
    theta=np.radians(90.-Alts)

    u_gain_per = (1. - avg_power)*0.12 + (theta/90.)*(theta/90.)*2. + 0.1
    u_gain = gain * u_gain_per #assumed to be 10% 
        
    #then centers it around the max flux (hopefully the centre of the pulse
    shiftby=int(np.argmax(profile))-int(num_bins)/2 
    if shiftby < 0.0:
        shiftby = shiftby + int(num_bins)
    profile = np.append(profile[shiftby:], profile[:shiftby])

    #plot the profile so the pulse width can be determined by eye
    print 'Please examine the plot by eye to determine the pulse width'
    plt.plot(profile)
    plt.axis([0, num_bins, 0, max(profile)])
    plt.title(pulsar)
    plt.show(block=False)
    #width min and max input
    min_bin = raw_input("Input the first bin of the pulse:  ")
    max_bin = raw_input("Input the last bin of the pulse:  ")
    plt.close()
    #calc pulse width
    pulse_width_bins = float(max_bin)-float(min_bin)+1.0 #in bins
    off_pulse_width_bins = float(num_bins)-pulse_width_bins
    pulse_width_frac = pulse_width_bins/float(num_bins)


    #calc off-pulse average
    off_pulse_total = 0.
    for o in range(len(profile)):
        if not (int(min_bin) < o < int(max_bin)):
            off_pulse_total = off_pulse_total + float(profile[o])
    off_pulse_mean = off_pulse_total / off_pulse_width_bins

    #adjust profile to be around the off-pulse mean
    for p in range(len(profile)):
        profile[p] = profile[p] - off_pulse_mean

    profile_uncert = 500. #uncertainty approximation as none is given
        
    #calc off-pulse sigma and it's uncertainty
    sigma_total = 0.
    u_sigma_total = 0.
    for o in range(len(profile)):
        if not (int(min_bin) < o < int(max_bin)):
            sigma_total = sigma_total + math.pow(float(profile[o]),2)
            u_sigma_total = u_sigma_total + 2 * float(profile[o]) * profile_uncert
    sigma = math.sqrt(abs(sigma_total / off_pulse_width_bins))
    u_sigma = profile_uncert / (2 * math.sqrt(abs(sigma_total * off_pulse_width_bins)))

    #calc area under the pulse in it's input arbritatry units from here on denoted as p.
    p_total = 0.
    u_p_total = 0. #p uncertainty
    for p in range(len(profile)):
        if int(min_bin) <= p <= int(max_bin):
            p_total = p_total + profile[p]
            u_p_total = math.sqrt(math.pow(u_p_total,2) + math.pow(profile_uncert,2))
            
    #the equivalent width (assumes top hat pulsar) in bins and it's uncertainty
    w_equiv_bins = p_total / max(profile)
    w_equiv_ms = w_equiv_bins / float(num_bins) * float(period) # in ms
    u_w_equiv_bins = math.sqrt(math.pow(u_p_total / max(profile),2) + \
                               math.pow(p_total * profile_uncert / math.pow(max(profile),2),2))
    u_w_equiv_ms = u_w_equiv_bins / float(num_bins) * float(period) # in ms
                               
    #calc signal to noise ration and it's uncertainty
    sn = p_total / (pulse_width_bins * sigma)
    u_sn = math.sqrt( math.pow( u_p_total / sigma , 2)  +  math.pow( p_total * u_sigma / math.pow(sigma,2) ,2) )

    #final calc of the mean fluxdesnity
    S_mean = sn * t_sys / ( gain * math.sqrt(2. * float(time_detection) * bandwidth)) *\
             math.sqrt( w_equiv_bins / (num_bins - w_equiv_bins)) * 1000.
    #constants to make uncertainty calc easier
    S_mean_cons = t_sys / ( math.sqrt(2. * float(time_detection) * bandwidth)) *\
             math.sqrt( w_equiv_bins / (num_bins - w_equiv_bins)) * 1000. 
    u_S_mean = math.sqrt( math.pow(S_mean_cons * u_sn / gain , 2)  +\
                          math.pow(sn * S_mean_cons * u_gain / math.pow(gain,2) , 2) )  

    print "SN " + str(sn)
    #print "Sky temp " + str(sky_temp) + " K"
    print "T_sys " + str(t_sys) + " K"
    print "Gain " + str(gain) + " K/Jy"
    print "Smean " + str(S_mean) + ' +/- ' + str(u_S_mean) + ' mJy'

    #calc obstype
    if (maxfreq - minfreq) == 23:
        obstype = 1
    else:
        obstype = 2
        
    #calc scattering 
    scat_height = max(profile) / 2.71828
    scat_bins = 0
    for p in profile:
        if p > scat_height:
            scat_bins = scat_bins + 1
    scattering = float(scat_bins + 1) * float(period) /1000. #in s
    u_scattering = 1. * float(period) /1000. # assumes the uncertainty is one bin

    #calc sub-bands
    subbands = 1
    for b in range(len(channels)):
        if b == 0:
            continue
        if not (channels[b] - channels[b-1]) == 1:
            subbands = subbands + 1
    
    #get cal id
    cal_list = client.calibrator_list(web_address, auth)
    cal_already_created = False
    for c in cal_list:
        if ( c[u'observationid'] == int(obsid) ) and ( c[u'caltype'] == calibrator_type ):
            cal_already_created = True
            cal_db_id = c[u'id']
    if not cal_already_created:
        cal_db_id = client.calibrator_create(web_address, auth,
                                              observationid = str(args.cal_id),
                                              caltype = calibrator_type)[u'id']
    
    try:
        client.detection_create(web_address, auth, 
                               observationid = int(obsid),
                               pulsar = str(pulsar), 
                               subband = str(subbands), 
                               incoherent = incoh,
                               observation_type = int(obstype),
                               calibrator = int(cal_db_id),
                               startcchan = int(minfreq), stopcchan = int(maxfreq), 
                               flux = float("{0:.2f}".format(S_mean)),
                               flux_error = float("{0:.2f}".format(u_S_mean)),
                               width = float("{0:.2f}".format(w_equiv_ms)),
                               width_error = float("{0:.2f}".format(u_w_equiv_ms)),
                               scattering = float("{0:.5f}".format(scattering)), 
                               scattering_error = float("{0:.5f}".format(u_scattering)),
                               dm = float(dm))
    except:
        print "Detection already on database so updating the values"
        client.detection_update(web_address, auth, 
                               observationid = int(obsid),
                               pulsar = str(pulsar), 
                               subband = str(subbands), 
                               incoherent = incoh,
                               observation_type = int(obstype),
                               calibrator = int(cal_db_id),
                               startcchan = int(minfreq), stopcchan = int(maxfreq), 
                               flux = float("{0:.2f}".format(S_mean)),
                               flux_error = float("{0:.2f}".format(u_S_mean)),
                               width = float("{0:.2f}".format(w_equiv_ms)),
                               width_error = float("{0:.2f}".format(u_w_equiv_ms)),
                               scattering = float("{0:.5f}".format(scattering)), 
                               scattering_error = float("{0:.5f}".format(u_scattering)),
                               dm = float(dm))
                           
    print "Observation submitted to database"



if args.pulsar and not args.bp_file:  
    #uploads files to database if there's the no calc option
    #checks if the observation is on the database
    try:
        temp_dict = client.detection_get(web_address, auth, observationid = str(obsid))
    except:
        #no obsid so creats a blank one and assumes the subbands are continuous
        client.detection_create(web_address, auth, 
                                observationid = str(obsid),
                                pulsar = str(pulsar),
                                subband = 1,
                                incoherent = incoh,
                                observation_type = '1')  
        temp_dict = client.detection_get(web_address, auth, str(obsid))  
        subbands=1
    
    
    for t in range(len(temp_dict)):
        if pulsar == temp_dict[t][u'pulsar']:
            subbands = temp_dict[t][u'subband']

#Archive files
if args.archive:
    print "Uploading archive file to database"
    client.detection_file_upload(web_address, auth, 
                                observationid = str(obsid),
                                pulsar = str(pulsar), 
                                subband = int(subbands),
                                incoherent = incoh,
                                filetype = 1,
                                filepath = str(args.archive))

if args.timeseries:
    print "Uploading timeseries file to database"
    client.detection_file_upload(web_address, auth,
                                observationid = str(obsid),
                                pulsar = str(pulsar), 
                                subband = int(subbands),
                                incoherent = incoh,
                                filetype = 2,
                                filepath = str(args.timeseries))

if args.diagnostics:
    diag = args.diagnostics
    for i in range(3):
        if diag[i] != 'None':
            if i == 0:
                cp(str(args.diagnostics[0]) ,str(obsid) + "_" + str(pulsar) + ".prepfold.ps")
                d_file_loc = str(obsid) + "_" + str(pulsar) + ".prepfold.ps"
            if i == 1:
                cp(str(args.diagnostics[1]),str(obsid) + "_" + str(pulsar) + ".prof.ps")
                d_file_loc = str(obsid) + "_" + str(pulsar) + ".prof.ps"
            if i ==2:
                cp(str(args.diagnostics[2]),str(obsid) + "_" + str(pulsar) + ".freq.vs.phase.ps")
                d_file_loc = str(obsid) + "_" + str(pulsar) + ".freq.vs.phase.ps"
            print "Uploading diagnostic file to database"
            client.detection_file_upload(web_address, auth, 
                                observationid = str(obsid),
                                pulsar = str(pulsar), 
                                subband = int(subbands),
                                incoherent = incoh,
                                filetype = 3,
                                filepath = str(d_file_loc))
            os.system("rm " + d_file_loc)
    
if args.calibration:
    cal_list = client.calibrator_list(web_address, auth)
    cal_already_created = False
    for c in cal_list:
        if ( c[u'observationid'] == int(args.cal_id) ) and ( c[u'caltype'] == calibrator_type ):
            cal_already_created = True
            cal_db_id = c[u'id']
    if not cal_already_created:
        cal_db_id = client.calibrator_create(web_address, auth,
                                              observationid = str(args.cal_id),
                                              caltype = calibrator_type)#[u'id']
    
    if args.andre:
        cp(str(args.calibration),str(args.cal_id) + "_" + str(obsid) + ".bin")
        cal_file_loc = str(args.cal_id) + "_" + str(obsid) + ".bin"
    else:
        cp(str(args.calibration),str(args.cal_id) + "_" + str(obsid) + ".tar")
        cal_file_loc = str(args.cal_id) + "_" + str(obsid) + ".tar"
    
    print "Uploading calibration solution to database"
    client.calibrator_file_upload(web_address, auth, 
                                   observationid = str(args.cal_id),
                                   caltype = calibrator_type, 
                                   filepath = str(cal_file_loc))
    os.system("rm " + cal_file_loc)
    
    #result = client.detection_find_calibrator(web_address, auth,detection_obsid = 35)
    
    
                                
    


if args.u_archive or args.u_timeseries or args.u_diagnostics:
    #runs all needed jobs to create all files
    with open(str(obsid) + '_' + str(pulsar) + '.batch','w') as batch_file:
        batch_line = "#!/bin/bash -l\n" +\
                     "#SBATCH --job-name=submit_to_databse\n" +\
                     "#SBATCH --output=" + str(obsid) + '_' + str(pulsar) + ".out\n" +\
                     "#SBATCH --export=NONE\n" +\
                     "#SBATCH --partition=workq\n" +\
                     "#SBATCH --time=6:50:00\n" +\
                     "#SBATCH --gid=mwaops\n" +\
                     "#SBATCH --account=mwaops\n" +\
                     "#SBATCH --nodes=1\n" +\
                     "ncpus=20\n"+\
                     "export OMP_NUM_THREADS=$ncpus\n"  +\
                     "psrcat -e " + str(pulsar) + " > " +str(pulsar) + ".par\n" +\
                     "fits=(" + fits_files_loc + ")\n"
        batch_file.write(batch_line)            
        if args.archive:
            batch_line = "ar_loc=" + str(args.archive)[:-3] + "\n"
            batch_file.write(batch_line)
        elif args.timeseries:
            batch_line = "ar_loc=" + str(args.timeseries)[:-3] + "\n"
            batch_file.write(batch_line)
        elif args.u_diagnostics or args.u_archive:
            batch_line = "aprun -b -cc none -d $ncpus dspsr -E " +\
                            str(pulsar) + ".par -b " + str(num_bins) + " -A -cont -O "+\
                            str(obsid) + "_" + str(pulsar) + " " + "${fits}\n" +\
                         'psraddstring="psradd -o '+ str(obsid) + "_" + str(pulsar) + '.ar "\n' +\
                         'for ((i=0;i<${#fits[@]};i++)); do psraddstring=${psraddstring}" "' +\
                            str(obsid) + "_" + str(pulsar) + '_$(expr $i).ar ; done\n' +\
                         "aprun -b -cc none -d $ncpus $psraddstring\n" +\
                         "rm " + str(obsid) + "_" + str(pulsar) + "_[0-9].ar\n" +\
                         "rm " + str(obsid) + "_" + str(pulsar) + "_[0-9][0-9].ar\n" +\
                         "ar_loc=" + str(obsid) + "_" + str(pulsar) + "\n"
            batch_file.write(batch_line)
        
        if args.u_diagnostics:
            batch_line = "pav -CDFTp -N1,1 -g " + str(obsid) + "_" + str(pulsar) + ".prof.ps/cps " +\
                             "${ar_loc}.ar\n" +\
                         'psrplot -pG -jCDTp -j "B ' +str(num_bins) + '" -D '+ str(obsid) + "_" +\
                             str(pulsar) + ".freq.vs.phase.ps/cps ${ar_loc}.ar\n" 
            batch_file.write(batch_line)
        if args.u_timeseries:
            batch_line = "aprun -b -cc none -d $ncpus dspsr -E " + str(pulsar) + ".par -b " + str(num_bins)\
                            + " -cont -s -K ${fits}\n" +\
                         'psraddstring="psradd -o '+ str(obsid) + "_" + str(pulsar) + '.ts.ar "\n' +\
                         "ts=(pulse*.ar)\n" +\
                         'for ((i=0;i<${#ts[@]};i++)); do psraddstring=${psraddstring}" "' +\
                            '${ts[i]} ; done\n' +\
                         "aprun -b -cc none -d $ncpus $psraddstring\n" +\
                         "rm pulse*.ar\n" 
            batch_file.write(batch_line)
        
    submit_line = 'sbatch ' + str(obsid) + '_' + str(pulsar) + '.batch'
    submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
    #TODO upload to database automatically once I'm confident and work out if there's an easier way then just running the script once the jobs are done, maybe add some checks


"""
Program tests:
python submit_to_database.py  ../1120799848/fold/1120799848_PSR_0837-4135.pfd.bestprof
python submit_to_database.py --bp_file /group/mwaops/incoh_census_psr/PSR_Prof/1121173352nsch_PSR_1534-5334.pfd.bestprof  --u_archive --u_timeseries --u_diagnostics
python submit_to_database.py --bp_file ../PSR_Prof/1152636328_PSR_1943-1237.pfd.bestprof  --u_archive --u_timeseries --u_diagnostics --incoh
#^ no fits
python submit_to_database.py --bp_file ../PSR_Prof/1139324488_PSR_0837+0610.pfd.bestprof --u_archive --u_timeseries --u_diagnostics --incoh

"""
