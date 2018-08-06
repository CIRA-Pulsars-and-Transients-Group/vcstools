#! /usr/bin/env python 

"""
Author: Nicholas Swainston
Creation Date: /05/2016

The MWA Pulsar Database was created by David Pallot and he wrote the database call functions.

This code is used to submit observations to the MWA Pulsar Database. The goal is to calculate all need values (flux density, width, scattering) for each observation and  submit them to the pulsar database without having to manually input them.
"""

__author__ = 'Nicholas Swainston'
__date__ = '2016-05-12'


import os
import argparse
import numpy as np
import subprocess
import sys
from shutil import copyfile as cp
import math
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import glob
from astropy.table import Table
from astropy.time import Time

import ephem
from mwapy import ephem_utils

#TODO check if I need the below two imports
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()

from mwapy.pb import primary_beam
from mwapy.pb import primarybeammap_tant as pbtant
import mwapy.pb.primarybeammap as pbl
from mwa_pulsar_client import client
import mwa_metadb_utils as meta
import find_pulsar_in_obs

web_address = 'mwa-pawsey-volt01.pawsey.ivec.org'
auth = ('mwapulsar','veovys9OUTY=')

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
    import requests
    path = 'https://{0}/{1}/'.format(addr, 'psrcat')
    payload = {'name': pulsar, 'format': 'json'}
    r = requests.post(url = path, auth = auth, data = payload)
    r.raise_for_status()
    return r.json()

    
def sex2deg( ra, dec):
    """
    sex2deg( ra, dec)
    ra - the right ascension in HH:MM:SS
    dec - the declination in DD:MM:SS
    
    Convert sexagesimal coordinates to degrees.
    """ 
    from astropy.coordinates import SkyCoord
    from astropy import units as u
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
    from mwapy.pb import primary_beam
    from astropy.time import Time
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
    avg_power = np.mean(Powers)
    return avg_power

def get_pulsar_ra_dec(pulsar):
    #Gets the ra and dec from the output of PSRCAT
    cmd = ['psrcat', '-c', 'Raj', pulsar]
    output = subprocess.Popen(cmd,stdout=subprocess.PIPE).communicate()[0]
    temp = []
    lines = output.split('\n')
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


def from_power_to_gain(powers,cfreq,n,incoh=False):
    from astropy.constants import c,k_B
    from math import sqrt

    obswl = c.value/cfreq
    #for incoherent
    if incoh:
        coeff = obswl**2*16*sqrt(n)/(4*np.pi*k_B.value) #TODO add a coherent option
    else:
        coeff = obswl**2*16*n/(4*np.pi*k_B.value)
    print "Wavelength",obswl,"m"
    print "Gain coefficient:",coeff
    SI_to_Jy = 1e-26
    return (powers*coeff)*SI_to_Jy


def get_Trec(tab,obsfreq):
    Trec = 0.0
    for r in range(len(tab)-1):
        if tab[r][0]==obsfreq:
            Trec = tab[r][1]
        elif tab[r][0] < obsfreq < tab[r+1][0]:
            Trec = ((tab[r][1] + tab[r+1][1])/2)
    if Trec == 0.0:
        print "ERROR getting Trec"
    return Trec


class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()  
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


def sigmaClip(data, alpha=3, tol=0.1, ntrials=10):
    x = np.copy(data)
    oldstd = np.nanstd(x)
    
    for trial in range(ntrials):
        median = np.nanmedian(x)

        lolim = median - alpha * oldstd
        hilim = median + alpha * oldstd
        x[x<lolim] = np.nan
        x[x>hilim] = np.nan

        newstd = np.nanstd(x)
        tollvl = (oldstd - newstd) / newstd

        if tollvl <= tol:
            print "Took {0} trials to reach tolerance".format(trial+1)
            return oldstd, x

        if trial == ntrials:
            print "Reached number of trials without reaching tolerance level"
            return oldstd, x

        oldstd = newstd


def enter_exit_calc(time_detection, time_obs, metadata, start = None, stop = None):
    """
    time_detection: the time in seconds of the dectection from the bestprof file
    time_obs: the time in seconds of the dectection from the metadata
    start: input start time of the obs (0 is the true start)
    end: input end time of the obs
    metadata: list from the function get_obs_metadata
    """
    obsid,ra_obs,dec_obs,time_obs,delays,centrefreq,channels = metadata
    #check if the obs time is entire obs. The meta data will round down to the nearest 200 seconds
    #(someitmmes 100 depending on the obs type) 
    if 0. < (time_detection - time_obs) < 199.:
        #all available fits files used
        enter = 0.
        exit = time_detection
    else:
        #check if there is an input start and stop time
        if time_detection > (time_obs *5.):
            #some of the comissioning data has terrible metadata 
            #going to assume it's for the entire data set
            enter = 0.
            exit = float(time_detection)
        #start == None means no default used
        elif (start == None) and (stop is not None):
            print "No start time input so assumed it's 0"
            enter = 0.
            exit = stop
            exit = float(exit)
        elif (start is not None) and (stop) == None:
            print "No stop time input so assumed it's the end of the observation"
            exit = float(time_detection)
            enter = args.start
            enter = float(enter)
        elif (start is not None) and (stop is not None):
            enter = start
            exit = stop
            enter = float(enter)
            exit = float(exit)
        else:
            #find_pulsar_in_obs wrapping to use it to find start and end
            find_pulsar_in_obs.grab_pulsaralog([args.pulsar])
            catDIR = 'temp.csv'
            catalog = Table.read(catDIR)
            enter, exit = find_pulsar_in_obs.get_beam_power([obsid,ra_obs,dec_obs,
                                           time_obs,delays, centrefreq,channels],
                                           catalog)
            enter *= float(time_obs) 
            exit *= float(time_obs)
            os.remove('temp.csv')
            if os.path.exists("{0}_analytic_beam.txt".format(obsid)):
                os.remove("{0}_analytic_beam.txt".format(obsid))
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
    return enter, exit



def flux_cal_and_sumbit(time_detection, time_obs, metadata, bestprof_data,
                        pul_ra, pul_dec, incoh,
                        start = None, stop = None,
                        trcvr = "/group/mwaops/PULSAR/MWA_Trcvr_tile_56.csv"):
    """
    time_detection: the time in seconds of the dectection from the bestprof file
    time_obs: the time in seconds of the dectection from the metadata
    start: input start time of the obs (0 is the true start)
    end: input end time of the obs
    metadata: list from the function get_obs_metadata
    bestprof_data: list from the function get_from_bestprof
    trcvr: the file location of antena temperatures
    """
    #calculate the start and stop time of the detection
    enter, exit = enter_exit_calc(time_detection, time_obs, metadata, start = None, stop = None)
    obsdur = enter - exit

    #unpack data
    obsid, pulsar, dm, period, period_uncer, time_detection, profile, num_bins = bestprof_data
    obsid,ra_obs,dec_obs,time_obs,delays,centrefreq,channels = metadata
    
    
    #Gain calc
    #get antena temperatures    
    trec_table = Table.read(trcvr,format="csv")
    
    ntiles = 128#TODO actually we excluded some tiles during beamforming, so we'll need to account for that here
    
    print "Calculating beam power and antena temperature..."
    bandpowers = get_beam_power(metadata,enter,exit, [[pul_ra, pul_dec]],
                                centeronly=False,dt=100,option="a") 
    
    beamsky_sum_XX,beam_sum_XX,Tant_XX,beam_dOMEGA_sum_XX,\
     beamsky_sum_YY,beam_sum_YY,Tant_YY,beam_dOMEGA_sum_YY =\
     pbtant.make_primarybeammap(obsid, delays, centrefreq*1e6, 'analytic', plottype='None')
    
    #TODO can be inaccurate for coherent but is too difficult to simulate
    tant = (Tant_XX + Tant_YY) / 2.

    print "Tant: " + str(tant)
    t_sys_table = tant + get_Trec(trec_table,centrefreq)
    
    print "Converting to gain from power..."
    gain = from_power_to_gain(bandpowers,centrefreq*1e6,ntiles,incoh)
    #print 'Frequency',centrefreq*1e6,'Hz'
    
    t_sys = np.mean(t_sys_table)
    avg_power = np.mean(bandpowers)
    print "Average Power: " + str(avg_power)

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
    shiftby=-int(np.argmax(profile))+int(num_bins)/2 
    profile = np.roll(profile,shiftby)

    sigma, flagged_profile  = sigmaClip(profile, alpha=3, tol=0.05, ntrials=10)
    """
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
    """

    #adjust profile to be around the off-pulse mean
    off_pulse_mean = np.nanmean(flagged_profile)   
    profile -= off_pulse_mean
    flagged_profile -= off_pulse_mean

    profile_uncert = 500. #uncertainty approximation as none is given

    #calc off-pulse sigma uncertainty
    pulse_width_bins = 1
    off_pulse_width_bins = 1
    p_total = 0.
    u_p_total = 0. #p uncertainty
    sigma_total = 0.
    u_sigma_total = 0.
    for p in range(len(profile)):
        if math.isnan(flagged_profile[p]):
            #May increase the pulse width if there are large noise spikes
            pulse_width_bins += 1    
            p_total = p_total + profile[p]
            u_p_total = math.sqrt(math.pow(u_p_total,2) + math.pow(profile_uncert,2))
        else:
            off_pulse_width_bins += 1
            sigma_total = sigma_total + math.pow(float(profile[p]),2)
            u_sigma_total = u_sigma_total + 2 * float(profile[p]) * profile_uncert

    u_sigma = profile_uncert / (2 * math.sqrt(abs(sigma_total * off_pulse_width_bins)))
    profile_uncert = 500. #uncertainty approximation as none is given
        
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
    if not incoh:
        cal_list = client.calibrator_list(web_address, auth)
        cal_already_created = False
        for c in cal_list:
            if ( c[u'observationid'] == int(args.cal_id) ) and ( c[u'caltype'] == calibrator_type ):
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
    else:
        #submits without the cal_id
        try:
            client.detection_create(web_address, auth, 
                                   observationid = int(obsid),
                                   pulsar = str(pulsar), 
                                   subband = str(subbands), 
                                   incoherent = incoh,
                                   observation_type = int(obstype),
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
                                   startcchan = int(minfreq), stopcchan = int(maxfreq), 
                                   flux = float("{0:.2f}".format(S_mean)),
                                   flux_error = float("{0:.2f}".format(u_S_mean)),
                                   width = float("{0:.2f}".format(w_equiv_ms)),
                                   width_error = float("{0:.2f}".format(u_w_equiv_ms)),
                                   scattering = float("{0:.5f}".format(scattering)), 
                                   scattering_error = float("{0:.5f}".format(u_scattering)),
                                   dm = float(dm))
                               
        print "Observation submitted to database"
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    This code is used to submit observations to the MWA Pulsar Database. The goal is to calculate all need values (flux density, width, scattering) for each observation and  submit them to the pulsar database without having to manually input them. It can also submit diagnostic files such as pulse profiles and calibrations. This code assumes that the observation is coherent and correlated using the RTS so please use --incoh and --andre if your observation is incoherent or correlated using Andre's Tools respectively.
    Two common uses are to simply upload a tarred calibration file (this can be done before a detection is uploaded)\n
    > submit_to_database.py -c <calibration tar file location> -o <obsid> --cal_id <calibrator obsid>\n
    Another is to upload the pulsars parameters, such as flux density and width, to the database. To do this the bestprof file is needed \n
    > submit_to_database.py -b <bestprof file location>\n
    An example of uploading a diagnostic file is:
    > submit_to_database.py  -o <obsid> -p <pulsar> --ppps <PRESTO prepfold output post script file location> --cal_id <calid>
    """, formatter_class=SmartFormatter)
    parser.add_argument('-o','--obsid',type=str,help='The obsid.')
    parser.add_argument('--incoh',action='store_true',help='Used for incoherent detections to accurately calculate gain. Default is coherent.')
    parser.add_argument('--andre',action='store_true',help="Used for calibrations done using Andre Offrina's tools. Default is RTS.")
    parser.add_argument('-p','--pulsar',type=str,help='The pulsar J name.')
    parser.add_argument('-f','--fits_files',type=str,help='The fits files location to be used in any Detection Calculation processing. Recommended to end in *.fits and surrounded by quotation marks.', default = "/group/mwaops/vcs/${obsid}/fits/*fits")
    parser.add_argument('--cal_id',type=str,help='The obsid of the calibrator.')

    calcargs = parser.add_argument_group('Detection Calculation Options', 'All the values of a pulsar detection for the MWA pulsar databse can be calculated by this code using either a .besprof file or a (Soon to be implimented dspsr equivalent). If neither are included then all values will be left as null and the observation(/detection) can still be used to upload files without the calculation being performed')
    calcargs.add_argument('-b','--bestprof',type=str,help='The location of the .bestprof file. Using this option will cause the code to calculate the needed paramters to be uploaded to the database (such as flux density, width and scattering). Using this option can be used instead of inputing the obsid and pulsar.')
    calcargs.add_argument('--start',type=str,help='The start time of the detection in seconds. For example: 0', default = None)
    calcargs.add_argument('--stop',type=str,help='The stop time of the detection in seconds. For example: 1000', default = None)
    calcargs.add_argument('--trcvr',type=str,help='File location of the reciever temperatures to be used',default = "/group/mwaops/PULSAR/MWA_Trcvr_tile_56.csv")

    uploadargs = parser.add_argument_group('Upload options', 'The different options for each file type that can be uploaded to the pulsar database. Will cause an error if the wrong file type is being uploaded.')
    uploadargs.add_argument('-a','--archive',type=str,help="The dspsr archive file location to be uploaded to the database. Expects a single file that is the output of dspsr using the pulsar's ephemeris.")
    uploadargs.add_argument('--single_pulse_series',type=str,help='The single pulse series file location to be uploaded to the database. Expects a single file that is the output of dspsr in single pulse mode (the -s option).')
    uploadargs.add_argument('--ppps',type=str,help="The Presto Prepfold PostScript file location to be uploaded to the database.")
    uploadargs.add_argument('-i','--ippd',type=str,help="The Intergrates Pulse Profile given by Dspsr.")
    uploadargs.add_argument('-w','--waterfall',type=str,help="A waterfall plot of pulse phase vs frequency using dspsr's psrplot.")
    uploadargs.add_argument('-c','--calibration',type=str,help='The calibration solution file location to be uploaded to the database. Expects a single file so please zip or tar up the bandpass calibrations, the DI Jones matrices, the flagged_channels.txt file, the flagged_tiles.txt file, the rts.in file and the source file.')

    dspsrargs = parser.add_argument_group('Dspsr Calculation Options', "Requires the --fits_files. These options will send off dspsr jobs to process the needed files that can be uploaded to the database. The files will be uploaded automatically when the dspsr scripts are tested more") #TODO remove when I'm confident with dspsr
    dspsrargs.add_argument('--u_archive', action='store_true',help='The archive to be processed with dspsr and then uploaded to the database')
    dspsrargs.add_argument('--u_single_pulse_series', action='store_true',help='The single pulse series to be processed with dspsr and then uploaded to the database')
    dspsrargs.add_argument('--u_ppps', action='store_true', help="The Presto Prepfold PostScript file will be process by sending off PRESTO jobs.")
    dspsrargs.add_argument('--u_ippd', action='store_true', help="The Intergrates Pulse Profile file will be process by sending off DSPSR jobs.")
    dspsrargs.add_argument('--u_waterfall', action='store_true', help="A waterfall plot of pulse phase vs frequency file will be process by sending off DSPSR jobs.")
    args=parser.parse_args()



    #defaults for incoh and calibrator type
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

    if args.fits_files:
        fits_files_loc = args.fits_files
    else:
        fits_files_loc = '/scratch2/mwaops/vcs/'+str(obsid)+'/fits/*.fits'

    #get info from .bestprof file
    if args.bestprof:
        bestprof_data = get_from_bestprof(args.bestprof)
        obsid, pulsar, dm, period, period_uncer, time_detection, profile, num_bins = bestprof_data
    elif args.obsid and args.pulsar:
        num_bins = 128 #used in dspsr calculations
        obsid = args.obsid
        pulsar = args.pulsar
    elif args.obsid and args.calibration and not (args.archive or args.single_pulse_series or args.ppps or args.ippd  or args.waterfall or args.u_archive or args.u_single_pulse_series or args.u_ppps or args.u_ippd  or args.u_waterfall):
        obsid = args.obsid
    else:
        print "Please us either --obsid and --pulsar or --bestprof"
        sys.exit(0)
        
    if args.pulsar or args.bestprof:
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
    metadata = meta.get_common_obs_metadata(obsid)
    obsid,ra_obs,dec_obs,time_obs,delays,centrefreq,channels = metadata
    minfreq = float(min(channels))
    maxfreq = float(max(channels))

    bandwidth = 30720000. #In Hz

    #calc obstype
    if (maxfreq - minfreq) == 23:
        obstype = 1
    else:
        obstype = 2


    if args.bestprof:
        #Does the flux calculation and submits the results to the MWA pulsar database
        flux_cal_and_sumbit(time_detection, time_obs, metadata, bestprof_data,
                            pul_ra, pul_dec, incoh,
                            start = args.start, stop = args.stop, trcvr = args.trcvr)

    if args.pulsar and not args.bestprof:  
        #calc sub-bands
        subbands = 1
        for b in range(len(channels)):
            if b == 0:
                continue
            if not (channels[b] - channels[b-1]) == 1:
                subbands = subbands + 1
                
        if not incoh:
            cal_list = client.calibrator_list(web_address, auth)
            cal_already_created = False
            for c in cal_list:
                if ( c[u'observationid'] == int(args.cal_id) ) and ( c[u'caltype'] == calibrator_type ):
                    cal_already_created = True
                    cal_db_id = c[u'id']
            if not cal_already_created:
                cal_db_id = client.calibrator_create(web_address, auth,
                                                      observationid = str(args.cal_id),
                                                      caltype = calibrator_type)[u'id']
        elif not args.cal_id:
            cal_db_id = None
            
        #uploads files to database if there's the no calc option
        #checks if the observation is on the database
        try:
            temp_dict = client.detection_get(web_address, auth, observationid = str(obsid))
        except:
            client.detection_create(web_address, auth, 
                                    observationid = str(obsid),
                                    pulsar = str(pulsar),
                                    calibrator = int(cal_db_id),
                                    subband = int(subbands),
                                    incoherent = incoh,
                                    startcchan = int(minfreq), stopcchan = int(maxfreq), 
                                    observation_type = int(obstype))  
            temp_dict = client.detection_get(web_address, auth, observationid =str(obsid))  
        
        if not temp_dict:
            #no obsid so creats a blank one and assumes the subbands are continuous
            client.detection_create(web_address, auth, 
                                    observationid = str(obsid),
                                    pulsar = str(pulsar),
                                    calibrator = int(cal_db_id),
                                    subband = int(subbands),
                                    incoherent = incoh,
                                    startcchan = int(minfreq), stopcchan = int(maxfreq), 
                                    observation_type = int(obstype))  
            temp_dict = client.detection_get(web_address, auth, observationid = str(obsid)) 
        
        pulsar_dict_check = False
        for t in range(len(temp_dict)):
            if pulsar == temp_dict[t][u'pulsar']:
                subbands = temp_dict[t][u'subband']
                pulsar_dict_check = True
        
        if not pulsar_dict_check:
            client.detection_create(web_address, auth, 
                                    observationid = str(obsid),
                                    pulsar = str(pulsar),
                                    calibrator = int(cal_db_id),
                                    subband = int(subbands),
                                    incoherent = incoh,
                                    startcchan = int(minfreq), stopcchan = int(maxfreq), 
                                    observation_type = int(obstype))  
            temp_dict = client.detection_get(web_address, auth, observationid = str(obsid))  

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

    if args.single_pulse_series:
        print "Uploading single_pulse_series file to database"
        client.detection_file_upload(web_address, auth,
                                    observationid = str(obsid),
                                    pulsar = str(pulsar), 
                                    subband = int(subbands),
                                    incoherent = incoh,
                                    filetype = 2,
                                    filepath = str(args.single_pulse_series))

    if args.ppps:
        cp(str(args.ppps) ,str(obsid) + "_" + str(pulsar) + ".prepfold.ps")
        d_file_loc = str(obsid) + "_" + str(pulsar) + ".prepfold.ps"
        print "Uploading Presto Prepfold PostScript file to database"
        client.detection_file_upload(web_address, auth, 
                            observationid = str(obsid),
                            pulsar = str(pulsar), 
                            subband = int(subbands),
                            incoherent = incoh,
                            filetype = 3,
                            filepath = str(d_file_loc))
        os.system("rm " + d_file_loc)
        
    if args.ippd:
        cp(str(args.ippd),str(obsid) + "_" + str(pulsar) + ".prof.ps")
        d_file_loc = str(obsid) + "_" + str(pulsar) + ".prof.ps"
        print "Uploading Intergrates Pulse Profile file to database"
        client.detection_file_upload(web_address, auth, 
                            observationid = str(obsid),
                            pulsar = str(pulsar), 
                            subband = int(subbands),
                            incoherent = incoh,
                            filetype = 3,
                            filepath = str(d_file_loc))
        os.system("rm " + d_file_loc)
        
    if args.waterfall:
        cp(str(args.waterfall),str(obsid) + "_" + str(pulsar) + ".freq.vs.phase.ps")
        d_file_loc = str(obsid) + "_" + str(pulsar) + ".freq.vs.phase.ps"
        print "Uploading waterfall file to database"
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
            cp(str(args.calibration),str(args.cal_id) + "_andre_calibrator.bin")
            cal_file_loc = str(args.cal_id) + "_andre_calibrator.bin"
        else:
            cp(str(args.calibration),str(args.cal_id) + "_rts_calibrator.tar")
            cal_file_loc = str(args.cal_id) + "_rts_calibrator.tar"
        
        print "Uploading calibration solution to database"
        client.calibrator_file_upload(web_address, auth, 
                                       observationid = str(args.cal_id),
                                       caltype = calibrator_type, 
                                       filepath = str(cal_file_loc))
        os.system("rm " + cal_file_loc)
        
        #result = client.detection_find_calibrator(web_address, auth,detection_obsid = 35)
        
        
                                    
        


    if args.u_archive or args.u_single_pulse_series or args.u_ppps or args.u_ippd  or args.u_waterfall:
        if not glob.glob(fits_files_loc):
            print "No fits files in given location. Use -f to provide file location."
            quit()
        
        #runs all needed jobs to create all files
        from job_submit import submit_slurm
        dspsr_batch = "{0}_{1}".format(obsid, pulsar)
        commands = []
        n_omp_threads = 20
        commands.append("ncpus={0}".format(n_omp_threads))
        commands.append("export OMP_NUM_THREADS={0}".format(n_omp_threads))
        commands.append("psrcat -e {0} > {0}.par".format(pulsar))
        commands.append("fits=({0})".format(fits_files_loc))
        if args.archive:
            commands.append("ar_loc=" + str(args.archive)[:-3])
        elif args.single_pulse_series:
            commands.append("ar_loc=" + str(args.single_pulse_series)[:-3])
        elif args.u_ippd  or args.u_waterfall or args.u_archive:
            commands.append("srun -n 1 -c $ncpus dspsr -U 600 -E {0}.par -b {1} -A -cont -O {2}_{0} ".format(pulsar, num_bins, obsid) + "${fits}")
            commands.append('psraddstring="psradd -o {0}_{1}.ar "'.format(obsid,pulsar))
            commands.append('for ((i=0;i<${{#fits[@]}};i++)); do psraddstring=${{psraddstring}}" "{0}_{1}_$(expr $i).ar ; done'.format(obsid,pulsar))
            commands.append("srun -n 1 -c $ncpus $psraddstring")
            commands.append("rm {0}_{1}_[0-9].ar".format(obsid,pulsar))
            commands.append("rm {0}_{1}_[0-9][0-9].ar".format(obsid,pulsar))
            commands.append("ar_loc={0}_{1}".format(obsid,pulsar))
        if args.u_ippd:
            commands.append("pav -CDFTp -N1,1 -g {0}_{1}.prof.ps/cps ".format(obsid,pulsar) +\
                                "${ar_loc}.ar")
        if args.u_waterfall:
            commands.append('psrplot -pG -jCDTp -j "B {0}" -D {1}_{2}.freq.vs.phase.ps/cps '.\
                                format(num_bins,obsid,pulsar)+ "${ar_loc}.ar") 
        if args.u_ppps:
            commands.append("psrcat -e {0} > {0}.eph".format(pulsar))
            commands.append("srun -n 1 -c $ncpus prepfold -ncpus $ncpus -o {0} -topo -runavg -noclip -par {1}.eph -nsub 256 ".format(obsid, pulsar) + "${fits}")
            commands.append("{0}.eph".format(pulsar))
        if args.u_single_pulse_series:
            commands.append("srun -n 1 -c $ncpus dspsr -U 600 -E {0}.par -b {1} -cont -s -K ".\
                                format(pulsar,num_bins) + "${fits}")
            commands.append('psraddstring="psradd -o {0}_{1}.ts.ar "'.format(obsid,pulsar))
            commands.append("ts=(pulse*.ar)")
            commands.append('for ((i=0;i<${#ts[@]};i++)); do psraddstring=${psraddstring}" "${ts[i]} ; done')
            commands.append("srun -n 1 -c $ncpus $psraddstring")
            commands.append("rm pulse*.ar")
        commands.append("{0}.par".format(pulsar))
        job_id = submit_slurm(dspsr_batch, commands,
                              batch_dir="./",
                              slurm_kwargs={"time": "6:50:00", "partition": "workq"},
                              submit=True)
    """
    Program tests:
    python submit_to_database.py  ../1120799848/fold/1120799848_PSR_0837-4135.pfd.bestprof
    python submit_to_database.py --bestprof /group/mwaops/incoh_census_psr/PSR_Prof/1121173352nsch_PSR_1534-5334.pfd.bestprof  --u_archive --u_single_pulse_series --u_diagnostics
    python submit_to_database.py --bestprof ../PSR_Prof/1152636328_PSR_1943-1237.pfd.bestprof  --u_archive --u_single_pulse_series --u_diagnostics --incoh
    #^ no fits
    python submit_to_database.py --bestprof ../PSR_Prof/1139324488_PSR_0837+0610.pfd.bestprof --u_archive --u_single_pulse_series --u_diagnostics --incoh

    """
