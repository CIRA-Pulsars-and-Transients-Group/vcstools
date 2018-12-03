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
import textwrap as _textwrap
import warnings

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
import find_pulsar_in_obs as fpio

web_address = 'https://mwa-pawsey-volt01.pawsey.org.au'

class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _split_lines(self, text, width):
        text = _textwrap.dedent(self._whitespace_matcher.sub(' ', text).strip())
        return _textwrap.wrap(text, width)

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
    

def get_sigma_from_bestprof(file_loc):
    with open(file_loc,"rb") as bestprof:
        lines = bestprof.readlines()
        sigma = float(lines[13].split('~')[-1].split()[0])
        if lines[13].startswith("# Prob(Noise)") and isinstance(sigma, float):
            return sigma
        else:
            print 'Invalid sigma in {}. Exiting'.format(file_loc)
            sys.exit()


def get_from_bestprof(file_loc):
    with open(file_loc,"rb") as bestprof:
        lines = bestprof.readlines()
        #assumes the input fits files names begins with the obsid
        obsid = lines[0][22:32]
        try:
            obsid = int(obsid)
        except:
            obsid = None
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


def from_power_to_gain(powers,cfreq,n,coh=True):
    from astropy.constants import c,k_B
    from math import sqrt

    obswl = c.value/cfreq
    #for incoherent
    if coh:
        coeff = obswl**2*16*n/(4*np.pi*k_B.value)
    else:
        coeff = obswl**2*16*sqrt(n)/(4*np.pi*k_B.value)
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



def sigmaClip(data, alpha=3, tol=0.1, ntrials=10):
    x= np.copy(data)
    oldstd = np.nanstd(x)
    #removes the warning for the x[x<lolim] command when it encounters a nan
    warnings.simplefilter("ignore")
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
            warnings.simplefilter("default")
            return oldstd, x

        if trial == ntrials:
            print "Reached number of trials without reaching tolerance level"
            warnings.simplefilter("default")
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
            names_ra_dec = fpio.grab_source_alog(pulsar_list=[args.pulsar])
            source_ob_power = fpio.get_beam_power_over_time([obsid,ra_obs,dec_obs,time_obs,delays,
                                                             centrefreq,channels], names_ra_dec)
            enter, exit = fpio.beam_enter_exit(source_ob_power[0], time_obs)
            enter *= float(time_obs) 
            exit *= float(time_obs)
            input_detection_time = exit - enter
        if not int(input_detection_time) == int(time_detection):
            print "WARNING: Input detection time does not equal the dectetion time of the .bestprof file"
            print "Input (metadata) time: " + str(input_detection_time)
            print "Bestprof time: " + str(time_detection)

            #make sure calculation uses Bestprof time
            exit = enter + time_detection
    return enter, exit


def zip_calibration_files(base_dir, cal_obsid, source_file):
    """
    Checkes that all the expected calibration files are where they should be 
    and returns the file location of the zipped file
    """
    import tarfile

    #Bandpass calibrations
    bandpass_files = glob.glob("{0}/BandpassCalibration_node0*.dat".format(base_dir))
    if len(bandpass_files) != 24:
        print "Bandpass files not found. Exiting"
        quit()

    #DI Jones matricies
    DIJ_files = glob.glob("{0}/DI_JonesMatrices_node0*.dat".format(base_dir))
    if len(DIJ_files) != 24:
        print "DI Jones matricies not found. Exiting"
        quit()

    #flagged txt files
    if not ( os.path.isfile("{0}/flagged_channels.txt".format(base_dir))
            and os.path.isfile("{0}/flagged_tiles.txt".format(base_dir)) ):
        print "Flag files not found. Exiting"
        quit()

    #rts.in
    if not os.path.isfile("{0}/rts_{1}.in".format(base_dir, cal_obsid)):
        print "No rts in file. Exiting"
        quit()

    #source file
    if not os.path.isfile(source_file):
        print "No source file. Exiting"
        quit()

    #zip the files
    zip_file_location = "{0}/{1}_rts_calibrator.zip".format(base_dir, cal_obsid)
    out = tarfile.open(zip_file_location, mode='w')
    for bf in bandpass_files:
        out.add(bf, arcname=bf.split("/")[-1])
    for DIJ in DIJ_files:
        out.add(DIJ, arcname=DIJ.split("/")[-1])
    out.add("{0}/flagged_channels.txt".format(base_dir), arcname="flagged_channels.txt")
    out.add("{0}/flagged_tiles.txt".format(base_dir), arcname="flagged_tiles.txt")
    out.add("{0}/rts_{1}.in".format(base_dir, cal_obsid), arcname="rts_{0}.in".format(cal_obsid))
    out.add(source_file, arcname=source_file.split("/")[-1])
    out.close()

    return zip_file_location


def flux_cal_and_sumbit(time_detection, time_obs, metadata, 
                        bestprof_data, bestprof_loc,
                        pul_ra, pul_dec, coh, auth,
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
    enter, exit = enter_exit_calc(time_detection, time_obs, metadata, start=start, stop=stop)
    obsdur = enter - exit

    #unpack data
    obsid, pulsar, dm, period, period_uncer, time_detection, profile, num_bins = bestprof_data
    obsid,ra_obs,dec_obs,time_obs,delays,centrefreq,channels = metadata
    
    
    #Gain calc
    #get antena temperatures    
    trec_table = Table.read(trcvr,format="csv")
    
    ntiles = 128#TODO actually we excluded some tiles during beamforming, so we'll need to account for that here
    #starting from enter for the bestprof duration
    bandpowers = fpio.get_beam_power_over_time([obsid,ra_obs,dec_obs,time_detection,
                                                delays,centrefreq,channels],
                                               np.array([["source",pul_ra, pul_dec]]),
                                               dt=100, start_time=int(enter))
    bandpowers = np.mean(bandpowers)
    print "Calculating antena temperature..." 
    #represses print statements and warnings
    sys.stdout = open(os.devnull, 'w') 
    warnings.simplefilter("ignore")
    
    beamsky_sum_XX,beam_sum_XX,Tant_XX,beam_dOMEGA_sum_XX,\
     beamsky_sum_YY,beam_sum_YY,Tant_YY,beam_dOMEGA_sum_YY =\
     pbtant.make_primarybeammap(int(obsid), delays, centrefreq*1e6, 'analytic', plottype='None')
    
    #turns prints and warnings back on
    sys.stdout = sys.__stdout__
    warnings.simplefilter("default")

    #TODO can be inaccurate for coherent but is too difficult to simulate
    tant = (Tant_XX + Tant_YY) / 2.

    print "Tant: " + str(tant)
    t_sys_table = tant + get_Trec(trec_table,centrefreq)
    
    print "Converting to gain from power..."
    gain = from_power_to_gain(bandpowers,centrefreq*1e6,ntiles,coh)
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

    #sigma calc using profile, inaccuare for high SN and just used for the uncertainty calculation
    sigma, flagged_profile  = sigmaClip(profile, alpha=3, tol=0.05, ntrials=10)
    
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

    #grabbing SN from PRESTO bestprof, can be inaccurate but better than our estimates
    sn = get_sigma_from_bestprof(bestprof_loc)

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
    if coh:
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
    else:
        cal_db_id = None
 
    try:
        client.detection_create(web_address, auth, 
                               observationid = int(obsid),
                               pulsar = str(pulsar), 
                               subband = str(subbands), 
                               coherent = coh,
                               observation_type = int(obstype),
                               calibrator = int(cal_db_id),
                               startcchan = int(minfreq), stopcchan = int(maxfreq), 
                               flux = float("{0:.2f}".format(S_mean)),
                               flux_error = float("{0:.2f}".format(u_S_mean)),
                               width = float("{0:.2f}".format(w_equiv_ms)),
                               width_error = float("{0:.2f}".format(u_w_equiv_ms)),
                               scattering = float("{0:.5f}".format(scattering)), 
                               scattering_error = float("{0:.5f}".format(u_scattering)),
                               dm = float(dm),
                               version = 1)
    except:
        print "Detection already on database so updating the values"
        client.detection_update(web_address, auth, 
                               observationid = int(obsid),
                               pulsar = str(pulsar), 
                               subband = str(subbands), 
                               coherent = coh,
                               observation_type = int(obstype),
                               calibrator = int(cal_db_id),
                               startcchan = int(minfreq), stopcchan = int(maxfreq), 
                               flux = float("{0:.2f}".format(S_mean)),
                               flux_error = float("{0:.2f}".format(u_S_mean)),
                               width = float("{0:.2f}".format(w_equiv_ms)),
                               width_error = float("{0:.2f}".format(u_w_equiv_ms)),
                               scattering = float("{0:.5f}".format(scattering)), 
                               scattering_error = float("{0:.5f}".format(u_scattering)),
                               dm = float(dm),
                               version = 1)
    print "Observation submitted to database"
    return subbands


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=LineWrapRawTextHelpFormatter, description=_textwrap.dedent("""
    This code is used to submit pulsar detections and observation calibrations to the MWA Pulsar Database. It can calculate all need values (flux density, width, scattering) for each detection and submit them to the pulsar database without having to manually input them. It can also submit diagnostic files such as pulse profiles and calibrations. This code assumes that the observation is coherent and correlated using the RTS so please use --incoh and --andre if your observation is incoherent or correlated using Andre's Tools respectively.

    A common use case is to simply upload a calibration file (this can be done before a detection is uploaded)
    > submit_to_database.py --cal_dir_to_tar <calibration directory> --srclist <srclist> -o <obs ID> -O <calibrator obs ID>
    Another is to record a pulsar detection and its paramters, such as flux density and pulse width, to the database. To do this the bestprof file is needed 
    > submit_to_database.py -b <bestprof file location>
    Diagnostic files (see Upload Options) such as PRESTO plots can be uploaded using:
    > submit_to_database.py  -o <obs ID> -O <calibrator obs ID> -p <pulsar> --ppps <PRESTO prepfold output post script file location>
    The diagnostic files can be created using the Dspsr Calculation Options (not robustly tested).
    """))
    parser.add_argument('-o','--obsid',type=str,help='The observation ID (eg. 1221399680).')
    parser.add_argument('-O', '--cal_id',type=str,help='The observation ID of the calibrator.')
    parser.add_argument('-p','--pulsar',type=str,help='The pulsar J name.')
    parser.add_argument('--incoh',action='store_true',help='Used for incoherent detections to accurately calculate gain. The default is coherent.')
    parser.add_argument('--andre',action='store_true',help="Used for calibrations done using Andre Offrina's tools. Default is RTS.")
        
    calcargs = parser.add_argument_group('Detection Calculation Options', 'All the values of a pulsar detection (such as flux density, width, scattering) can be calculated by this code using either a .besprof file or a soon to be implemented DSPSR equivalent and automatically uploaded to the MWA Pulsar Database. Analysis files can still be uploaded before this step but this will leave the pulsar values as nulls.')
    calcargs.add_argument('-b','--bestprof',type=str,help='The location of the .bestprof file. Using this option will cause the code to calculate the needed parameters to be uploaded to the database (such as flux density, width and scattering). Using this option can be used instead of inputting the observation ID and pulsar name.')
    calcargs.add_argument('--start',type=str,help="The start time of the detection in seconds. For example if the detection begins at the start of the observation you would use --start 0. If this option isn't used, the code will calculate when the pulsar entered the beam", default = None)
    calcargs.add_argument('--stop',type=str,help="The stop time of the detection in seconds. For example if the detection ended 60 minutes in use --end 600. If this option isn't used, the code will calculate when the pulsar exited the beam", default = None)
    calcargs.add_argument('--trcvr',type=str, default = "/group/mwaops/PULSAR/MWA_Trcvr_tile_56.csv", help='File location of the receiver temperatures to be used. Only required if you do not want to use the default values located in %(default)s.')

    uploadargs = parser.add_argument_group('Upload Options', 'The different options for each file type that can be uploaded to the pulsar database. Will cause an error if the wrong file type is being uploaded.')
    uploadargs.add_argument('--cal_dir_to_tar',type=str,help='The calibration directory of a calibration solution that you would like to tar and upload to the database (eg. /group/mwaops/vcs/1221832280/cal/1221831856/rts). Must be used with --srclist so the correct source list is uploaded. If the calibration files are in the default positions then they will be tared and uploaded.')
    uploadargs.add_argument('--srclist',type=str,help='Used with --cal_dir to indicate the source list file location. eg /group/mwaops/vcs/1221832280/cal/1221831856/vis/srclist_pumav3_EoR0aegean_EoR1pietro+ForA_1221831856_patch1000.txt.')
    uploadargs.add_argument('-c','--calibration',type=str,help='The calibration solution file location to be uploaded to the database. Expects a single file so please zip or tar up the bandpass calibrations, the DI Jones matrices, the flagged_channels.txt file, the flagged_tiles.txt file, the rts.in file and the source file.')
    uploadargs.add_argument('-a','--archive',type=str,help="The DSPSR archive file location to be uploaded to the database. Expects a single file that is the output of DSPSR using the pulsar's ephemeris.")
    uploadargs.add_argument('--single_pulse_series',type=str,help='The single pulse series file location to be uploaded to the database. Expects a single file that is the output of DSPSR in single pulse mode (the -s option).')
    uploadargs.add_argument('--ppps',type=str,help="The Presto Prepfold PostScript file location to be uploaded to the database. Expects a single file that is the output of PRESTO's prepfold script.")
    uploadargs.add_argument('-i','--ippd',type=str,help="The Intergrated Pulse Profile (sometimes called a pulse profile) file location. Expects a single file that is the output of DSPSR's pav script.")
    uploadargs.add_argument('-w','--waterfall',type=str,help="The file location of a waterfall plot of pulse phase vs frequency. Expects a single file that is the output of DSPSR's psrplot.")
        
    dspsrargs = parser.add_argument_group('DSPSR Calculation Options', "Requires the --fits_files option. These options are all boolean flags that when used will send off DSPSR jobs to process the needed files that can be uploaded to the database. The files will be uploaded automatically when the DSPS scripts are tested more") #TODO remove when I'm confident with dspsr
    dspsrargs.add_argument('-f','--fits_files',type=str,help='The fits files location to be used in any Detection Calculation processing. Recommended to end in *.fits and surrounded by quotation marks.', default = "/group/mwaops/vcs/${obsid}/fits/*fits")
    dspsrargs.add_argument('--u_archive', action='store_true',help='Used to create an archive file using DSPSR.')
    dspsrargs.add_argument('--u_single_pulse_series', action='store_true',help='Used to create a single pulse archive using DSPSR.')
    dspsrargs.add_argument('--u_ppps', action='store_true', help="Used to create a Presto Prepfold PostScript file using PRESTO")
    dspsrargs.add_argument('--u_ippd', action='store_true', help="Used to create an Intergrated Pulse Profile (sometimes called a pulse profile) using DSPSR")
    dspsrargs.add_argument('--u_waterfall', action='store_true', help="Used to create a waterfall plot of pulse phase vs frequency using DSPSR.")
    args=parser.parse_args()

    if 'MWA_PULSAR_DB_USER' in os.environ and 'MWA_PULSAR_DB_PASS' in os.environ:
        auth = (os.environ['MWA_PULSAR_DB_USER'],os.environ['MWA_PULSAR_DB_PASS'])
    else:
        auth = ('mwapulsar','veovys9OUTY=')
        print "No MWA Pulsar Database username and password found so using the defaults."
        print 'Please add ". ~/.MWA_Pulsar_DB_profile" to your .bashrc where .MWA_Pulsar_DB_profile contains: '
        print 'export MWA_PULSAR_DB_USER="<username>"'
        print 'export MWA_PULSAR_DB_PASS="<password>"'
        print 'replacing <username> <password> with your MWA Pulsar Database username and password.'

    #defaults for incoh and calibrator type
    if args.incoh:
        coh = False
        calibrator_type = None
    else:
        coh = True
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
        #The obsid can only be obtained from the bestprof file if the 
        #fits file name start with the obsid
        if obsid is None and args.obsid:
            obsid = args.obsid
        elif obsid is None:
            print "Please use --obsid. Exiting"
            sys.exit(0)
    elif args.obsid and args.pulsar:
        num_bins = 128 #used in dspsr calculations
        obsid = args.obsid
        pulsar = args.pulsar
    elif args.obsid and (args.calibration or args.cal_dir_to_tar) and not\
         (args.archive or args.single_pulse_series or args.ppps or args.ippd\
          or args.waterfall or args.u_archive or args.u_single_pulse_series\
          or args.u_ppps or args.u_ippd  or args.u_waterfall):
        #uploaded calibrator
        obsid = args.obsid
    else:
        print "Please us either (--obsid and --pulsar) or --bestprof"
        quit()
        
    if args.pulsar or args.bestprof:
        #Checks to see if the pulsar is already on the database
        pul_list_dict = client.pulsar_list(web_address, auth)
        pul_list_str = ''
        for p in pul_list_dict:
            pul_list_str = pul_list_str + p[u'name']
        if pulsar in pul_list_str:
            print 'This pulsar is already on the database'
            
            #gets Ra and DEC from PSRCAT
            pulsar_ra_dec = fpio.get_psrcat_ra_dec(pulsar_list=[pulsar])
            pulsar_name, pul_ra, pul_dec = pulsar_ra_dec[0]
        else:
            print 'Congratulations you have detected ' + pulsar + ' for the first time with the MWA'
            #gets Ra and DEC from PSRCAT
            pulsar_ra_dec = fpio.get_psrcat_ra_dec(pulsar_list=[pulsar])
            pulsar_name, pul_ra, pul_dec = pulsar_ra_dec[0]

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
        subbands = flux_cal_and_sumbit(time_detection, time_obs, metadata, 
                            bestprof_data, args.bestprof,
                            pul_ra, pul_dec, coh, auth,
                            start = args.start, stop = args.stop, trcvr = args.trcvr)

    if args.cal_dir_to_tar:
        if not args.srclist:
            print "You must use --srclist to define the srclist file location. Exiting"
            quit()
        args.calibration = zip_calibration_files(args.cal_dir_to_tar, args.cal_id, args.srclist)


    if args.pulsar and not args.bestprof:  
        #calc sub-bands
        subbands = 1
        for b in range(len(channels)):
            if b == 0:
                continue
            if not (channels[b] - channels[b-1]) == 1:
                subbands = subbands + 1
                
        if coh:
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
                                    coherent = coh,
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
                                    coherent = coh,
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
                                    coherent = coh,
                                    startcchan = int(minfreq), stopcchan = int(maxfreq), 
                                    observation_type = int(obstype))  
            temp_dict = client.detection_get(web_address, auth, observationid = str(obsid))  

    #Upload analysis files to the database
    if args.bestprof:
        print "Uploading bestprof file to database"
        cp(str(args.bestprof) ,str(obsid) + "_" + str(pulsar) + ".bestprof")
        d_file_loc = str(obsid) + "_" + str(pulsar) + ".bestprof"
        print "Uploading Presto Prepfold PostScript file to database"
        client.detection_file_upload(web_address, auth, 
                            observationid = str(obsid),
                            pulsar = str(pulsar), 
                            subband = int(subbands),
                            coherent = coh,
                            filetype = 5,
                            filepath = str(d_file_loc))
        os.system("rm " + d_file_loc)
        
    if args.archive:
        print "Uploading archive file to database"
        client.detection_file_upload(web_address, auth, 
                                    observationid = str(obsid),
                                    pulsar = str(pulsar), 
                                    subband = int(subbands),
                                    coherent = coh,
                                    filetype = 1,
                                    filepath = str(args.archive))

    if args.single_pulse_series:
        print "Uploading single_pulse_series file to database"
        client.detection_file_upload(web_address, auth,
                                    observationid = str(obsid),
                                    pulsar = str(pulsar), 
                                    subband = int(subbands),
                                    coherent = coh,
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
                            coherent = coh,
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
                            coherent = coh,
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
                            coherent = coh,
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
            commands.append("srun -n 1 -c $ncpus dspsr -U 600 -E {0}.par -b {1} -A -cont -O {2}_{0} {3}".format(pulsar, num_bins, obsid, fits_files_loc))
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
            commands.append("rm {0}.eph".format(pulsar))
        if args.u_single_pulse_series:
            commands.append("srun -n 1 dspsr -U 600 -E {0}.par -b {1} -cont -s -K ".\
                                format(pulsar,num_bins) + "${fits}")
            commands.append('psraddstring="psradd -o {0}_{1}.ts.ar "'.format(obsid,pulsar))
            commands.append("ts=(pulse*.ar)")
            commands.append('for ((i=0;i<${#ts[@]};i++)); do psraddstring=${psraddstring}" "${ts[i]} ; done')
            commands.append("srun -n 1 -c $ncpus $psraddstring")
            commands.append("rm pulse*.ar")
        commands.append("rm {0}.par".format(pulsar))
        job_id = submit_slurm(dspsr_batch, commands,
                              batch_dir="./",
                              slurm_kwargs={"time": "6:50:00", "partition": "workq"},
                              submit=True)
    
