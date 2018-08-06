#! /usr/bin/env python

"""
Author: Nicholas Swainston
Creation Date: /03/2016

Some of the orignal code was created by Bradley Meyers

This code is used to list the sources within the beam of observations IDs or using --obs_for_source list all the observations for each source. The sources can be input serval ways: using a list of pulsar names (--pulsar), using a complete catalogue file of pulsars (--dl_PSRCAT) or RRATs (--RRAT and --dl_RRAT), using a compatable catalogue (--in_cat with the help of --names and --coordstype) or using a RA and DEC coordinate (--coords). The observation IDs can be input (--obsid) or gathered from a directory (--FITS_dir). The default is to search all observation IDs from http://mwa-metadata01.pawsey.org.au/metadata/ that have voltages and list every known pulsar from PSRCAT in each observation ID.

Two most comon uses for this code is to search a single observation ID for pulsars like so:
> find_pulsar_in_obs.py -o 1099414416
which will output a text file called 1099414416_analytic_beam.txt
or to search all observation IDs for a single pulsar like so:
> find_pulsar_in_obs.py -p J0534+2200 --obs_for_source
which will output a text file called J0534+2200_analytic_beam.txt
"""

__author__ = 'Nicholas Swainston'
__date__ = '2016-03-21'

#TODO: (BWM) might be worth trying to only load the necessary parts of modules rather than loading an ENTIRE module. The speed up will be minimal overall, but it's also just a bit neater.
# You can also just load module/part of modules internally within functions. So if a module is only used once, just load it in the local function rather than up here (globally).
import os
import sys
import math
import argparse
import urllib
import urllib2
import json
import subprocess
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astropy.time import Time
from mwapy.pb import primary_beam
import ephem
from mwapy import ephem_utils,metadata
import psycopg2
from contextlib import closing
from ConfigParser import SafeConfigParser
import mwa_metadb_utils as meta



def sex2deg( ra, dec):
    """
    Convert sexagesimal coordinates to degrees.

    sex2deg( ra, dec)
    Args:
        ra: the right ascension in HH:MM:SS
        dec: the declination in DD:MM:SS
    """
    c = SkyCoord( ra, dec, frame='icrs', unit=(u.hourangle,u.deg))

    # return RA and DEC in degrees in degrees
    return [c.ra.deg, c.dec.deg]


def deg2sex( ra, dec):
    """
    Convert decimal coordingates into sexagesimal strings, i.e. hh:mm:ss.ss and dd:mm:ss.ss

    deg2sex( ra, dec)
    Args:
        ra: the right ascension in degrees
        dec: the declination in degrees
    """

    c = SkyCoord( ra, dec, frame='icrs', unit=(u.deg,u.deg))
    coords = c.to_string('hmsdms')
    coords = coords.replace('h',':').replace('d',':').replace('m',':').replace('s','')

    # return RA and DEC in "hh:mm:ss.ssss dd:mm:ss.ssss" form
    return coords


def grab_FRBalog():
    """
    Creates a catalogue csv file using data from http://www.astronomy.swin.edu.au/pulsar/frbcat/table.php?format=html

    grab_GCalog()
    """
    rratalog_website = 'http://www.astronomy.swin.edu.au/pulsar/frbcat/table.php?format=html'

    print "Retrieveing base FRBalog from {0}".format(rratalog_website)
    os.remove('table.php?format=html')
    os.system( 'wget {0}'.format( rratalog_website))

    print "Converting to CSV format for easier use..."
    txt_file = "table.php?format=html"
    csv_file = "frbalog.csv"


    with open(txt_file,"rb") as in_txt:
        lines = in_txt.readlines()

        header = ['NAME','Raj','Decj']
        data = []
        for l in lines[7:]:
            ltemp = l.strip('<td>').split('</td>')
            #print ltemp
            if len(ltemp) > 2:
                ratemp = ltemp[7].lstrip('<td>')
                dectemp = ltemp[8].lstrip('<td>')
                temp = [ltemp[0].lstrip('<td>')[:9],ratemp,dectemp]
                print temp
                data.append(temp)

    #loop to format ra and dec
    for i in range(len(data)):
        data[i][0] = data[i][0].replace('*','')

        if data[i][1].endswith(":"):
            data[i][1]=data[i][1]+'00'

        if len(data[i][1])==5:
            data[i][1]=data[i][1]+':00'

        if len(data[i][1])==7:
            data[i][1]=data[i][1]+'0'


        if len(data[i][2])==2 or (len(data[i][2])==3 and \
                          data[i][2].startswith('-')):
            data[i][2]=data[i][2]+':00:00'

        if len(data[i][2])==5 or len(data[i][2])==6:
            data[i][2]=data[i][2]+':00'

        if len(data[i][2])==3 and data[i][2].endswith(':'):
            data[i][2]=data[i][2]+'00:00'

        if data[i][2].startswith('-') and data[i][2].endswith(':0'):
            data[i][2]=data[i][2]+'0'
            print data[i][2]

        if len(data[i][2])==7 and data[i][2].endswith(':'):
            data[i][2]=data[i][2]+'00'

        #error fix in the database
        if '.' in data[i][2][:7]:
            data[i][2] = data[i][2][:7].replace('.',':') + data[i][2][7:]

    with open(csv_file,"wb") as out_csv:
        out_csv.write(','.join(header)+'\n')
        for d in data:
            out_csv.write(','.join(d)+'\n')

    return


def grab_GCalog():
    """
    Creates a catalogue csv file using data from http://physwww.physics.mcmaster.ca/~harris/mwgc.dat

    grab_GCalog()
    """
    rratalog_website = 'http://physwww.physics.mcmaster.ca/~harris/mwgc.dat'

    print "Retrieveing base GCalog from {0}".format( rratalog_website)
    os.remove('mwgc.dat')
    os.system( 'wget {0}'.format( rratalog_website))

    print "Converting to CSV format for easier use..."
    txt_file = "mwgc.dat"
    csv_file = "gcalog.csv"


    with open(txt_file,"rb") as in_txt:
        lines = in_txt.readlines()
        lines = lines[71:229] #TODO may have to check if this changes

        header = ['ID','RA','DEC']
        data = []
        for l in lines[1:]:
            ratemp = l[25:37].rstrip().replace(' ',':')
            dectemp = l[38:51].rstrip().replace(' ',':')
            temp = [l[1:10].rstrip(),ratemp,dectemp]
            print temp
            data.append(temp)

    #loop to format ra and dec
    for i in range(len(data)):
        data[i][0] = data[i][0].replace('*','')

        if data[i][1].endswith(":"):
            data[i][1]=data[i][1]+'00'

        if len(data[i][1])==5:
            data[i][1]=data[i][1]+':00'

        if len(data[i][1])==7:
            data[i][1]=data[i][1]+'0'


        if len(data[i][2])==2 or (len(data[i][2])==3 and \
                          data[i][2].startswith('-')):
            data[i][2]=data[i][2]+':00:00'

        if len(data[i][2])==5 or len(data[i][2])==6:
            data[i][2]=data[i][2]+':00'

        if len(data[i][2])==3 and data[i][2].endswith(':'):
            data[i][2]=data[i][2]+'00:00'

        if data[i][2].startswith('-') and data[i][2].endswith(':0'):
            data[i][2]=data[i][2]+'0'
            print data[i][2]

        if len(data[i][2])==7 and data[i][2].endswith(':'):
            data[i][2]=data[i][2]+'00'

    with open(csv_file,"wb") as out_csv:
        out_csv.write(','.join(header)+'\n')
        for d in data:
            out_csv.write(','.join(d)+'\n')

    return


def grab_RRATalog(jlist=None):
    """
    Creates a catalogue csv file using data from http://astro.phys.wvu.edu/rratalog/rratalog.txt

    grab_RRATalog()
    """
    rratalog_website = 'http://astro.phys.wvu.edu/rratalog/rratalog.txt'

    print "Retrieveing base RRATalog from {0}".format( rratalog_website)
    os.remove('rratalog.txt')
    os.system( 'wget {0}'.format( rratalog_website))

    print "Converting to CSV format for easier use..."
    txt_file = "rratalog.txt"
    if jlist==None:
        csv_file = "rratalog.csv"
    else:
        csv_file = "temp.csv"


    with open(txt_file,"rb") as in_txt:
        lines = in_txt.readlines()

        header = lines[0].strip().replace(" ", '\t').split('\t')
        htemp = []
        for h in header:
            if h != '':
                htemp.append(h)
        htemp[13] = 'PulseWidth'
        header = htemp[0:14]
        data = []
        for l in lines[1:]:
            columns = l.strip().replace(" ", '\t').split('\t')
            temp = []
            if jlist == None or (columns[0] in jlist):
                for entry in columns:
                    if entry not in ['', ' ', '\t']:
                        temp.append(entry.replace('--',''))
                data.append(temp[0:14])

    #loop to format ra and dec
    for i in range(len(data)):
        data[i][0] = data[i][0].replace('*','')

        if data[i][4].endswith(":"):
            data[i][4]=data[i][4]+'00'

        if len(data[i][4])==5:
            data[i][4]=data[i][4]+':00'

        if len(data[i][4])==7:
            data[i][4]=data[i][4]+'0'


        if len(data[i][5])==2 or (len(data[i][5])==3 and \
                          data[i][5].startswith('-')):
            data[i][5]=data[i][5]+':00:00'

        if len(data[i][5])==5 or len(data[i][5])==6:
            data[i][5]=data[i][5]+':00'

        if len(data[i][5])==3 and data[i][5].endswith(':'):
            data[i][5]=data[i][5]+'00:00'

        if data[i][5].startswith('-') and data[i][5].endswith(':0'):
            data[i][5]=data[i][5]+'0'
            print data[i][5]

        if len(data[i][5])==7 and data[i][5].endswith(':'):
            data[i][5]=data[i][5]+'00'

    with open(csv_file,"wb") as out_csv:
        out_csv.write(','.join(header)+'\n')
        for d in data:
            out_csv.write(','.join(d)+'\n')
    return


def grab_pulsaralog(jlist=None):
    """
    Uses PSRCAT and returns every pulsar in the catalouge in a csv file with the requested
    paramaters. Removes pulsars without any RA or DEC recorded

    grab_pulsaralog(jlist=None)
    Args:
        jlist: A space seperated string of pulsar names eg: J0534+2200 J0538+2817.
               (default: uses all pulsars)
    """
    params = ['Jname', 'Raj', 'Decj', 'P0', 'P1', 'DM']
    #If more paramaters are needed add them above
    #The proper motion is not accounted for as it is assumed that the beam is not accurate
    #enought to be necessary
    pulsars = [[]]
    for p in params:
        #Gets the output of PSRCAT for each pparameter for each pulsar as a list
        cmd = ['psrcat', '-c', p]
        if jlist != None:
            for j in jlist:
                cmd.append(j)
        output = subprocess.Popen(cmd,stdout=subprocess.PIPE).communicate()[0]
        if output.startswith("WARNING: PSR"):
            print "Pulsar not on psrcat. Please use the --RRAT option if it's an RRAT or -c to use a position"
            quit()
        temp = []
        lines = output.split('\n')
        for l in lines[4:-1]:
            columns = l.split()
            if len(columns) > 1:
                temp.append([columns[1]])
        if p == params[0]:
            pulsars=temp
        else:
            pulsars = [pulsars[x] + temp[x] for x in range(len(pulsars))]


    i = 0
    while i < len(pulsars):
        if '*' in pulsars[i][1]:
            pulsars.remove(pulsars[i][:])
        if '*' in pulsars[i][2]:
            pulsars.remove(pulsars[i][:])
        else:
            pulsars[i][3] = pulsars[i][3].replace('*','')
            pulsars[i][4] = pulsars[i][4].replace('*','')
            pulsars[i][5] = pulsars[i][5].replace('*','')

            if pulsars[i][1].endswith(":"):
                pulsars[i][1]=pulsars[i][1]+'00'

            if len(pulsars[i][1])==5:
                pulsars[i][1]=pulsars[i][1]+':00'

            if len(pulsars[i][1])==7:
                pulsars[i][1]=pulsars[i][1]+'0'


            if len(pulsars[i][2])==2 or (len(pulsars[i][2])==3 and \
                              pulsars[i][2].startswith('-')):
                pulsars[i][2]=pulsars[i][2]+':00:00'

            if len(pulsars[i][2])==5 or len(pulsars[i][2])==6:
                pulsars[i][2]=pulsars[i][2]+':00'

            if len(pulsars[i][2])==3 and pulsars[i][2].endswith(':'):
                pulsars[i][2]=pulsars[i][2]+'00:00'

            if pulsars[i][2].startswith('-') and pulsars[i][2].endswith(':0'):
                pulsars[i][2]=pulsars[i][2]+'0'
                print pulsars[i][5]

            if len(pulsars[i][2])==7 and pulsars[i][2].endswith(':'):
                pulsars[i][2]=pulsars[i][2]+'00'

            i = i + 1

    if jlist != None:
        with open('temp.csv',"wb") as out_csv:
            out_csv.write(','.join(params)+'\n')
            for r in pulsars:
                for c in r[:-1]:
                    out_csv.write(c + ",")
                out_csv.write(r[-1])
                out_csv.write("\n")
    else:
        if os.path.exists('pulsaralog.csv'):
            os.remove('pulsaralog.csv')
        with open('pulsaralog.csv',"wb") as out_csv:
            out_csv.write(','.join(params)+'\n')
            for r in pulsars:
                for c in r[:-1]:
                    out_csv.write(c + ",")
                out_csv.write(r[-1])
                out_csv.write("\n")
    return


def calcFWHM(freq):
    """
    Calculate the FWHM for the beam assuming ideal response/gaussian-like profile. Will eventually be depricated.

    calcFWHM(freq)
    Args:
        freq: observation frequency in MHz
    """
    c = 299792458.0                 # speed of light (m/s)
    Dtile = 4.0                     # tile size (m) - incoherent beam
    freq = freq * 1e6               # convert from MHz to Hz
    fwhm = 1.2 * c / (Dtile * freq) # calculate FWHM using standard formula

    return fwhm


def singles_source_search(ra, dec):
    """
    Used to creates a 30 degree box around the source to make searching for obs_ids more efficient

    singles_source_search(ra, dec)
    Args:
        ra: ra of source in degrees
        dec: dec of source in degrees
    """
    ra = float(ra)
    dec = float(dec)
    box_size = 30.
    m_o_p = False # moved over (north or south) pole

    dec_top = dec + box_size
    if dec_top > 90.:
        dec_top = 90.
        m_o_p = True

    dec_bot = dec - box_size
    if dec_top < -90.:
        dec_top = -90.
        m_o_p = True

    if m_o_p:
        OBSID = []
        temp = meta.getmeta(service='find', params={'mode':'VOLTAGE_START','limit':10000,'minra':0.,\
                                               'maxra':360.,'mindec':dec_bot,'maxdec':dec_top})
        for row in temp:
            OBSID.append(row[0])
    else:
        ra_low = ra - 30. - box_size #30 is the how far an obs would drift in 2 hours(used as a max)
        ra_high = ra + box_size
        if ra_low < 0.:
            ra_new = 360 + ra_low
            OBSID = []
            temp = meta.getmeta(service='find', params={'mode':'VOLTAGE_START','limit':10000,'minra':ra_new,\
                                               'maxra':360.,'mindec':dec_bot,'maxdec':dec_top})
            for row in temp:
                OBSID.append(row[0])
            temp = meta.getmeta(service='find', params={'mode':'VOLTAGE_START','limit':10000,'minra':0.,\
                                               'maxra':ra_high,'mindec':dec_bot,'maxdec':dec_top})
            for row in temp:
                OBSID.append(row[0])
        elif ra_high > 360:
            ra_new = ra_high - 360
            OBSID = []
            temp = meta.getmeta(service='find', params={'mode':'VOLTAGE_START','limit':10000,'minra':ra_low,\
                                               'maxra':360.,'mindec':dec_bot,'maxdec':dec_top})
            for row in temp:
                OBSID.append(row[0])
            temp = meta.getmeta(service='find', params={'mode':'VOLTAGE_START','limit':10000,'minra':0.,\
                                               'maxra':ra_new,'mindec':dec_bot,'maxdec':dec_top})
            for row in temp:
                OBSID.append(row[0])
        else:
            OBSID =[]
            temp = meta.getmeta(service='find', params={'mode':'VOLTAGE_START','limit':10000,'minra':ra_low,\
                                               'maxra':ra_high,'mindec':dec_bot,'maxdec':dec_top})
            for row in temp:
                OBSID.append(row[0])
    return OBSID


def beam_enter_exit(min_power, powers, imax, dt, duration):
    """
    Calculates when the source enters and exits the beam

    beam_enter_exit(min_power, powers, imax, dt):
        min_power: power cut off from get_bem_power(_obsforsource)
        powers: list of powers fo the duration every dt
        imax: an int that will give the max power of the file (max = powers[imax])
        dt: the time interval of how often powers are calculated
    """
    #calc num of files inn obs
    file_num = duration / 200
    if (float(duration) % 200.) != 0.0:
        file_num = file_num + 1

    #Creates cut offs around min_power
    low_limit = min_power - 0.02
    high_limit = min_power + 0.02

    #creates lists of times before and after the max
    before_max = []
    for i in range(0,imax):
        before_max.append([powers[i],i*dt])

    after_max = []
    for i in range(imax,len(powers)):
        after_max.append([powers[i],i*dt])

    #Shortens the list to powers near min_power
    before_list = []
    before_check = False
    for b in before_max:
        if low_limit < b[0] < high_limit:
            before_list.append(b)
            before_check = True

    after_list = []
    after_check = False
    for a in after_max:
        if low_limit < a[0] < high_limit:
            after_list.append(a)
            after_check = True


    #assumes min power is at start of before_list and end of after_list. Then extract the time and
    #covert it into a fraction of obs and adding a file either side incase of ccode error or bright sources
    if before_check:
        before = before_list[0]
        before_time = float(before[1])
        if before_time > 199.:
            enter = float((before_time - 200.)/float(duration))
        else:
            enter = 0.0
    else:
        enter = 0.0

    if after_check:
        after = after_list[-1]
        after_time = float(after[1])

        if after_time < (float(duration) - 199.):
            exit = float((after_time + 200.)/float(duration))
        else:
            exit = 1.0
    else:
        exit = 1.0

    return [enter,exit]


def cal_on_database_check(obsid):
    from mwa_pulsar_client import client
    web_address = 'mwa-pawsey-volt01.pawsey.ivec.org'
    auth = ('mwapulsar','veovys9OUTY=')
    detection_list = client.detection_list(web_address, auth)
    
    cal_used = False
    cal_avail = False
    for d in detection_list:
        if int(d[u'observationid']) == int(obsid):
            if d[u'calibrator'] is not None:
                cal_used = True
                #TODO add a check if there is a cal file option
    
    #No cal
    check_result = 'N'
    if cal_avail:
        #Cal available
        check_result = 'A'
    elif cal_used:
        #Cal used
        check_result = 'U'
    
    return check_result

def get_beam_power(obsid_data,
                   sources, 
                   coord1 = 'Raj',
                   coord2 = 'Decj',
                   names = 'Jname',
                   dt=296,
                   centeronly=True,
                   verbose=False,
                   min_power=0.3,
                   option = 'a',
                   output = "./",
                   cal_check = False):
    """
    Calulates the power (gain at coordinate/gain at zenith) for each source and if it is above
    the min_power then it outputs it to the text file.

    get_beam_power(obsid_data,
                   sources, coord1, coord2, names,
                   dt=296,
                   centeronly=True,
                   verbose=False,
                   min_power=0.3)
    Args:
        obsid_data: [obsid,ra, dec, time, delays,centrefreq, channels] obsid data
        sources: astropy table catalogue of sources
        coord1: the coordinate type for sources (usualy RA)
        coord2: the coordinate type for sources (usualy Dec)
        names: the name given to the sources (usualy JName)
        dt: time step in seconds for power calculations (default 296)
        centeronly: only calculates for the centre frequency (default True)
        verbose: prints extra data to (default False)
        min_power: minimum power that the code will print a text file for (default (0.3)
    """
    print "Calculating beam power"
    print obsid_data
    obsid,ra, dec, time, delays,centrefreq, channels = obsid_data
    starttimes=np.arange(0,time,dt)
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

    RAs, Decs = sex2deg(sources[coord1],sources[coord2])
    if len(RAs)==0:
        sys.stderr.write('Must supply >=1 source positions\n')
        return None
    if not len(RAs)==len(Decs):
        sys.stderr.write('Must supply equal numbers of RAs and Decs\n')
        return None
    for itime in xrange(Ntimes):
        obstime = Time(midtimes[itime],format='gps',scale='utc')
        observer.date = obstime.datetime.strftime('%Y/%m/%d %H:%M:%S')
        LST_hours = observer.sidereal_time() * ephem_utils.HRS_IN_RADIAN

        HAs = -RAs + LST_hours * 15
        Azs, Alts = ephem_utils.eq2horz(HAs, Decs, mwa.lat)
        # go from altitude to zenith angle

        #Decide on beam model
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
        elif option == 'd':
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
        elif option == 'e':
            beam_string = "full_EE"
            theta=np.array([np.radians(90-Alts)])
            phi=np.array([np.radians(Azs)])

            for ifreq in xrange(len(frequencies)):
                rX,rY=primary_beam.MWA_Tile_full_EE(theta, phi,
                                                     freq=frequencies[ifreq], delays=delays,
                                                     zenithnorm=True,
                                                     power=True)
                PowersX[:,itime,ifreq]=rX
                PowersY[:,itime,ifreq]=rY
            print '{0:.2f}'.format(100.*float(itime)/float(Ntimes))+"% complete for obsid: "+str(obsid)
    #Power [#sources, #times, #frequencies]
    Powers=0.5*(PowersX+PowersY)

    outputfile = output
    primary_beam_output_file = outputfile + str(obsid) + '_' + beam_string + '_beam.txt'
    if os.path.exists(primary_beam_output_file):
        os.remove(primary_beam_output_file)
    with open(primary_beam_output_file,"wb") as out_list:
        out_list.write('All of the sources that the ' + beam_string + ' beam model calculated a power of '\
                       + str(min_power) + ' or greater for observation ID: ' + str(obsid) + '\n' +\
                       'Observation data :RA(deg): ' + str(ra) + ' DEC(deg): ' + str(dec)+' Duration(s): ' \
                       + str(time) + '\n')
        if cal_check:
            #checks the MWA Pulsar Database to see if the obsid has been 
            #used or has been calibrated
            print "Checking the MWA Pulsar Databse for the obsid: {0}".format(obsid)
            cal_check_result = cal_on_database_check(obsid)
            out_list.write("Calibrator Availability: {0}\n".format(cal_check_result))
        out_list.write('Source  Time of max power in observation    File number source entered the '\
                       + 'beam    File number source exited the beam       Max Power \n')
        counter=0
        for sourc in Powers:
            counter = counter + 1
            if max(sourc) > min_power:
                for imax in range(len(sourc)):
                    if sourc[imax] == max(sourc):
                        max_time = midtimes[imax]
                        enter, exit = beam_enter_exit(min_power, sourc, imax, dt, time)
                        out_list.write(str(sources[names][counter - 1])  + ' ' + \
                                str(int(max_time) - int(obsid)) + ' ' + str(enter) + ' ' + str(exit) +\
                                ' ' + str(max(sourc)[0]) + "\n")
        print "A list of sources for the observation ID: " + str(obsid) + \
              " has been output to the text file: " + str(obsid) + '_' + beam_string + '_beam.txt'
    return enter, exit


def get_beam_power_obsforsource(obsid_data,
                   sources,
                   coord1 = 'Raj',
                   coord2 = 'Decj',
                   names = 'Jname',
                   dt=296,
                   centeronly=True,
                   verbose=False,
                   min_power=0.3,
                   option = 'a',
                   output = "./",
                   cal_check = False):
    """
    Calulates the power (gain at coordinate/gain at zenith) for each source and if it is above
    the min_power then it outputs it to the text file.

    get_beam_power(obsid_data,
                   sources, coord1, coord2, names,
                   dt=296,
                   centeronly=True,
                   verbose=False,
                   min_power=0.3)
    Args:
        obsid_data: [obsid,ra, dec, time, delays,centrefreq, channels] obsid data
        sources: astropy table catalogue of sources
        coord1: the coordinate type for sources (usualy RA)
        coord2: the coordinate type for sources (usualy Dec)
        names: the name given to the sources (usualy JName)
        dt: time step in seconds for power calculations (default 296)
        centeronly: only calculates for the centre frequency (default True)
        verbose: prints extra data to (default False)
        min_power: minimum power that the code will print a text file for (default (0.3)
    """
    print "Calculating beam power"
    Powers= []
    for ob in obsid_data:
        obsid,ra, dec, time, delays,centrefreq, channels = ob

        starttimes=np.arange(0,time,dt)
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

        RAs, Decs = sex2deg(sources[coord1],sources[coord2])
        if len(RAs)==0:
            sys.stderr.write('Must supply >=1 source positions\n')
            return None
        if not len(RAs)==len(Decs):
            sys.stderr.write('Must supply equal numbers of RAs and Decs\n')
            return None
        for itime in xrange(Ntimes):
            obstime = Time(midtimes[itime],format='gps',scale='utc')
            observer.date = obstime.datetime.strftime('%Y/%m/%d %H:%M:%S')
            LST_hours = observer.sidereal_time() * ephem_utils.HRS_IN_RADIAN

            HAs = -RAs + LST_hours * 15
            Azs, Alts = ephem_utils.eq2horz(HAs, Decs, mwa.lat)
            # go from altitude to zenith angle

            #Decide on beam model
            if option == 'a':
                beam_string = "_analytic"
                theta=np.radians(90-Alts)
                phi=np.radians(Azs)

                for ifreq in xrange(len(frequencies)):
                    rX,rY=primary_beam.MWA_Tile_analytic(theta, phi,
                                                         freq=frequencies[ifreq], delays=delays,
                                                         zenithnorm=True,
                                                         power=True)
                    PowersX[:,itime,ifreq]=rX
                    PowersY[:,itime,ifreq]=rY
            elif option == 'd':
                beam_string = "_advanced"
                theta=np.array([np.radians(90-Alts)])
                phi=np.array([np.radians(Azs)])

                for ifreq in xrange(len(frequencies)):
                    rX,rY=primary_beam.MWA_Tile_advanced(theta, phi,
                                                         freq=frequencies[ifreq], delays=delays,
                                                         zenithnorm=True,
                                                         power=True)
                    PowersX[:,itime,ifreq]=rX
                    PowersY[:,itime,ifreq]=rY
            elif option == 'e':
                beam_string = "_full_EE"
                theta=np.array([np.radians(90-Alts)])
                phi=np.array([np.radians(Azs)])

                for ifreq in xrange(len(frequencies)):
                    rX,rY=primary_beam.MWA_Tile_full_EE(theta, phi,
                                                         freq=frequencies[ifreq], delays=delays,
                                                         zenithnorm=True,
                                                         power=True)
                    PowersX[:,itime,ifreq]=rX
                    PowersY[:,itime,ifreq]=rY
                print '{0:.2f}'.format(100.*float(itime)/float(Ntimes))+"% complete for obsid: "\
                                        +str(obsid)
            #temp_power [#sources, #times, #frequencies]
        temp_power=0.5*(PowersX+PowersY)
        counter = 0
        for sourc in temp_power:
            counter = counter + 1
            if max(sourc) > min_power:
                print obsid

                for imax in range(len(sourc)):
                    if sourc[imax] == max(sourc):
                        max_time = float(midtimes[imax]) - float(obsid)
                        Powers.append([sources[names][counter-1], obsid, time, max_time, sourc, imax, max(sourc)])

    for sourc in sources:
        outputfile = output
        primary_beam_output_file = outputfile + str(sourc[names]) + beam_string + '_beam.txt'
        if os.path.exists(primary_beam_output_file):
            os.remove(primary_beam_output_file)
        with open(outputfile + str(sourc[names]) + beam_string + '_beam.txt',"wb") as out_list:
            out_list.write('All of the observation IDs that the ' + beam_string + ' beam model calculated a power of '\
                           + str(min_power) + ' or greater for the source: ' + str(sourc[names]) + '\n' +\
                           'Obs ID   Duration  Time during observation that the power was at a'+\
                           ' maximum    File number source entered    File number source exited  Max power')
            if cal_check:
                out_list.write("   Cal\n")
            else:
                out_list.write("\n")
                
            for p in Powers:
                if str(p[0]) == str(sourc[names]):
                    enter, exit = beam_enter_exit(min_power, p[4], p[5], dt, time)
                    out_list.write(str(p[1]) + ' ' + str(p[2]) + ' ' + str(p[3]) + ' ' + str(enter) +\
                                   ' ' + str(exit) + ' ' + str(p[6][0]))
                    if cal_check:
                        #checks the MWA Pulsar Database to see if the obsid has been 
                        #used or has been calibrated
                        print "Checking the MWA Pulsar Databse for the obsid: {0}".format(p[1])
                        cal_check_result = cal_on_database_check(p[1])
                        out_list.write(" {0}\n".format(cal_check_result))
                    else:
                        out_list.write("\n")


    print "A list of observation IDs that containt: " + str(sourc[names]) + \
              " has been output to the text file: " + str(sourc[names]) + beam_string + '_beam.txt'
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    This code is used to list the sources within the beam of observations IDs or using --obs_for_source list all the observations for each source. The sources can be input serval ways: using a list of pulsar names (--pulsar), using a complete catalogue file of pulsars (--dl_PSRCAT) or RRATs (--RRAT and --dl_RRAT), using a compatable catalogue (--in_cat with the help of --names and --coordstype) or using a RA and DEC coordinate (--coords). The observation IDs can be input (--obsid) or gathered from a directory (--FITS_dir). The default is to search all observation IDs from http://mwa-metadata01.pawsey.org.au/metadata/ that have voltages and list every known pulsar from PSRCAT in each observation ID.
    """)
    parser.add_argument('--obs_for_source',action='store_true',help='Instead of listing all the sources in each observation it will list all of the observations for each source. For increased efficiency it will only search OBSIDs within the primary beam.')
    parser.add_argument('--output',type=str,help='Chooses a file for all the text files to be output to. The default is your current directory', default = './')
    parser.add_argument('-b','--beam',type=str,help='Decides the beam approximation that will be used. Options: "a" the analytic beam model (2012 model, fast and reasonably accurate), "d" the advanced beam model (2014 model, fast and slighty more accurate) or "e" the full EE model (2016 model, slow but accurate). " Default: "a"',default = "a")
    parser.add_argument('-m','--min_power',type=float,help='The minimum fraction of the zenith normalised power that a source needs to have to be recorded. Default 0.3', default=0.3)
    parser.add_argument("-V", "--version", action="store_true", help="Print version and quit")
    #impliment an option for the accurate beam later on and maybe the old elipse approximation if I can make it accurate

    #source options
    sourargs = parser.add_argument_group('Source options', 'The different options to control which sources are used. Default is all known pulsars.')
    sourargs.add_argument('-p','--pulsar',type=str, nargs='*',help='Searches for all known pulsars. This is the default. To search for individual pulsars list their Jnames in the format " -p J0534+2200 J0630-2834"')
    sourargs.add_argument('--RRAT',action='store_true',help='Searches for all known RRATs.')
    sourargs.add_argument('--GC',action='store_true',help='Searches for all known Globular Clusters.')
    sourargs.add_argument('--FRB',action='store_true',help='Searches for all known FRBs.')
    #TODO Eventually impliment to search for FRBs and a search for all mode
    sourargs.add_argument('--dl_RRAT',action='store_true',help='Download the RRATalog from http://astro.phys.wvu.edu/rratalog/ and uses this as the source catalogue.')
    sourargs.add_argument('--dl_PSRCAT',action='store_true',help='Download the Puslar alog from http://www.atnf.csiro.au/research/pulsar/psrcat/ and uses this as the source catalogue.')
    sourargs.add_argument('--in_cat',type=str,help='Location of source catalogue, must be readable by astropy.table.Table (i.e. csv, txt, votable, fits) . Default: for pulsars pulsaralog.csv from grab_pulsaralog.py and for RRATs rratalog.csv from grab_RRATalog.py')
    sourargs.add_argument('--source_names',type=str,help='String containing the column name for the source names in the input catalogue (--in_cat). If there is no such column, use: --names=-1 and the output text file will be labelled using the coordinates in degrees: <longitudinal>_<latitudinal>.txt. Default: "Jname".')
    sourargs.add_argument('--coord_names',type=str,help='String containing the two column labels of the source coordinates for the input catalouge (--in_cat). i.e.: "x,y" or "long,lat". If not provided, assumes that the coordinates are "Raj,Decj". Must be enterered as: "coord1,coord2".')
    sourargs.add_argument('-c','--coords',type=str,help='String containing coordinates in "RA,DEC". This will list the OBS IDs that contain this coordinate. Must be enterered as: "hh:mm:ss.ss,+dd:mm:ss.ss".')
    #finish above later and make it more robust to incclude input as sex or deg and perhaps other coordinte systmes

    #observation options
    obargs = parser.add_argument_group('Observation ID options', 'The different options to control which observation IDs are used. Default is all observation IDs with voltages.')
    obargs.add_argument('--FITS_dir',type=str,help='Location of FITS files on system. If not chosen will search the database for metadata.')
    obargs.add_argument('-o','--obsid',type=str,nargs='*',help='Input several OBS IDs in the format " -o 1099414416 1095506112". If this option is not input all OBS IDs that have voltages will be used')
    obargs.add_argument('--all_volt',action='store_true',help='Includes observation IDs even if there are no raw voltages in the archive. Some incoherent observation ID files may be archived even though there are raw voltage files. The default is to only include files with raw voltage files.')
    obargs.add_argument('--cal_check',action='store_true',help='Check the MWA Pulsar Database to check if the obsid has every succesfully detected a pulsar and if it has a calibration solution.')
    args=parser.parse_args()


    if args.version:
        try:
            import version
            print(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            print("Couldn't import version.py - have you installed vcstools?")
            print("ImportError: {0}".format(ie))
            sys.exit(0)

    #Parser default control
    if args.dl_RRAT:
        grab_RRATalog()

    if args.dl_PSRCAT:
        grab_pulsaralog()

    #defaults for the catalouge directory
    if args.in_cat:
        catDIR = args.in_cat
    else:
        if args.RRAT:
            catDIR = 'rratalog.csv'
        elif args.GC:
            catDIR = 'gcalog.csv'
        elif args.FRB:
            catDIR = 'frbalog.ccsv'
        elif args.pulsar:
            if args.pulsar == None:
                catDIR = 'pulsaralog.csv'
            if (len(args.pulsar) ==1) and (not args.obs_for_source):
                answer = raw_input('You are only searching for one pulsar so it is recommened that you use' +\
                                  ' --obs_for_source. Would you like to use --obs_for_source. (Y/n)')
                if (answer == 'Y') or (answer == 'y') or (answer == 'yes') or (answer == ''):
                    args.obs_for_source = True
                    print "Using option --obs_for_source"
                else:
                    print "Not using option --obs_for_source"
            if args.pulsar != None:
                #converts the list of pulsars into a string so they can be used as an agrument
                jlist = args.pulsar
                if args.RRAT:
                    grab_RRATalog(jlist)
                else:
                    grab_pulsaralog(jlist)
                catDIR = 'temp.csv'
        else:
            catDIR = 'pulsaralog.csv'

    #defaults for the coords types
    if args.coord_names:
        c1, c2 = args.coord_names.split(',')
    else:
        if args.pulsar:
            c1, c2 = ['Raj', 'Decj']
        if args.RRAT or args.GC:
            c1, c2 = ['RA','DEC']
        else:
            c1, c2 = ['Raj', 'Decj']

    #defaults for the fits dirs
    if args.FITS_dir:
        fitsDIR = args.FITS_dir
    else:
        fitsDIR = '/data_01/pulsar/fitsfiles/'

    #sets the column name for the sources in the catalouge defending on different defaults
    if args.source_names and args.source_names=='-1':
        name_col='-1'
    elif args.source_names:
        name_col = args.source_names
    else:
        if args.RRAT:
            name_col = 'Name'
        elif args.pulsar:
            name_col = 'Jname'
        elif args.coords:
            name_col = '-1'
        elif args.GC:
            name_col = 'ID'
        elif args.FRB:
            name_col = 'NAME'
        else:
            name_col = 'Jname'



    #main code
    #get cataloge
    if args.coords:
        #creates a table for a single coordinate
        racor, deccor = args.coords.split(',')
        x,y = sex2deg(racor,deccor)
        name = str(round(x,3))+'_'+str(round(y,3))
        catalog = Table([[name],[racor], [deccor]], names=(name_col,c1, c2))
    else:
        print "Creating catalogue from file", catDIR, "..."
        try:
            catalog = Table.read( catDIR)
        except IOError as e:
            print "No file {0} found. Using grab_pulsars.py to creat directory.".format( e.strerror)
            if args.RRAT:
                grab_RRATalog()
                catalog = Table.read('rratalog.csv')
            elif args.GC:
                grab_GCalog()
                catalog = Table.read('gcalog.csv')
            elif args.FRB:
                grab_FRBalog()
                catalog = Table.read('frbalog.csv')
            else:
                grab_pulsaralog()
                catalog = Table.read('pulsaralog.csv')

    header = catalog.colnames


    #get obs IDs
    if args.obsid:
        OBSID = args.obsid
    elif args.FITS_dir:
        print "Creating list of observation IDs in given FITS directory"
        obsIDs = os.walk( fitsDIR).next()[1]
        if fitsDIR == '/data_01/pulsar/fitsfiles/':
            obsIDs.remove('1133_drift')
    else:
        #queries the database for a list of all the OBS IDs with voltages
        #TODO for a one -p search limit the obsid list
        OBSID = []
        print "Obtaining list of observation IDs that have recorded voltages from http://mwa-metadata01.pawsey.org.au/metadata/"
        if args.pulsar and (len(args.pulsar) ==1): #if there is a single pulsar simply search around it
            ras, decs = sex2deg(catalog[c1][0],catalog[c2][0])
            OBSID = singles_source_search(ras, decs)
        elif args.coords:
            ras, decs = sex2deg(catalog[c1][0],catalog[c2][0])
            OBSID = singles_source_search(ras, decs)
        else:
            temp = meta.getmeta(service='find', params={'mode':'VOLTAGE_START','limit':10000})
            for row in temp:
                OBSID.append(row[0])


    #for source in obs
    #gets all of the basic meta data for each observation ID
    #prepares metadata calls
    cord = []
    for ob in OBSID:
        beam_meta_data, full_meta = meta.get_common_obs_metadata(ob, return_all = True)
        #obsid,ra_obs,dec_obs,time_obs,delays,centrefreq,channels

        #check for raw volatge files
        filedata = full_meta[u'files']
        keys = filedata.keys()
        check = False
        for k in keys:
            if '.dat' in k:
                check = True
        if check or args.all_volt:
            if args.obs_for_source:
                if args.obsid and (len(args.obsid) == 1):
                    cord = beam_meta_data
                else:
                    cord.append(beam_meta_data)
            else:
                cord = beam_meta_data
                
                if args.beam == 'e':
                    dt = 300
                else:
                    dt = 100
                get_beam_power(cord, catalog, coord1 = c1, coord2 = c2, names = name_col,
                               dt=300, centeronly=True, verbose=False, min_power=args.min_power,
                               option = args.beam, output = args.output,
                               cal_check = args.cal_check)
        else:
            print('No raw voltage files for %s' % ob)

    #chooses the beam type and whether to list the source in each obs or the obs for each source
    #more options will be included later
    if args.obs_for_source:
        if args.beam:
            get_beam_power_obsforsource(cord, catalog, coord1 = c1, coord2 = c2, 
                                        names = name_col, dt=300, centeronly=True,
                                        verbose=False, min_power=args.min_power,
                                        option=args.beam, output = args.output,
                                        cal_check = args.cal_check)
        else:
            get_beam_power_obsforsource(cord, catalog, coord1 = c1, coord2 = c2, 
                                        names = name_col, dt=100,centeronly=True,
                                        verbose=False, min_power=args.min_power,
                                        output = args.output, cal_check = args.cal_check)



    print "The code is complete and all results have been output to text files"

    #remove csv file

    if args.pulsar != None:
        if os.path.exists('temp.csv'):
            os.remove('temp.csv')
    else:
        if os.path.exists(catDIR):
            os.remove(catDIR)
