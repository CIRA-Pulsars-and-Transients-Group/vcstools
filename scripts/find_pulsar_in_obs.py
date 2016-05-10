#! /usr/bin/env python 

"""
Author: Nicholas Swainston
Creation Date: /03/2016

Some of the orignal code was created by Bradley Meyers

This code is used to list the sources within the beam of observations IDs or using --obs_for_source list all the observations for each source. The sources can be input serval ways: using a list of pulsar names (--pulsar), using a complete catalogue file of pulsars (--dl_PSRCAT) or RRATs (--RRAT and --dl_RRAT), using a compatable catalogue (--in_cat with the help of --names and --coordstype) or using a RA and DEC coordinate (--coords). The observation IDs can be input (--obsid) or gathered from a directory (--FITS_dir). The default is to search all observation IDs from http://mwa-metadata01.pawsey.org.au/metadata/ that have voltages and list every known pulsar from PSRCAT in each observation ID.
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
    

def deg2sex( ra, dec):
    """
    deg2sex( ra, dec)
    ra - the right ascension in degrees
    dec - the declination in degrees
    
    Convert decimal coordingates into sexagesimal strings, i.e. hh:mm:ss.ss and dd:mm:ss.ss
    """
    
    c = SkyCoord( ra, dec, frame='icrs', unit=(u.deg,u.deg))
    coords = c.to_string('hmsdms')
    coords = coords.replace('h',':').replace('d',':').replace('m',':').replace('s','')
    
    # return RA and DEC in "hh:mm:ss.ssss dd:mm:ss.ssss" form
    return coords
    
    
def grab_RRATalog():
    """
    grab_RRATalog()
    
    Creates a catalogue csv file using data from http://astro.phys.wvu.edu/rratalog/rratalog.txt
    """
    rratalog_website = 'http://astro.phys.wvu.edu/rratalog/rratalog.txt'

    print "Retrieveing base RRATalog from {0}".format( rratalog_website)
    os.system( 'rm rratalog.txt')
    os.system( 'wget {0}'.format( rratalog_website))

    print "Converting to CSV format for easier use..."
    txt_file = "rratalog.txt"
    csv_file = "rratalog.csv"


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
            for entry in columns:
                if entry not in ['', ' ', '\t']:
                    temp.append(entry.replace('--','')) 
            data.append(temp[0:14])


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
    

    print ','.join(header)
    with open(csv_file,"wb") as out_csv:
        out_csv.write(','.join(header)+'\n')
        for d in data:
            out_csv.write(','.join(d)+'\n')
    return

    
def grab_pulsaralog(jlist=None):
    """
    grab_pulsaralog(jlist=None)
    jlist - A space seperated string of pulsar names eg: J0534+2200 J0538+2817. (default: uses all pulsars)
    
    Uses PSRCAT and returns every pulsar in the catalouge in a csv file with the requested paramaters. Removes pulsars without any RA or DEC recorded
    """
    params = ['Jname', 'Raj', 'Decj', 'P0', 'P1', 'DM']
    #If more paramaters are needed add them above
    #The proper motion is not accounted for as it is assumed that the beam is not accurate enought to be necessary
    pulsars = [[]]
    if jlist != None:
        for p in params:
            #Gets the output of PSRCAT for each pparameter for each pulsar as a list
            cmd = ['psrcat', '-c', p]
            for j in jlist:
                cmd.append(j)
            output = subprocess.Popen(cmd,stdout=subprocess.PIPE).communicate()[0]
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
    else:
        for p in params:
            #Gets the output of PSRCAT for each pparameter for each pulsar as a list
            cmd = ['psrcat', '-c', p]
            output = subprocess.Popen(cmd,stdout=subprocess.PIPE).communicate()[0]
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
        os.system( 'rm -f pulsaralog.csv')
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
    calcFWHM(freq)
    freq - observation frequency in MHz
    
    Calculate the FWHM for the beam assuming ideal response/gaussian-like profile. Will eventually be depricated.
    """
    c = 299792458.0                 # speed of light (m/s)
    Dtile = 4.0                     # tile size (m) - incoherent beam
    freq = freq * 1e6               # convert from MHz to Hz
    fwhm = 1.2 * c / (Dtile * freq) # calculate FWHM using standard formula
 
    return fwhm


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
  
    
def singles_source_search(ra, dec):
    """
    singles_source_search(ra, dec)
    ra = ra of source in degrees
    dec = dec of source in degrees
    
    Used to creat a a 30 degree box around the source to make searching for obs_ids more efficient
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
        temp = getmeta(service='find', params={'mode':'VOLTAGE_START','limit':10000,'minra':0.,\
                                               'maxra':360.,'mindec':dec_bot,'maxdec':dec_top})
        for row in temp:
            OBSID.append(row[0])
    else:
        ra_low = ra - 30. - box_size #30 is the how far an obs would drift in 2 hours(used as a max)
        ra_high = ra + box_size
        if ra_low < 0.:
            ra_new = 360 + ra_low
            OBSID = []
            temp = getmeta(service='find', params={'mode':'VOLTAGE_START','limit':10000,'minra':ra_new,\
                                               'maxra':360.,'mindec':dec_bot,'maxdec':dec_top})
            for row in temp:
                OBSID.append(row[0])
            temp = getmeta(service='find', params={'mode':'VOLTAGE_START','limit':10000,'minra':0.,\
                                               'maxra':ra_high,'mindec':dec_bot,'maxdec':dec_top})
            for row in temp:
                OBSID.append(row[0])
        elif ra_high > 360:
            ra_new = ra_high - 360
            OBSID = []
            temp = getmeta(service='find', params={'mode':'VOLTAGE_START','limit':10000,'minra':ra_low,\
                                               'maxra':360.,'mindec':dec_bot,'maxdec':dec_top})
            for row in temp:
                OBSID.append(row[0])
            temp = getmeta(service='find', params={'mode':'VOLTAGE_START','limit':10000,'minra':0.,\
                                               'maxra':ra_new,'mindec':dec_bot,'maxdec':dec_top})
            for row in temp:
                OBSID.append(row[0])
        else:
            OBSID =[]
            temp = getmeta(service='find', params={'mode':'VOLTAGE_START','limit':10000,'minra':ra_low,\
                                               'maxra':ra_high,'mindec':dec_bot,'maxdec':dec_top})
            for row in temp:
                OBSID.append(row[0])
    return OBSID


def calcbox( ras, dec, time):
    """
    calcbox( ras, dec, time)
    ras - right acension at the start of the observation in degrees
    dec - declination at the start of the observation in degrees
    time - duration of the observation in seconds
    
    Calculate the end ra after the duration of the observation and makes a box from 
    the intial to final ra with a height of 20 degrees, 10 up and 10 down.
    """
    adjustRA = False
    RAloop = False 
    
    # check that the initial RA is sensible (i.e. 360=0, not negative and not > 360)    
    if ras.is_integer() and int(ras) == 360:
        ras = 0.0
    elif int(1e15 * (ras - 360.0)) >= 1:
        ras = ras - 360.0
    elif int(1e15 * ras) < 0:
        ras = ras + 360.0  


    # define the central DEC and the top and bottom corners of the LHS of the box
    # as data is drift scan, DEC does not change from start to finish
    dec_top = dec + 10.
    dec_bot = dec - 10.
    
    
    # if the botDEC is suddenly more negative than -90 (i.e. abs(botDEC)>90) then we've gone past the pole.
    # the DEC then needs to become -180+abs(botDEC) and the RA is changed by 180 degrees    
    if dec_bot < -90.0:
        dec_bot = -180.0 + math.fabs(dec_bot)
        adjustRA = True
    
    #calc end RA
    raf = ras + (time / 60.0**2) * 15.0  
    
    
    # check to see if the RA needs to be adjusted due to passing over the pole   
    # if RA>180 then new RA = old RA -180, otherwise it's oldRA +180 
    if adjustRA is True:
        if int(1e15 * (raf - 180.0)) >= 1:
            raf = raf - 180.0
        else:
            raf = raf + 180.0


    # check that final RA is sensible (i.e. 360=0, not negative and not > 360)
    if raf.is_integer() and int(raf) == 360:
        raf = 0.0
    elif int(1e15 * (raf - 360.0)) >= 1:
        raf = raf - 360.0
    elif int(1e15 * (raf)) < 0:
        eraf = raf + 360.0
      
        
    # list of start and end coordinates describing "observed rectangle"
    box = [ras, raf, dec_top, dec_bot] 
    return box
    
def beamcheck(beam_input, ras, decs):
    """
    beamcheck(beam_input, ras, decs)
    beam_input - a list of [obs ID, RA at the centre of the beam, DEC at the centre of the beam,
                 duration of the observation]
    ras - right ascension of the source in degrees
    dec - declination of the source in degrees
    
    Creates a box and a cicle at the start and end of the obs then checks if the source is inside it
    """
    found = False
    
    #extract some inital data
    source = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')
    obs, ra_beam, dec_beam, time_beam, delays, centrefreq, channels = beam_input #in degrees from metadata
    ra_intial, ra_final, dec_top, dec_bot = calcbox( ra_beam, dec_beam, time_beam)
    
    #create astropy coordinates for the start and end of the file
    centrebeam_start = SkyCoord(ra=ra_intial*u.degree, dec=dec_beam*u.degree, frame='icrs')
    centrebeam_end = SkyCoord(ra=ra_final*u.degree, dec=dec_beam*u.degree, frame='icrs')
    
    #checks if the source is within 10 degrees of the start and end of the file
    angdif_start = centrebeam_start.separation(source).degree
    angdiff_start = float(angdif_start)
    angdif_end = centrebeam_end.separation(source).degree
    angdiff_end = float(angdif_end)
    
    #loop is inaccurate for beams near the south pole
    #check the circle beam at the start and end of the observation and the rectangle connecting them
    if ra_intial > ra_final:
        if (angdiff_start < 10.) or (angdiff_end < 10.) or \
           ( ( ((ras > ra_intial) and (ras < 360.)) or ((ras > 0.) and (ras < ra_final)) ) and \
           ((dec_top > decs) and (dec_bot < decs)) ):
            found = True
        #TODO: just fyi, if you are comparing a value to two limits, i.e. want to know if x>1 and x<10
        # the equivalent in  python is actually just: if 1<x<10: *do stuff*, rather than having to use a million "and"
    else:
        if (angdiff_start < 10.) or (angdiff_end < 10.) or \
           ( ((ras > ra_intial) and (ras < ra_final)) and ((dec_top > decs) and (dec_bot < decs)) ):
            found = True
    return found
  
    
def circlebeam(beam, coord1, coord2, sourcelist, names):
    """
    circlebeam(beam, coord1, coord2, sourcelist, names)
    beam - a list of [obs ID, RA at the centre of the beam, DEC at the centre of the beam,
           duration of the observation]
    coord1 - name of the first coordinate catalogue column (normaly Raj)
    coord2 - name of the second coordiante catalogue column (normaly Decj)
    sourcelist - a table of sources readable by astropy.table.Table
    names - the name of the catalogue column of the source names (normaly Jname)
    
    Lists all of the pulsars in the beam using a simple circle of radius of 10 
    degrees and the entire PSRCAT catalouge. 
    """
    sourceinbeam = {}
    obs = beam[0]
    print "Finding pulsars within 10 degrees of the beam centre of " + str(obs)
    #will loop for time later
    temp = []
    for s in sourcelist:
        #chooses the colomn name of the sources 
        if names=='-1':
            x,y = sex2deg(s[coord1],s[coord2])
            name = str(round(x,3))+'_'+str(round(y,3))
        else:
            name=s[str(names)]
    
        ras, decs = sex2deg(s[coord1],s[coord2])  # source coordinates
        ras = float(ras)
        decs = float(decs)
        found = beamcheck(beam, ras, decs)
        if found:
            temp.append(name)
            print str(name)
    #write a simple file with the list of pulsar names that can be expanded 
    #later to include their paramaters
    if temp: #check the file isn't going to be empty
        outputfile = str(args.output)
        os.system( 'rm -f ' + outputfile + str(obs) + '_circle_beam.txt')
        with open(outputfile + str(obs) + '_circle_beam.txt',"wb") as out_list:
            out_list.write('All of the sources found within a simple 10 degree radius '\
            + 'circle from the centre of the beam for observation ID: ' + str(obs) + '\n' +\
            'RA(deg): ' + str(beam[1]) + '   DEC(deg): ' + str(beam[2]) + '   Duration(s): ' + str(beam[3]) + '\n')
            for row in temp:
                 out_list.write(row + "\n")
        print "A list of sources for the observation ID: " + str(obs) + \
                  " has been output to the text file: " + str(obs) + '_circle_beam.txt'
    return


def circlebeam_obsforsource(beam, coord1, coord2, sourcelist, names):
    """
    circlebeam(beam, coord1, coord2, sourcelist, names)
    beam - a list of [obs ID, RA at the centre of the beam, DEC at the centre of the beam,
           duration of the observation]
    coord1 - name of the first coordinate catalogue column (normaly Raj)
    coord2 - name of the second coordiante catalogue column (normaly Decj)
    sourcelist - a table of sources readable by astropy.table.Table
    names - the name of the catalogue column of the source names (normaly Jname)
    
    Lists all of the observation IDs that have a the list of sources within a simple 
    circle of radius of 10 degrees. 
    """
    sourceinbeam = {}
    for s in sourcelist:
        #chooses the colomn name of the sources 
        if names=='-1':
            x,y = sex2deg(s[coord1],s[coord2])
            name = str(round(x,3))+'_'+str(round(y,3))
        else:
            name=s[str(names)]
    
        print "Finding observation IDs that contain " + str(name) + "  within 10 degrees of the beam centre"
        ras, decs = sex2deg(s[coord1],s[coord2])  # source coordinates
        ras = float(ras)
        decs = float(decs)
        temp = []
        for b in beam:
            obs = b[0]
            found = beamcheck(b, ras, decs)
            if found:
                temp.append(obs)
                print str(obs)
            
        #write a simple file with the list of obs for each source that 
        #can be expanded later to include their paramaters
        outputfile = str(args.output)
        os.system( 'rm -f ' + outputfile + str(name) + '_circle_beam.txt')
        with open(outputfile + str(name) + '_circle_beam.txt',"wb") as out_list:
            out_list.write('All of the observation IDs that found ' + str(name) +\
            ' within a simple 10 degree radius circle from the centre of the beam.\n' +\
            'RA(deg): ' + str(ras) + '   DEC(deg): ' + str(decs) + '\n')
            for row in temp:
                 out_list.write(str(row) + "\n")
        print "A list of observation ID for the source: " + str(name) + \
              " has been output to the text file: " + str(name) + '_circle_beam.txt'                   
    return
   
   
def beam_enter_exit(min_power, powers, imax, dt, duration):
    """
    beam_enter_exit(min_power, powers, imax, dt):
    min_power - power cut off from get_bem_power(_obsforsource)
    powers - list of powers fo the duration every dt
    imax - an int that will give the max power of the file (max = powers[imax])
    dt - the time interval of how often powers are calculated
    
    Calculates when the source enters and exits the beam 
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
    #covert it into a file number and adding a file either side incase of ccode error or bright sources
    if before_check:
        before = before_list[0]
        before_time = before[1]
        before_file = int(before_time / 200)
        if before_file > 0:
            enter = before_file - 1
        else:
            enter = before_file
    else:
        enter = 0
        
    if after_check:
        after = after_list[-1]
        after_time = after[1]
        after_file = int(after_time / 200)
        
        if after_file < file_num:
            exit = after_file + 1
        else:
            exit = after_file
    else:
        exit = file_num
    
    return [enter,exit]
    
    
def get_beam_power(obsid_data,
                   sources, coord1, coord2, names,
                   dt=296,
                   centeronly=True,
                   verbose=False,
                   min_power=0.3):
    """
    obsid_data = [obsid,ra, dec, time, delays,centrefreq, channels]
    sources=[names,coord1,coord2] #astropy table coloumn names

    Calulates the power (gain at coordinate/gain at zenith) for each source and if it is above
    the min_power then it outputs it to the text file.

    """
    print "Calculating beam power"
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
        theta=np.radians(90-Alts)
        phi=np.radians(Azs)
        
        for ifreq in xrange(len(frequencies)):
            rX,rY=primary_beam.MWA_Tile_analytic(theta, phi,
                                                 freq=frequencies[ifreq], delays=delays,
                                                 zenithnorm=True,
                                                 power=True)
            PowersX[:,itime,ifreq]=rX
            PowersY[:,itime,ifreq]=rY

    #Power [#sources, #times, #frequencies]
    Powers=0.5*(PowersX+PowersY)
    
    outputfile = str(args.output)
    os.system( 'rm -f ' + outputfile + str(obsid) + '_analytic_beam.txt')
    with open(outputfile + str(obsid) + '_analytic_beam.txt',"wb") as out_list:
        out_list.write('All of the sources that the analytic beam model calculated a power of ' +\
                        str(min_power) + ' or greater for observation ID: ' + str(obsid) + '\n' +\
                       'Observation data :RA(deg): ' + str(ra) + ' DEC(deg): ' + str(dec)+' Duration(s): ' \
                       + str(time) + '\n' + \
                       'Source      Time of max power in observation    File number source entered the '\
                       + 'beam    File number source exited the beam\n')
        counter=0
        for sourc in Powers:
            counter = counter + 1
            if max(sourc) > min_power:                
                for imax in range(len(sourc)):
                    if sourc[imax] == max(sourc):
                        max_time = midtimes[imax]
                        enter, exit = beam_enter_exit(min_power, sourc, imax, dt, time)
                        out_list.write(str(sources[names][counter - 1]) + ' ' + \
                                str(int(max_time) - int(obsid)) + ' ' + str(enter) + ' ' + str(exit) + "\n") 
        print "A list of sources for the observation ID: " + str(obsid) + \
              " has been output to the text file: " + str(obsid) + '_analytic_beam.txt'             
    return 


def get_beam_power_obsforsource(obsid_data,
                   sources, coord1, coord2, names,
                   dt=296,
                   centeronly=True,
                   verbose=False,
                   min_power=0.3):
    """
    obsid_data = [obsid,ra, dec, time, delays,centrefreq, channels]
    sources=[names,coord1,coord2] #astropy table coloumn names

    Calulates the power (gain at coordinate/gain at zenith) for each source and if it is above
    the min_power then it outputs it to the text file.
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
            theta=np.radians(90-Alts)
            phi=np.radians(Azs)
            
            for ifreq in xrange(len(frequencies)):
                rX,rY=primary_beam.MWA_Tile_analytic(theta, phi,
                                                     freq=frequencies[ifreq], delays=delays,
                                                     zenithnorm=True,
                                                     power=True)
                PowersX[:,itime,ifreq]=rX
                PowersY[:,itime,ifreq]=rY
        #temp_power [#sources, #times, #frequencies]
        temp_power=0.5*(PowersX+PowersY)
        counter = 0
        for sourc in temp_power:
            counter = counter + 1
            if max(sourc) > min_power:
                print obsid
                
                for imax in range(len(sourc)):
                    if sourc[imax] == max(sourc):
                        max_time = midtimes[imax] - obsid
                        Powers.append([sources[names][counter-1], obsid, time, max_time, sourc, imax])
    
    for sourc in sources:    
        outputfile = str(args.output)
        os.system( 'rm -f ' + outputfile + str(sourc[names]) + '_analytic_beam.txt')
        with open(outputfile + str(sourc[names]) + '_analytic_beam.txt',"wb") as out_list:
            out_list.write('All of the observation IDs that the analytic beam model calculated a power of '\
                           + str(min_power) + ' or greater for the source: ' + str(sourc[names]) + '\n' +\
                           'Obs ID     Duration  Time during observation that the power was at a'+\
                           ' maximum    File number source entered    File number source exited\n')
            for p in Powers:
                if str(p[0]) == str(sourc[names]):
                    enter, exit = beam_enter_exit(min_power, p[4], p[5], dt, time)
                    out_list.write(str(p[1]) + ' ' + str(p[2]) + ' ' + str(p[3]) + ' ' + str(enter) +\
                                   ' ' + str(exit) + "\n")
           
    print "A list of observation IDs that containt: " + str(sourc[names]) + \
              " has been output to the text file: " + str(sourc[names]) + '_analytic_beam.txt'             
    return 


parser = argparse.ArgumentParser(description="""
This code is used to list the sources within the beam of observations IDs or using --obs_for_source list all the observations for each source. The sources can be input serval ways: using a list of pulsar names (--pulsar), using a complete catalogue file of pulsars (--dl_PSRCAT) or RRATs (--RRAT and --dl_RRAT), using a compatable catalogue (--in_cat with the help of --names and --coordstype) or using a RA and DEC coordinate (--coords). The observation IDs can be input (--obsid) or gathered from a directory (--FITS_dir). The default is to search all observation IDs from http://mwa-metadata01.pawsey.org.au/metadata/ that have voltages and list every known pulsar from PSRCAT in each observation ID.
""")
parser.add_argument('--obs_for_source',action='store_true',help='Instead of listing all the sources in each observation it will list all of the observations for each source.')
parser.add_argument('--output',type=str,help='Chooses a file for all the text files to be output to. To use the current file use "./". The default is /scratch2/mwaops/pulsar/incoh_census/analytic_beam_output/', default = '/scratch2/mwaops/pulsar/incoh_census/analytic_beam_output/')
parser.add_argument('-b','--beam',type=str,help='Decides the beam approximation that will be used. Options: "c" a simple circular beam approximation of radius 10 degrees and "a" an analytic beam model. Default: "a"')
#impliment an option for the accurate beam later on and maybe the old elipse approximation if I can make it accurate

#source options
sourargs = parser.add_argument_group('Source options', 'The different options to control which sources are used. Default is all known pulsars.')
sourargs.add_argument('-p','--pulsar',type=str, nargs='*',help='Searches for all known pulsars. This is the default. To search for individual pulsars list their Jnames in the format " -p J0534+2200 J0630-2834"')
sourargs.add_argument('--RRAT',action='store_true',help='Searches for all known RRATs.')
#Eventually impliment to search for FRBs and a search for all mode
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
args=parser.parse_args()
  
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
    if args.pulsar:
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
            grab_pulsaralog(jlist)
            catDIR = 'temp.csv'
    else:
        catDIR = 'pulsaralog.csv'

#defaults for the coords types
if args.coord_names:
    c1, c2 = args.coord_names.split(',')
else:
    if args.RRAT:
        c1, c2 = ['RA','DEC']
    if args.pulsar:
        c1, c2 = ['Raj', 'Decj']
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
    if args.pulsar:
        name_col = 'Jname'
    if args.coords:
        name_col = '-1'
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
        grab_pulsaralog()
        catalog = Table.read('pulsaralog.csv')
        
header = catalog.colnames
if args.pulsar != None:
    os.system('rm -f temp.csv')
    

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
    elif args.coords and len(args.coords) == 1:
        ras, decs = sex2deg(catalog[c1][0],catalog[c2][0])
        OBSID = singles_source_search(ras, decs)
    else:
        temp = getmeta(service='find', params={'mode':'VOLTAGE_START','limit':10000})
        for row in temp:
            OBSID.append(row[0])


#for source in obs
#gets all of the basic meta data for each observation ID
#prepares metadata calls
if args.all_volt: #drops the d.filetype = 11
    sql_meta = ('select a.starttime, a.stoptime-a.starttime as duration, m.ra_pointing, m.dec_pointing, r.frequencies, d.filename from mwa_setting as a '
                'inner join rf_stream as r on a.starttime = r.starttime '
                'inner join schedule_metadata as m on a.starttime = m.observation_number '
                'inner join data_files as d on a.starttime = d.observation_num '
                'where a.starttime = %s')
else:
    sql_meta = ('select a.starttime, a.stoptime-a.starttime as duration, m.ra_pointing, m.dec_pointing, r.frequencies, d.filename from mwa_setting as a '
                'inner join rf_stream as r on a.starttime = r.starttime '
                'inner join schedule_metadata as m on a.starttime = m.observation_number '
                'inner join data_files as d on a.starttime = d.observation_num '
                'where d.filetype = 11 and a.starttime = %s')

sql_delay = ('select xdelaysetting from obsc_recv_cmds where observation_number = %s')

#downloads password so only mwaops group can access the archive
password_parser = SafeConfigParser()
password_parser.read('/scratch2/mwaops/pulsar/incoh_census/beam_code/MWA_metadata.ini')

cord = []
delays = []
try:
    with closing(psycopg2.connect(database='mwa', user='mwa_ro', password=password_parser.get('MWA_admin','password') , host='mwa-metadata01', port='5432')) as conn:
        print ' '
except:#incase of file finding or permission errors use webservice
    print 'Error using admin account. Using slower webservice instead.'
    for ob in OBSID:
        print "Obtaining metadata from http://mwa-metadata01.pawsey.org.au/metadata/ for OBS ID: " + str(ob)
        beam_meta_data = getmeta(service='obs', params={'obs_id':ob})
        ra = beam_meta_data[u'metadata'][u'ra_pointing']
        dec = beam_meta_data[u'metadata'][u'dec_pointing']
        time = beam_meta_data[u'stoptime'] - beam_meta_data[u'starttime'] #gps time 
        skytemp = beam_meta_data[u'metadata'][u'sky_temp']
        delays = beam_meta_data[u'rfstreams'][u'0'][u'delays']

        channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
        minfreq = float(min(channels))
        maxfreq = float(max(channels))
        centrefreq = 1.28e6 * (minfreq + (maxfreq-minfreq)/2) #in MHz
        
        #check for raw volatge files
        filedata = beam_meta_data[u'files']
        keys = filedata.keys()
        check = False 
        for k in keys:
            if '.dat' in k:
                check = True
        if check or args.all_volt:
            if args.obs_for_source:
                if args.obsid:
                    if len(args.obsid) == 1:
                        cord = [[ob, ra, dec, time, delays,centrefreq, channels]]
                else:
                    cord.append([ob, ra, dec, time, delays,centrefreq, channels])
            else:
                cord = [ob, ra, dec, time, delays,centrefreq, channels]
                if args.beam == 'c':
                    circlebeam(cord, c1, c2, catalog, name_col)
                elif args.beam == 'a':    #center only means it isn't in picket fence mode
                    get_beam_power(cord, catalog, c1, c2, name_col, dt=100,centeronly=True, verbose=False)
                else: #TODO impliment a picket fence mode
                    get_beam_power(cord, catalog, c1, c2, name_col, dt=100,centeronly=True, verbose=False)
        else:
            print('No raw voltage files for %s' % ob) 
            
else:
    with closing(psycopg2.connect(database='mwa', user='mwa_ro', password=password_parser.get('MWA_admin','password') , host='mwa-metadata01', port='5432')) as conn:
        print 'Admin access to http://mwa-metadata01.pawsey.org.au/metadata/ obtained'
        for ob in OBSID:
            print "Obtaining metadata from http://mwa-metadata01.pawsey.org.au/metadata/ for OBS ID: " + str(ob)
            with closing(conn.cursor()) as cur:
                cur.execute(sql_meta, (ob,))
                meta_result = cur.fetchall()
                if not meta_result:
                    print('No raw voltage files for %s' % ob) 
                    continue
               
                cur.execute(sql_delay, (ob,))
                delay_result = cur.fetchall()
                if not delay_result:
                    raise Exception('could not get delay meta data for %s' % ob)
                ra = meta_result[0][2]
                dec = meta_result[0][3]
                time = meta_result[0][1] #duration
                channels = meta_result[0][4]
                delays = delay_result[0][0][0]
                minfreq = float(min(channels))
                maxfreq = float(max(channels))
                centrefreq = 1.28e6 * (minfreq + (maxfreq-minfreq) / 2) #in Hz
            
                #instead of downloading all of the obs id first, if not in obs_for_source mode, 
                #downlads one obs at a time
                if args.obs_for_source:
                    if args.obsid:
                        if len(args.obsid) == 1:
                            cord = [[ob, ra, dec, time, delays,centrefreq, channels]]
                    else:
                        cord.append([ob, ra, dec, time, delays,centrefreq, channels])
                else:
                    cord = [ob, ra, dec, time, delays,centrefreq, channels]
                    if args.beam == 'c':
                        circlebeam(cord, c1, c2, catalog, name_col)
                    elif args.beam == 'a':    #center only means it isn't in picket fence mode
                        get_beam_power(cord, catalog, c1, c2, name_col, dt=100,centeronly=True, verbose=False)
                    else: #TODO impliment a picket fence mode
                        get_beam_power(cord, catalog, c1, c2, name_col, dt=100,centeronly=True, verbose=False)

#chooses the beam type and whether to list the source in each obs or the obs for each source
#more options will be included later
if args.obs_for_source:
    if args.beam == 'c':
        circlebeam_obsforsource(cord, c1, c2, catalog, name_col)
    elif args.beam == 'a':    
        get_beam_power_obsforsource(cord, catalog, c1, c2, name_col, dt=100,centeronly=True, verbose=False)
    else:
        get_beam_power_obsforsource(cord, catalog, c1, c2, name_col, dt=100,centeronly=True, verbose=False)



print "The code is complete and all results have been output to text files"   

                 
"""
Program tests
python find_pulsar_in_obs.py -p --obsid 1099414416
python find_pulsar_in_obs.py -p J0534+2200 J0538+2817 --obsid 1099414416
python find_pulsar_in_obs.py -p J0630-2834 --obsid 1067285064 10689221844 1101491208 1102270216
python find_pulsar_in_obs.py -p J0534+2200 J0538+2817 --obsid 1099414416 --obs_for_source
python find_pulsar_in_obs.py -p J1921+2153 J1932+1056 J1935+1616 --obsid 1095506112 --obs_for_source
python find_pulsar_in_obs.py -p J0630-2834 J0742-2822 --obsid 1067285064 --obs_for_source

python find_pulsar_in_obs.py -c 05:34:31.973,+22:00:52.06 --obsid 1099414416 --obs_for_source               

Found in my third year project
J0437-4715 J0534+2200 J0630-2834 J0742-2822 J0835-4510 J0953+0755 J1731-4744 J1752-2806 J1900-2600 J1921+2153 J1932+1059 J1935+1616 J2048-1616 J2145-0750                 
"""
