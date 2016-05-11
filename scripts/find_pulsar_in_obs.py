#! /usr/bin/env python 

"""
Author: Nicholas Swainston
Creation Date: /03/2016

Some of the orignal code was created by Bradley Meyers

This code is used to list the sources within the beam of observations IDs or using --obs_for_source list all the observations for each source. The sources can be input serval ways: using a list of pulsar names (--pulsar), using a complete catalogue file of pulsars (--dl_PSRCAT) or RRATs (--RRAT and --dl_RRAT), using a compatable catalogue (--in_cat with the help of --names and --coordstype) or using a RA and DEC coordinate (--coords). The observation IDs can be input (--obsid) or gathered from a directory (--FITS_dir). The default is to search all observation IDs from http://mwa-metadata01.pawsey.org.au/metadata/ that have voltages and list every known pulsar from PSRCAT in each observation ID.

The code becomes inaccurate when observatons get close to the south pole
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



#adding grab_* scripts to main code
  



def sex2deg( ra, dec):
    """
    Convert sexagesimal coordinates to degrees.
    """ #TODO: (BWM) what format does the function expect? just brief example: i.e. "hh:mm:ss", "dd:mm:ss" 
    c = SkyCoord( ra, dec, frame='icrs', unit=(u.hourangle,u.deg))
    
    # return RA and DEC in degrees in degrees
    return [c.ra.deg, c.dec.deg]
    

def deg2sex( ra, dec):
    """
    Convert decimal coordingates into sexagesimal strings, i.e. hh:mm:ss.ss and dd:mm:ss.ss
    """
    
    c = SkyCoord( ra, dec, frame='icrs', unit=(u.deg,u.deg))
    coords = c.to_string('hmsdms')
    coords = coords.replace('h',':').replace('d',':').replace('m',':').replace('s','')
    
    # return RA and DEC in "hh:mm:ss.ssss dd:mm:ss.ssss" form
    return coords
    
    
def grab_pulsaralog(jlist=None):
    """
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
    Calculate the FWHM for the beam assuming ideal response/gaussian-like profile. Will eventually be depricated.
    """
    c = 299792458.0                 # speed of light (m/s)
    Dtile = 4.0                     # tile size (m) - incoherent beam
    freq = freq * 1e6               # convert from MHz to Hz
    fwhm = 1.2 * c / (Dtile * freq) # calculate FWHM using standard formula
 
    return fwhm


def getmeta(service='obs', params=None):
    """
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
  
    
def calcbox( ras, dec, time):
    #calcbox( params, prevRA, final):
    """
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
    
    """
    # calculate the end RA 
    sDEC = eDEC = startDEC
    if prevRA == None:
        # no previous RA, calculate where the center of the beam ends up
        # correct for the edge of the beam (+0.5*beam size)
        endRA = calcEndRA( float(startRAcenter), duration, sDEC, eDEC)
    else:
        # have previous RA (starting from end RA of last file)
        # the center of the beam is then prevRA-0.5*beam size
        # correct for edge of the beam (+0.5*beam size)
        endRA = calcEndRA( (float(prevRA) - (size / 2.0)), duration, sDEC, eDEC) + (size / 2.0)
    """
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
    Creates a box and a cicle at the start and end of the obs then checks if the source is inside it
    """
    found = False
    #TODO: (BWM) this is just a coding prferences, but try to space out big chunks of code.
    source = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')
    obs, ra_beam, dec_beam, time_beam = beam_input #in degrees from metadata
    ra_intial, ra_final, dec_top, dec_bot = calcbox( ra_beam, dec_beam, time_beam)
    
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
    outputfile = str(args.output)
    os.system( 'rm -f ' + outputfile + str(obs) + '_circle_beam.txt')
    with open(outputfile + str(obs) + '_circle_beam.txt',"wb") as out_list:
        out_list.write('All of the sources found within a simple 10 degree radius '\
        + 'circle from the centre of the beam for observation ID: ' + str(obs) + '\n')
        for row in temp:
             out_list.write(row + "\n")
    print "A list of sources for the observation ID: " + str(obs) + \
              " has been output to the text file: " + str(obs) + '_circle_beam.txt'
    return


def circlebeam_obsforsource(beam, coord1, coord2, sourcelist, names):
    """
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
                temp.append(name)
                print str(name)
            
        #write a simple file with the list of obs for each source that 
        #can be expanded later to include their paramaters
        outputfile = str(args.output)
        os.system( 'rm -f ' + outputfile + str(name) + '_circle_beam.txt')
        with open(outputfile + str(name) + '_circle_beam.txt',"wb") as out_list:
            out_list.write('All of the observation IDs that found ' + str(name) +\
            ' within a simple 10 degree radius circle from the centre of the beam.\n' +\
            'RA(deg): ' + b[1] + '   DEC(deg): ' + b[2] + '   Durcation(s): ' + b[3] + '\n')
            for row in temp:
                 out_list.write(str(row) + "\n")
        print "A list of observation ID for the source: " + str(name) + \
              " has been output to the text file: " + str(name) + '_circle_beam.txt'                   
    return


parser = argparse.ArgumentParser(description="""
This code is used to list the sources within the beam of observations IDs or using --obs_for_source list all the observations for each source. The sources can be input serval ways: using a list of pulsar names (--pulsar), using a complete catalogue file of pulsars (--dl_PSRCAT) or RRATs (--RRAT and --dl_RRAT), using a compatable catalogue (--in_cat with the help of --names and --coordstype) or using a RA and DEC coordinate (--coords). The observation IDs can be input (--obsid) or gathered from a directory (--FITS_dir). The default is to search all observation IDs from http://mwa-metadata01.pawsey.org.au/metadata/ that have voltages and list every known pulsar from PSRCAT in each observation ID.
""")
parser.add_argument('--obs_for_source',action='store_true',help='Instead of listing all the sources in each observation it will list all of the observations for each source.')
parser.add_argument('--output',type=str,help='Chooses a file for all the text files to be output to. To use the current file use "./". The default is /scratch2/mwaops/pulsar/incoh_census/circle_beam_output/', default = '/scratch2/mwaops/pulsar/incoh_census/circle_beam_output/')
parser.add_argument('-b','--beam',type=str,help='Decides the beam approximation that will be used. Options: "c" a simple circular beam approximation of radius 10 degrees. A more accurate beam model will be implimented at a later date. Default: circle beam')
#impliment an option for the accurate beam later on and maybe the old elipse approximation if I can make it accurate

#source options
sourargs = parser.add_argument_group('Source options', 'The different options to control which sources are used. Default is all known pulsars.')
sourargs.add_argument('-p','--pulsar',type=str, nargs='*',help='Searches for all known pulsars. This is the default. To search for individual pulsars list their Jnames in the format " -p J0534+2200 J0630-2834"')
sourargs.add_argument('--RRAT',action='store_true',help='Searches for all known RRATs.')
#Eventually impliment to search for FRBs and a search for all mode
sourargs.add_argument('--dl_RRAT',action='store_true',help='Download the RRATalog from http://astro.phys.wvu.edu/rratalog/ and uses this as the source catalogue.')
sourargs.add_argument('--dl_PSRCAT',action='store_true',help='Download the Puslar alog from http://www.atnf.csiro.au/research/pulsar/psrcat/ and uses this as the source catalogue.')
sourargs.add_argument('--in_cat',type=str,help='Location of source catalogue, must be readable by astropy.table.Table (i.e. csv, txt, votable, fits) . Default: for pulsars pulsaralog.csv from grab_pulsaralog.py and for RRATs rratalog.csv from grab_RRATalog.py')
sourargs.add_argument('--names',type=str,help='String containing the column name for the source names in the input catalogue (--in_cat). If there is no such column, use: --names=-1 and the output text file will be labelled using the coordinates in degrees: <longitudinal>_<latitudinal>.txt. Default: "Jname".')
sourargs.add_argument('--coordstype',type=str,help='String containing the two column labels of the source coordinates for the input catalouge (--in_cat). i.e.: "x,y" or "long,lat". If not provided, assumes that the coordinates are "Raj,Decj". Must be enterered as: "coord1,coord2".')
sourargs.add_argument('-c','--coords',type=str,help='String containing coordinates in "RA,DEC". This will list the OBS IDs that contain this coordinate. Must be enterered as: "hh:mm:ss.ss,+dd:mm:ss.ss".')
#finish above later and make it more robust to incclude input as sex or deg and perhaps other coordinte systmes

#observation options
obargs = parser.add_argument_group('Observation ID options', 'The different options to control which observation IDs are used. Default is all observation IDs with voltages.')
obargs.add_argument('--FITS_dir',type=str,help='Location of FITS files on system. If not chosen will search the database for metadata.')
obargs.add_argument('-o','--obsid',type=str,nargs='*',help='Input several OBS IDs in the format " -o 1099414416 1095506112". If this option is not input all OBS IDs that have voltages will be used')
obargs.add_argument('--all_volt',action='store_true',help='Includes observation IDs even if there are no raw voltages in the archive. Some incoherent observation ID files may be archived even though there are raw voltage files. The default is to only include files with raw voltage files.')
args=parser.parse_args()
  
if args.dl_RRAT:
    os.system('python /scratch2/mwaops/pulsar/incoh_census/beam_code/grab_RRATalog.py')

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
        else:
            #converts the list of pulsars into a string so they ccan be used as an agrument
            jlist = args.pulsar
            grab_pulsaralog(jlist)
            catDIR = 'temp.csv'
    else:
        catDIR = 'pulsaralog.csv'

#defaults for the coords types
if args.coordstype:
    c1, c2 = args.coordstype.split(',')
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
if args.names and args.names=='-1':
    name_col='-1'
elif args.names:
    name_col = args.names
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
    OBSID = []
    print "Obtaining list of observation IDs that have recorded voltages from http://mwa-metadata01.pawsey.org.au/metadata/"
    temp = getmeta(service='find', params={'mode':'VOLTAGE_START'})
    for row in temp:
        OBSID.append(row[0])

#for source in obs
#gets all of the basic meta data for each observation ID
cord = []
for ob in OBSID:
    print "Obtaining metadata from http://mwa-metadata01.pawsey.org.au/metadata/ for OBS ID: " + str(ob)
    beam_meta_data = getmeta(service='obs', params={'obs_id':ob})
    #check if there is a .dat file and therefore recorded voltages
    filedata = beam_meta_data[u'files']
    keys = filedata.keys()
    check = False 
    for k in keys:
        if '.dat' in k:
            check = True
    if check or args.all_volt:
        #if check succeedes then records the paramaters      
        ra = beam_meta_data[u'metadata'][u'ra_pointing']
        dec = beam_meta_data[u'metadata'][u'dec_pointing']
        time = beam_meta_data[u'stoptime'] - beam_meta_data[u'starttime'] #gps time
        
        #instead of downloading all of the obs id first, if not in obs_for_source mode, 
        #downlads one obs at a time
        if args.obs_for_source:
            cord.append([ob,ra,dec,time]) #in degrees
        else:
            cord = [ob,ra,dec,time]
            circlebeam(cord, c1, c2, catalog, name_col)
    else:
        print "No raw voltage file"

#chooses the beam type and whether to list the source in each obs or the obs for each source
#more options will be included later
if args.obs_for_source:
    if args.beam == 'c':
        circlebeam_obsforsource(cord, c1, c2, catalog, name_col)
    else:
        circlebeam_obsforsource(cord, c1, c2, catalog, name_col)



print "The code is complete and all results have been output to text files"   

                 
"""
Program tests
python find_pulsar_in_obs.py -p --obsid 1099414416
python find_pulsar_in_obs.py -p J0534+2200 J0538+2817 --obsid 1099414416
python find_pulsar_in_obs.py -p J0630-2834 --obsid 1067285064 10689221844 1101491208 1102270216
python find_pulsar_in_obs.py -p J0534+2200 J0538+2817 --obsid 1099414416 --obs_for_source
python find_pulsar_in_obs.py -c 05:34:31.973,+22:00:52.06 --obsid 1099414416 --obs_for_source               

Found in my third year project
J0437-4715 J0534+2200 J0630-2834 J0742-2822 J0835-4510 J0953+0755 J1731-4744 J1752-2806 J1900-2600 J1921+2153 J1932+1059 J1935+1616 J2048-1616 J2145-0750                 
"""
