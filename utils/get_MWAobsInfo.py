#! /usr/bin/env python

"""
Script to request obs info from MWA metadata database. 
The output .txt file from find_pulsar_in_obs.py --obs_for_source mode can be directly use as input for this script

Written by: Mengyao Xue
Originally created: 12 Sep 2018
Updated: 20 Nov 2018
"""

import argparse
import urllib
import urllib2
import json
import astropy
import os

def getmeta(service='obs', params=None):
    """
    Given a JSON web service ('obs', find, or 'con') and a set of parameters as
    a Python dictionary, return the RA and Dec in degrees from the Python dictionary.

    getmeta(service='obs', params=None)
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


def get_singleobs_meta(ob):
    beam_meta_data = getmeta(service='obs', params={'obs_id':ob})
    obn = beam_meta_data[u'obsname']
    ra = beam_meta_data[u'metadata'][u'ra_pointing']
    dec = beam_meta_data[u'metadata'][u'dec_pointing']
    dura = beam_meta_data[u'stoptime'] - beam_meta_data[u'starttime'] #gps time
    #dura = beam_meta_data[u'shifttime']
    mode = beam_meta_data[u'mode']
    gridpoint_n = beam_meta_data[u'metadata'][u'gridpoint_number']
    Tsky = beam_meta_data[u'metadata'][u'sky_temp']
    xdelays = str(beam_meta_data[u'rfstreams'][u"0"][u'xdelays']).replace(' ','')

    cord=[]#cord.append([ob,ra,dec,dura]) #in degrees
    minfreq = float(min(beam_meta_data[u'rfstreams'][u"0"][u'frequencies']))
    maxfreq = float(max(beam_meta_data[u'rfstreams'][u"0"][u'frequencies']))
    centrefreq = 1.28 * (minfreq + (maxfreq-minfreq)/2)
    return [obn,ra,dec,dura,centrefreq,Tsky,xdelays]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-o','--obid',type=str,help="Input Observation ID")
    parser.add_argument('-l','--list',type=str,help="Input Observation list file name (txt)")
    args=parser.parse_args()

    if args.obid:
        ob = args.obid
        [obn,ra,dec,dura,centrefreq,Tsky,xdelays]=get_singleobs_meta(ob)
        print "{0:<10}   {1:<15}   {2:<6}   {3:<5}{4} {5:<6}  {6}".format('Obs_ID', 'Obs_name', 'RA', 'Dec', 'Duration', 'CentreFreq', 'Tsky')
        print "{0}   {1:<15}   {2:6.2f}   {3:6.2f}   {4:>4}   {5}   {6:7.3f}".format(ob,obn,ra,dec,dura,centrefreq,Tsky)
        print "Delays Array"
        print xdelays
    elif args.list:
        txtfile = open(args.list).readlines()
        cord=[]
        print "{0:<10}   {1:<15}   {2:<6}   {3:<5}{4} {5:<6}  {6}".format('Obs_ID', 'Obs_name', 'RA', 'Dec', 'Duration', 'CentreFreq', 'Tsky')
        for line in txtfile:
            if line[0]!=('1'):
                continue
            words=line.split()
            ob=words[0]
            [obn,ra,dec,dura,centrefreq,Tsky,xdelays]=get_singleobs_meta(ob)
            print "{0}   {1:<15}   {2:6.2f}   {3:6.2f}   {4:>4}   {5}   {6:7.3f}".format(ob,obn,ra,dec,dura,centrefreq,Tsky)
    else:
        print "Please input obsevation id by setting -o or -l "
