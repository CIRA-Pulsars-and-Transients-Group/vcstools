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
from mwa_metadb_utils import getmeta 

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
    parser.add_argument('-o','--obsid',type=str,help="Input Observation ID")
    parser.add_argument('-l','--list',type=str,help="Input Observation list file name (txt)")
    args=parser.parse_args()

    if args.obsid:
        ob = args.obsid
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
