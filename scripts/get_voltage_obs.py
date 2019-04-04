#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import urllib
import urllib2
import json

import time

from mwa_metadb_utils import getmeta 

# Append the service name to this base URL, eg 'con', 'obs', etc.
BASEURL = 'http://ngas01.ivec.org/metadata/'



if __name__ == '__main__':
    start = 0
    obsinfo = getmeta(service='find', params={'mode':"VOLTAGE_START",'limit':400})
    for entry in obsinfo:
        if (entry[0] > start):
            fileinfo = getmeta(service='obs',params={'obs_id':str(entry[0]),'filetype':11})
            log_name = "report.log"
            with open(log_name, 'a') as log:
                string =  "ObsID: " + str(entry[0]) + " ObsName: "  + str(entry[1]) + " NumVCSFiles: " + str(len(fileinfo['files'])) + "\n"
                print string
                log.write(string)
            time.sleep(2)

    
  
