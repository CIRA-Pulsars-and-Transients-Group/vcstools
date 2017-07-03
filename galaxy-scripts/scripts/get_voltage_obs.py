#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import urllib
import urllib2
import json

import time

# Append the service name to this base URL, eg 'con', 'obs', etc.
BASEURL = 'http://ngas01.ivec.org/metadata/'


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

    
  
