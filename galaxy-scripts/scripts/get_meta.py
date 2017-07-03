#!/usr/bin/env python

import urllib
import urllib2
import json

# Append the service name to this base URL, eg 'con', 'obs', etc.
BASEURL = 'http://mwa-metadata01.pawsey.org.au/metadata/'


def getmeta(service='obs', params=None):
  """
  Function to call a JSON web service and return a dictionary:
  Given a JSON web service ('obs', find, or 'con') and a set of parameters as
  a Python dictionary, return a Python dictionary containing the result.
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
    from sys import argv

    obsinfo = getmeta(service='obs', params={'obs_id':str(argv[1])})
    print "R.A.: ", obsinfo['metadata']['ra_pointing']
    print "Dec.: ", obsinfo['metadata']['dec_pointing']
    print "Azimuth: ", obsinfo['metadata']['azimuth_pointing']
    print "Zenith: ", str(90 - obsinfo['metadata']['elevation_pointing'])
    #print sys.argv[1]




