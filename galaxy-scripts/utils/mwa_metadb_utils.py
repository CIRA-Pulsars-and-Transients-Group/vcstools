#!/usr/bin/env python

"""
Collection of database related utilities that are used throughout the VCS processing pipeline
"""

def getmeta(service='obs', params=None):
    """
    Function to call a JSON web service and return a dictionary:
    Given a JSON web service ('obs', find, or 'con') and a set of parameters as
    a Python dictionary, return a Python dictionary xcontaining the result.
    Taken verbatim from http://mwa-lfd.haystack.mit.edu/twiki/bin/view/Main/MetaDataWeb
    """
    import urllib
    import urllib2
    import json

    # Append the service name to this base URL, eg 'con', 'obs', etc.
    BASEURL = 'http://mwa-metadata01.pawsey.org.au/metadata/'


    if params:
        # Turn the dictionary into a string with encoded 'name=value' pairs
        data = urllib.urlencode(params)
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


def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def obs_max_min(obs_id):
    """
    Small function to query the database and return the times of the first and last file
    """
    obsinfo = getmeta(service='obs', params={'obs_id':obs_id})
    
    # Make a list of gps times excluding non-numbers from list
    times = [f[11:21] for f in obsinfo['files'].keys() if is_number(f[11:21])]
    obs_start = int(min(times))
    obs_end = int(max(times))
    return obs_start, obs_end


