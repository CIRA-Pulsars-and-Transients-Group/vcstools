#!/usr/bin/env python

from mwa_metadb_utils import getmeta 

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def print_minmax(obs_id):
    """
    Snippet function to show how to get the first and last files for a given observation
    """

    obsinfo = getmeta(service='obs', params={'obs_id':str(obs_id)})
    times=[file[11:21] for file in obsinfo['files'] if is_number(file[11:21])] #Make a list of gps times excluding non-numbers from list
    obs_start = int(min(times))
    obs_end = int(max(times))

    print "{0} file found with Observation ID {1}".format(len(obsinfo['files']), obs_id)
    print "First file: {0}".format(obs_start)
    print "Last file: {0}".format(obs_end)

if __name__ == '__main__':
    from sys import argv
    print_minmax(argv[1])




