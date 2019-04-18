#!/usr/bin/env python

from mwa_metadb_utils import getmeta 

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def print_info(obs_id):
    """
    Snippet function to show how to get the first and last files for a given observation
    """

    obsinfo = getmeta(service='obs', params={'obs_id':str(obs_id)})
    
    print "Obs ID:", obs_id
    print "Name:", obsinfo['obsname']
    print "Channels:", obsinfo['rfstreams']['0']['frequencies']
    print "Duration:", obsinfo['stoptime'] - obsinfo['starttime'], "seconds"
    #times=[file[11:21] for file in obsinfo['files'] if is_number(file[11:21])] #Make a list of gps times excluding non-numbers from list
    #obs_start = int(min(times))
    #obs_end = int(max(times))

    #print "{0} file found with Observation ID {1}".format(len(obsinfo['files']), obs_id)
    #print "First file: {0}".format(obs_start)
    #print "Last file: {0}".format(obs_end)

if __name__ == '__main__':
    from sys import argv
    print_info(argv[1])




