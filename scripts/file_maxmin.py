#!/usr/bin/env python

from mwa_metadb_utils import getmeta 

def print_minmax(obs_id):
    """
    Snippet function to show how to get the first and last files for a given observation
    """

    obsinfo = getmeta(service='obs', params={'obs_id':str(obs_id)})
    times=[file[11:21] for file in obsinfo['files'] if file[11:21].isnumeric()] #Make a list of gps times excluding non-numbers from list
    obs_start = int(min(times))
    obs_end = int(max(times))
    obs_dur = len(obsinfo['files'])

    return obs_start, obs_end, obs_dur

if __name__ == '__main__':
    from sys import argv
    obs_start, obs_end, obs_dur = print_minmax(argv[1])
    print("{0} file found with Observation ID {1}".format(obs_dur, argv[1]))
    print("First file: {0}".format(obs_start))
    print("Last file: {0}".format(obs_end))
    




