#!/usr/bin/env python

import logging, sys, os, glob, subprocess, string, re, urllib, math, time
from optparse import OptionParser,OptionGroup
import datetime,pytz

from mwapy import ephem_utils

# configure the logging
logging.basicConfig(format='# %(levelname)s:%(name)s: %(message)s')
logger=logging.getLogger('timeconvert')
logger.setLevel(logging.WARNING)


################################################################################
def main():

    usage="Usage: %prog [options] <file1> <file2>\n"
    usage+='\tConverts between GPS seconds, UT time, MJD\n\n'    
    usage+="""\tExample:
    \tplock[~/mwa]% bin/timeconvert.py --year=2012 --month=9 --day=30
    \tEntered 2012-09-30 00:00:00...
    \t2012-09-30 00:00:00 UTC
    \t20120930000000
    \tGPS 1032998416
    \tMJD 56200.00000
    \tJD 2456200.50000
    """

    parser = OptionParser(usage=usage)
    parser.add_option('--gps',dest='gps',default=None,
                      help='GPS seconds')
    parser.add_option('--datetime',dest='datetime',default=None,
                      help='Datetime string YYYYMMDDhhmmss (UTC)')
    parser.add_option('--mjd',dest='mjd',default=None,
                      help='MJD')
    parser.add_option('--year',dest='year',default=None,
                      help='Year')
    parser.add_option('--month',dest='month',default=None,
                      help='Month')
    parser.add_option('--day',dest='day',default=None,
                      help='Day')
    parser.add_option('--hour',dest='hour',default=None,
                      help='Hour')
    parser.add_option('--minute',dest='minute',default=None,
                      help='Minute')
    parser.add_option('--second',dest='second',default=None,
                      help='Second')

    (options, args) = parser.parse_args()
    t=None
    if options.gps is not None:
        try:            
            t=ephem_utils.MWATime(gpstime=float(options.gps))
            print "Entered GPStime=%s..." % options.gps
        except:
            logger.error('Error converting gps time %s\n\t%s' % (options.gps,sys.exc_info()[1]))
            sys.exit(0)
    elif options.datetime is not None:        
        try:
            t=ephem_utils.MWATime()
            t.datetimestring=options.datetime
            print "Entered datetime string=%s..." % options.datetime
        except:
            logger.error('Error converting datetime string %s\n\t%s' % (options.datetime,sys.exc_info()[1]))
            sys.exit(0)
    elif options.mjd is not None:
        try:
            t=ephem_utils.MWATime()
            t.MJD=float(options.mjd)
            print "Entered MJD=%s..." % options.mjd
        except:
            logger.error('Error converting MJD %s\n\t%s' % (options.mjd,sys.exc_info()[1]))
            sys.exit(0)
    elif options.year is not None and options.month is not None and options.day is not None:
        if options.hour is None:
            options.hour=0
        if options.minute is None:
            options.minute=0
        if options.second is None:
            options.second=0
        try:
            t=ephem_utils.MWATime(year=int(options.year),
                                  month=int(options.month),
                                  day=int(options.day),
                                  hour=int(options.hour),
                                  minute=int(options.minute),
                                  second=int(options.second))
            print "Entered %s..." % t.strftime('%Y-%m-%d %H:%M:%S')
        except:
            logger.error('Error converting time\n\t%s' % (sys.exc_info()[1]))
            sys.exit(0)
    if t is not None:
        print t.strftime('%Y-%m-%d %H:%M:%S %Z')
        print t.strftime('%Y%m%d%H%M%S')
        print 'GPS %d' % t.gpstime
        print 'MJD %.5f' % (t.MJD)
        print 'JD %.5f' % (t.MJD+2400000.5)
        print '%s LST' % (t.LST.strftime('%H:%M:%S'))
    sys.exit(0)
    
                            


######################################################################
if __name__=="__main__":
    main()
    
    
