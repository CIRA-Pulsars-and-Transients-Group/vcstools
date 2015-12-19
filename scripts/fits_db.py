#!/usr/bin/env python

import psycopg2
import sys

time =sys.argv[1] 
#time = str(1086725208)


conn = psycopg2.connect("dbname=mwa host=mwa-metadata01.pawsey.org.au user=MWA-guest password=guest port=5432")

cur = conn.cursor()

cur.execute("SELECT rf_stream.frequencies, rf_stream.ra, rf_stream.dec, rf_stream.azimuth, rf_stream.elevation FROM public.mwa_setting, public.rf_stream WHERE"
            " mwa_setting.starttime = rf_stream.starttime and mwa_setting.starttime = %s" % (time))
params = cur.fetchone()
freqs = params[0]



if freqs == range(freqs[0], freqs[0]+len(freqs)):
    centre = freqs[len(freqs)/2]
    print "Center Frequency ", str((centre * 1.28) - 0.64)
else:
    centre =  None
    print "Non contiguous frequency channels in observation"


#freqs = [str(freq) for freq in freqs]
#print ','.join(freqs)
