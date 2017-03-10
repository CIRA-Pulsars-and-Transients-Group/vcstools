#!/usr/bin/python
"""
Contains a single function that reads in a list of integers (channel numbers),
and returns a list of the same integers whose order depends on whether or not
any of the channels are above 128.
"""

def sfreq(freqs):

    if len(freqs) != 24:
        print "There are not 24 coarse chans defined for this obs. Got: %s" % freqs
        return

    #freqs.sort()   # It should already be sorted, but just in case...[SET] Commenting this out because sort() is ironically putting 2-digit channels out of order
    lowchans = [f for f in freqs if f <= 128]
    highchans = [f for f in freqs if f > 128]

    highchans.reverse()
    freqs = lowchans + highchans

    #print "lowchans", lowchans
    #print "highchans", highchans
    #print "freqs", freqs

    return freqs

