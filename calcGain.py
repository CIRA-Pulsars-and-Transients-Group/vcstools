#!/usr/bin/env python

import numpy as np
import sys
from astropy.constants import c,k_B

try:
	fname = sys.argv[1]
	freq = sys.argv[2]
	eta = sys.argv[3]
	tres = float(sys.argv[4])
	pres = float(sys.argv[5])
except:
	print "!!! Error: need 5 arguments !!!"
	print "filename (full path is best)"
	print "frequency (in Hz)"
	print "efficiency (as a fraction of 1)"
	print "ZA resoluion element (in degrees)"
	print "Azimuth resolution element (in degrees)"
	sys.exit(0)


omega_A = 0
with open(fname,"r") as f:
	for i in xrange(14):
		f.next()
	for row in f:
		words = row.split()
		theta,b = float(words[0]),float(words[-1])
		bsintheta = b*np.sin(np.radians(theta))

		omega_A += bsintheta * tres * pres

print "Omega_A = {0} sr".format(omega_A)

# the effective area is then
#  A_e = eta*(lambda**2/Omega_A)
print "calculating effective area"
area = eta*((c.value/freq)**2/omega_A)
print "A_e = {0} m^2".format(area)

# the gain is then A_e/2k_B * 1e-26 to convert to K/Jy
gain = (1e-26)*area/(2*k_B.value)
print "Gain: {0:.6f} K/Jy".format(gain)
