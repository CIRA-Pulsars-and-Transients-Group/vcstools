#!/usr/bin/env python

import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import sys

fname = sys.argv[1]
print "loading data"
theta,phi,beam = np.genfromtxt(fname,comments=("#","*"),skip_header=14,usecols=(0,1,8),unpack=True)
print "done"

print "plotting in polar coordinates"
fig = plt.figure()
ax = fig.add_subplot(111,polar=True)
ax.tricontourf(np.radians(phi),np.radians(theta),beam)
plt.savefig(fname.replace(".dat",".png"),bbox_inches='tight')

#BASIC TEST TO MAKE SURE DATA ARE LOADED ETC
#fig = plt.figure()
#ax = fig.add_subplot(111,polar=True)
#ax.scatter(theta,phi,c=beam,edgecolors='none')
#plt.show()




