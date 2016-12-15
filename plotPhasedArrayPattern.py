#!/usr/bin/env python

import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import sys

fname = sys.argv[1]
table = np.genfromtxt(fname,comments="#",skip_header=14,usecols=(0,1,8),max_rows=1e6)
theta = table[0]
phi = table[1]
power = table[2]

