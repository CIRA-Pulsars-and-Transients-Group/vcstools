#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import os
import glob

required_size=253440000
files_glob = "*.dat"
for file in sorted(glob.glob(files_glob)):
    statinfo = os.stat(file)
    if (statinfo.st_size != required_size):
        print "Bad file size for %s\n" % file


