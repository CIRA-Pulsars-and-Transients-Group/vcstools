#!/usr/bin/env python

"""
Script to work out flagging of channels in PRESTO and PSRCHIVE formats.

Created: Bradley Meyers
Version: 0.1
Date: 4 Nov 2016
"""


import sys
import argparse
import numpy as np


def rfifind_zap(max,nfine,nedges,zaps):
	"""
	Write the channels to be zapped in a PRESTO readable format.
	Edge channel removal will be done using ranges, individual channels will be prepended to the list.
	"""
	if (nedges is None) and (zaps is not None):
		return ",".join(zaps.astype(str)),len(zaps)
	elif (nedges is None) and (zaps is None):
		return "",0

	s = ""
	num = 0
	nintervals = 0
	while num < max:
		if num != 0: s = s + ","
		s = s + str(num) + ":" + str(num + (nedges-1))
		num += nfine-nedges
		s = s + "," + str(num) + ":" + str(num + (nedges-1))
		num += nedges
		nintervals += 2

	if zaps is None:
		
		return s,nintervals*nedges
	else:
		return ",".join(zaps.astype(str)) + "," + s,len(zaps)+nintervals*nedges


def paz_zap(max,nfine,nedges,zaps):
        """
        Write a kill file that can be fed into paz with the "-k" option
        """
	if (nedges is None) and (zaps is not None):
		return zaps.astype(int),len(zaps)
	elif (nedges is None) and (zaps is None):
		return [],0

        start = 0
        end = nedges
        chanlist = []
	nintervals = 0
        #chanlist.extend(np.arange(start,end))
        while end < max:
		if start == 0:
			chanlist.extend(np.arange(start,start+nedges))
			nintervals += 2
                start += (nfine-nedges)
                end += (nfine-nedges)
                chanlist.extend(np.arange(start,end))
                start += nedges
                end += nedges
                if start > max or end > max:
                        break
                else:
                        chanlist.extend(np.arange(start,end))
		nintervals += 2
	
	if zaps is None:
		return chanlist,nintervals*nedges
	else:
        	chanlist[0:0] = zaps.astype(int)
        	return chanlist,nintervals*nedges


def zap(r,p,channels,zapedges=False,nzap=20,zaps=None,verbose=False):
	"""
	Main function to handle zap formatting and output
	"""
	
	if channels[1] is not None:
		coarse = channels[0]
		fine = channels[1]
		total = coarse * fine
		if verbose:
			print "Total number of channels: {0}".format(total)
			print "Comprised of {0} coarse channels, each with {1} fine channels".format(coarse,fine)
			print "Will flag {0} edge channels from each coarse channel".format(nzap)
	else:
		fine = 128
		total = channels[0]
		if verbose:
			print "Total number of channels: {0}".format(total)
			print "Assuming there are 128 fine channels per coarse channel..."
			print "Will flag {0} edge channels from each coarse channel".format(nzap)

	
	if zaps is not None:
		zaps = np.array(zaps)
		if verbose:
			print "User defined channels to zap:"
			print zaps


	if verbose:
		if r and not p:
			print "Only formatting output for PRESTO's rfifind routine"
		elif p and not r:
			print "Only formatting output for PSRCHIVE's paz routine -- will create a kill file called {0}".format(p)
		elif r and p:
			print "Formatting output for both rfifind and paz:"
			print "Will print rfifind formatted list to screen (also return the string)"
			print "Will write a kill file for paz called {0}, which can be used with the -k option".format(p)
		else:
			# I don't see how this could really happen, but just in case!
			print "!!WARNING!! :: Somehow you didn't provide input AND removed the defaults!?"
			print "!!WARNING!! :: Exitting on sys.exit(1)"
			sys.exit(1)


	if r and not p:
		if zapedges:
			tmpzaps,nzapped = rfifind_zap(total,fine,nzap,zaps)
		else:
			tmpzaps,nzapped = rfifind_zap(total,fine,None,zaps)
		if verbose:
			print "Zapped {0}/{1} channels ({2:.3f}%)".format(nzapped,total,nzapped*100./total)
		print "\n\n"+tmpzaps
		

	elif p and not r:
		with open(p,'w') as filename:
			if zapedges:
				chanlist,nzapped = paz_zap(total,fine,nzap,zaps)
				for item in chanlist:
					filename.write(str(item)+"\n")
			else:
				chanlist,nzapped = paz_zap(total,fine,None,zaps)
				for item in chanlist:
                                        filename.write(str(item)+"\n")
			if verbose:
				print "Zapped {0}/{1} channels ({2:.3f}%)".format(nzapped,total,nzapped*100./total)
	
	elif r and p:
		if zapedges:
                        tmpzaps,nzapped = rfifind_zap(total,fine,nzap,zaps)
                else:
                        tmpzaps,nzapped = rfifind_zap(total,fine,None,zaps)
                if verbose:
                        print "Zapped {0}/{1} channels ({2:.3f}%)".format(nzapped,total,nzapped*100./total)
                print "\n\n"+tmpzaps

		with open(p,'w') as filename:
                        if zapedges:
                                chanlist,nzapped = paz_zap(total,fine,nzap,zaps)
                                for item in chanlist:
                                        filename.write(str(item)+"\n")
                        else:
                                chanlist,nzapped = paz_zap(total,fine,None,zaps)
                                for item in chanlist:
                                        filename.write(str(item)+"\n")


	else: 
		# safer to not return any channels to zap in the case where no option is given
		print ""
		return ""


## SETUP OPTIONS ##
parser = argparse.ArgumentParser(description="""Output the correctly formated MWA band edge channels and user-defined channels to remove for aliasing/RFI excision.\
						 If no options are given, returns an empty string.""",\
					formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# there should be two options for the channels, either:
# give me the total number of fine channels (and I'll assume there's 128 fine channels per coarse channel)
# OR 
# give me the number of coarse channels and the number of fine channels and I'll use that information to then calculate the intervals appropriate.
parser.add_argument("-N","--nchans",action='store',type=int,metavar='Nchan',help="TOTAL number of frequency channels in data set. Assumes 128 fine channels  per coarse channel.",default=None)
parser.add_argument("-C","--n_coarse_fine",action='store',type=int,nargs=2,metavar="Nchan",help="Number of coarse channels followed by the number of fine channels per coarse channel. This options overrides --nchans.",default=None)

# zap the edge fine channels?
parser.add_argument("-Z","--zapedges",action='store_true',help="Zap the fine channel edges of the coarse channels.")

# if so, how many?
parser.add_argument("-n","--nzap",type=int,action='store',metavar="Nedge",help="Number of edge channels to remove from each side of a coarse channel. If given, but --zapedges is not then this argument is ignored.",default=20)

# zap the middle channels?
parser.add_argument("-m","--middle",action='store_true',help="Flag the center fine channels for each coarse channel?")

# user-defined channels to zap (will be prioritised over edge channel zapping)
parser.add_argument("-z",type=str,action='store',metavar="chan",nargs="+",help="Individual channels to zap")

# what version of output do you want? PRESTO and/or PSRCHIVE
parser.add_argument("-r","--rfifind",action='store_true',help="Output to screen the edge channels in a format readable by the PRESTO rfifind routine",default=False)
parser.add_argument("-p","--paz",action='store',metavar="filename",type=str,help="Output channels for PSRCHIVE paz routine into the given file, using paz's '-k filename' option")

# be verbose about the actions takens
parser.add_argument("-v","--verbose",action='store_true',help="Use verbose mode: will tell you each step what I'm doing, but this will mean you can't easily just pipe the ouput with back-ticks")

args = parser.parse_args()


# if no options are passed, then just return an empty string 
# (i.e. don't zap anything) and exit gracefully
if len(sys.argv) == 1:
	zap(False,False,(None,None))
	sys.exit(0)

# check to make sure that the channels variable is formatted correctly, 
# no matter the input (should be a tuple/iterable)
if args.nchans and not args.n_coarse_fine:
	channels = (args.nchans,None)
else:
	channels = args.n_coarse_fine

# check to make sure that if nzap is given but user has not elected 
# to zap edges, we zap 0 edge channels
if args.nzap and not args.zapedges:
	args.nzap = 0

## DO THE THING! ##
zap(args.rfifind,args.paz,channels,args.zapedges,args.nzap,args.z,args.verbose)
