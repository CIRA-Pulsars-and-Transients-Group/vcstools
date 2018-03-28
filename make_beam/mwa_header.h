/***************************************************************************
 *   
 *   Based on cpsr2_header.h Copyright (C) 2002 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#ifndef __MWA_HEADER_h
#define __MWA_HEADER_h

#define MWA_HEADER_SIZE 4096
#define MWA_HEADER_VERSION 0.1

/* ************************************************************************

   The MWA header will be stored in the shared memory of the twopip machines
   and also at the beginning of each block of data written to disk. It will separate 
   data taken from each fiber (signal) thus allowing writes to disk to be concatenated.
   
   There is an individual header for each receiver - therefore each node will have 2 and the 
   2PiP system will have 4. 

   As 2PiP is currently an MPI process the nodes will need to know the keys for each command header

   It is an ASCII formatted list of fields written
   in "key value" pairs; each pair is separated by one or more newline
   characters.  Comments should be proceeded by '#' characters.  The last 
   string in the header should be the string, "DATA", followed by a single
   newline character. The size of data block in bytes needs to be known and is an essential
   element of this header as another header will be located at that point in the data stream.
   
   It is important that keywords (OBSERVATION for example) are not
   replicated in any other part of the file: in the beta version of
   the header reader, we simply examine the header for the first
   occurance of a keyword.  Hiding the replication behind a comment
   symbol is NOT acceptable.

   The header should contain at least the following information:
*/

#define MWA_HEADER_INIT \
"HDR_VERSION 0.1	# Version of this ASCII header\n" \
"MWA_CAPTURE_VERSION 0.1	# Version of the Data Acquisition Software\n" \
"MWA_SAMPLE_VERSION 0.1	# Version of the FFD FPGA Software\n" \
"\n" \
"VCSTOOL_VERSION " VERSION_BEAMFORMER "\n" \
"\n" \
"TELESCOPE	mwa		# telescope name\n" \
"\n" \
"SOURCE		unset		# name of the astronomical source\n" \
"RA		unset		# Right Ascension of the source\n" \
"DEC		unset		# Declination of the source\n" \
"\n" \
"FREQ		unset		# centre frequency on sky in MHz\n" \
"BW		unset		# bandwidth of in MHz (-ve lower sb)\n" \
/*
  

  Programs should initialize a cpsr2 header as follows:

  char mwa_header[MWA_HEADER_SIZE] = MWA_HEADER_INIT;

  You can also:

  strcpy (buffer, MWA_HEADER_INIT);

  It is recommended to use the "ascii_header_set/get" routines in
  order to manipulate the CPSR2 header block.  See
  test_cpsr2_header.c, or for example:


*/

#include "ascii_header.h"

#endif
