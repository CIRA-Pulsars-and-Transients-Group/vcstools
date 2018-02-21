/***************************************************************************
 *   
 *   Based on cpsr2_header.h Copyright (C) 2002 by Willem van Straten
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#ifndef __TWOPIP_HEADER_h
#define __TWOPIP_HEADER_h

#define TWOPIP_HEADER_SIZE 4096
#define TWOPIP_HEADER_VERSION 0.1

/* ************************************************************************

   The TWOPIP header will be stored in the shared memory of the twopip machines
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

#define TWOPIP_HEADER_INIT \
"TWOPIP_HEADER_VERSION 0.1	# Version of this ASCII header\n" \
"TWOPIP_CAPTURE_VERSION 0.1	# Version of the Data Acquisition Software\n" \
"TWOPIP_SAMPLE_VERSION 0.1	# Version of the FFD FPGA Software\n" \
"\n" \
"TELESCOPE	unset		# telescope name\n" \
"PRIMARY	unset		# primary node host id (0 | 1)\n" \
"UNIT	        unset		# primary node capture device (0 | 1)\n" \
"RECEIVER	unset		# receiver id (0..63) - as reported by the packets \n" \
"\n" \
"# time of the rising edge of the first time sample\n" \
"UTC_START	unset		# yyyy-mm-dd-hh:mm:ss.fs\n" \
"MJD_START	unset		# MJD equivalent to the start UTC\n" \
"\n" \
"OFFSET		unset		# bytes offset from the start MJD/UTC\n" \
"BLKSIZE	unset		# bytes in this block \n" \
"\n" \
"SOURCE		unset		# name of the astronomical source\n" \
"RA		unset		# Right Ascension of the source\n" \
"DEC		unset		# Declination of the source\n" \
"\n" \
"FREQ		unset		# centre frequency on sky in MHz\n" \
"BW		unset		# bandwidth of in MHz (-ve lower sb)\n" \
"CHBW		unset		# channel bandwidth \n" \
"TSAMP		unset		# sampling interval in microseconds\n" \
"NBIT		unset		# number of bits per sample\n" \
"NTILE		unset		# number of channels \n"
/*
  

  Programs should initialize a cpsr2 header as follows:

  char cpsr2_header[CPSR2_HEADER_SIZE] = CPSR2_HEADER_INIT;

  You can also:

  strcpy (buffer, CPSR2_HEADER_INIT);

  It is recommended to use the "ascii_header_set/get" routines in
  order to manipulate the CPSR2 header block.  See
  test_cpsr2_header.c, or for example:

  ------------------------------------------------------------------------

  char cpsr2_header[CPSR2_HEADER_SIZE] = CPSR2_HEADER_INIT;

  char* telescope_name = "parkes";
  ascii_header_set (cpsr2_header, "TELESCOPE", "%s", telescope_name);

  float bandwidth = 64.0;
  ascii_header_set (cpsr2_header, "BW", "%f", float);

  [...]

  double centre_frequency;
  ascii_header_get (cpsr2_header, "FREQ", "%lf", &centre_frequency);

*/

#include "ascii_header.h"

#endif
