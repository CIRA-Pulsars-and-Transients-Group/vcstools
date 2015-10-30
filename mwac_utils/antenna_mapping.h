/*
 * antenna_mapping.h
 *
 *  Created on: May 5, 2012
 *      Author: sord
 */

#ifndef ANTENNA_MAPPING_H_
#define ANTENNA_MAPPING_H_

#define NINPUT 256
/* Each pfb maps antennas in the following way:
 *
 *   antenna number is found as index
	 0 = 0
	 1 = 16
	 2 = 32
	 3 = 48
	 4 = 1
	 5 = 17
	 6 = 33
	 7 = 49
	 .
	 .

	 which can be described by

	 e.g

	 floor(index/4) + index%4*16 = antenna

	 e.g index 0

	 0 + 0 = 0;

	 index 1

	 0 + 1*16 = 16

	 index 2

	 0 + 2*16 = 32


	 index 3

	 0 + 3*16 = 48


	 index 4

	 1 + 0 = 1

	 index 5

	 1 + 1*16 = 17

	 etc ....

 *
 */


typedef struct mapping {

	int stn1;
	int stn2;
	int pol1;
	int pol2;

} map_t;

extern map_t corr_mapping[NINPUT][NINPUT];
extern int pfb_output_to_input[NINPUT];
extern int single_pfb_mapping[64];
extern int miriad_to_mwac[256];

#endif /* ANTENNA_MAPPING_H_ */
