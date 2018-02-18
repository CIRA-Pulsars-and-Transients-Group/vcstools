#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "fourbit.h"
#include "buffer_sizes.h"

/* stand alone 4 bit to 8 bit expander
 * works on a CPU using openmp. 
 *
 * (c) Steve Ord (steve.ord@gmail.com)
 */

/* probably quicker to expand to 32 bit numbers - so we get more samples per lookup
 * each 32bit number will be 4 eight bit numbers.
 *
 * In the LEDA/MWA context this represents the real and imaginary parts of 2 poloarisations
 *
 * The input to this is therefore 4 fourbit numbers. 
 *
 * So we get the index into the array using a ushort as an index
 *
 * One issue is that I <think> the 4 bit samples are 2's complement
 *
 * So I could simply shift to the right by 4 to keep the sign - but this "uses" up the 
 * dynamic range of the number - perhaps the best thing to do is:
 *
 * 
 * Actual expansion
 *
 * shift ORIGINAL left 4 and mask to get the sign get ANSWER
 *
 * AND ORIGINAL with 0x7 to remove the sign UNSIGNED
 *
 * ANSWER = UNSIGNED AND ANSWER
 *
 * This is 4 operations and no conditionals per sample - assume 10 ops in total
 *
 * So a 2GHz is 200 Msamples per core. But remember we have Ninput*BW samples. Pretty soon
 * you are in Gop land.
 *
 * I will build a lookup table that should give us a single operation for 4 samples.
 *
 * I will also divide the job over several cores to make things even quicker.
 *
 *
 * First stage is build the lookup table
 *
 *
 *
 * 
*/

static int8_t eight_bit_lookup[65536][4];

void build_eight_bit_lookup() {

    uint32_t index = 0;
    int value = 0;
    uint8_t original = 0;
    uint8_t answer = 0;
    uint8_t outval = 0;

    for (index=0; index <= 65535; index++) {
        for (outval = 0; outval < 4 ; outval++) {

            original = index >> (outval * 4);
            original = original & 0xf; // the sample
            if (original >= 0x8) {
            	value = original - 0x10;
            }
            else {
            	value = original;
            }
            answer = value & 0xff;
            eight_bit_lookup[index][outval] = answer;

        }
    }

    /* so inprinciple if we read in 4 x 4 bit samples we will get 4x8bit output samples with one lookup */

}

void float_to_4bit(int nsamps, float *in, int8_t *out) {

	float *input_ptr;
	int8_t *output_ptr;
	int8_t output_buf[2];
	int i=0,c=0;
	output_ptr = out;
	input_ptr =  in;

	for (i=0; i<nsamps; i=i+2) {
		for (c=0; c < 2; c++) {
			output_buf[c] = (int8_t) roundf(*input_ptr);
			input_ptr++;
		}
		*output_ptr = split_8bit(output_buf); // into 4bit samples
		output_ptr++;

	}
}

inline void expand_4bit(uint16_t *input, int8_t *output) {
	extern int8_t eight_bit_lookup[65536][4];
    memcpy(output,eight_bit_lookup[*input],4);
                   
}
inline int8_t split_8bit(int8_t *input) {

	/* this needs to take 2x8bit signed sample into 2 signed 4 bit numbers */

	int8_t sample = 0x0;
	int8_t output = 0x0;

	if (abs(*input) > 7) {	
		/* clipped */
		return 8;	
	}	
	else {
		output = *input & 0xf; // mask the top 4 off
	}

        input++; // increment the input sample

	if (abs(*input) > 7) {
		/* clipped */		
		return 8;
	}
	else {
		sample = 0x0;
		sample = *input << 4;
		output = output | sample;

		return output;
	}

		
}

void process_MWA_atomic_raw(void *input, void *output,size_t stride, size_t pfb) {
	/* this does the same job as the process_MWA_atomic - 	 *
	   BUT: Now 4 to 8 bit conversion is applied - this just turns the packet contents
         * No offsets need be calculated as these are 16 contig channels - but they do have to be turned
	 */

	int nstation = 64; // how many stations in a block
	int fine = 0;
	int nfine = 4;
	int grp = 0;
	int ngrp = 4;
	int time = 0;
	int ntime = 500;
	int ndim = 2;

	int8_t *out_ptr = (int8_t *) output;
	uint16_t *in_ptr = (uint16_t *) input;
	int8_t *current_out_ptr = NULL;

	

	for (grp = 0; grp < ngrp; grp++) { // the fine channel group
		for (time = 0; time < ntime; time++) { // the time sample

/* These are now set up for the 4 bit samples - compare to process_MWA_atomic */

			in_ptr = in_ptr + HEADER_WORDS; // Chop off header each packet is a single time step
			current_out_ptr = out_ptr + ((time * stride)/2); // this takes us to the beginning of the 16 channel block
			current_out_ptr += (grp * nfine * nstation * NPFB * ndim)/2; // this takes us to the beginning of the 4 channel block

/* this bit is really simple becuase we can block copy the samples now we do not have to expand them */
			for (fine = 0; fine<nfine; fine++) {// the fine channel
				// there are 64 complex 4 bit signals (128x4bit numbers therefore 64bytes) 
				current_out_ptr = current_out_ptr + ((ndim*pfb*64)/2);
				memcpy(current_out_ptr,in_ptr,64);
				in_ptr=in_ptr+32;
				current_out_ptr = current_out_ptr+64;
				
				// need to shift over the extra PFBs 
				if (pfb != (NPFB-1)) {
					// this is data from not from the last pfb - therefore an additional jump is needed
					current_out_ptr = current_out_ptr + (((NPFB - pfb - 1)*64*ndim/2));
					// so if the PFB is the last one there is no additional jump;
					// if the PFB is the penultimate one then we must shift by 1 pfb worth of output  
				}
					
			}
			in_ptr++;// checksum
		}
	}
}

void process_MWA_atomic(void *input, void *output,size_t stride, size_t pfb) {
	/* this does the same job as the process_MWA_buffer - but only does it to
	 * the 2000 packets in an atomic block.
	 *
	 * No offsets need be calculated as these are 16 contig channels - but they do have to be turned
	 */

	int station = 0; // which stations
	int nstation = 64; // how many stations in a block
	int fine = 0;
	int nfine = 4;
	int grp = 0;
	int ngrp = 4;
	int time = 0;
	int ntime = 500;
	int ndim = 2;

	int8_t *out_ptr = (int8_t *) output;
	uint16_t *in_ptr = (uint16_t *) input;
	int8_t *current_out_ptr = NULL;

	

	for (grp = 0; grp < ngrp; grp++) { // the fine channel group
		for (time = 0; time < ntime; time++) { // the time sample

			in_ptr = in_ptr + HEADER_WORDS; // Chop off header each packet is a single time step
			current_out_ptr = out_ptr + (time * stride); // this takes us to the beginning of the 16 channel block
			current_out_ptr += grp * nfine * nstation * NPFB * ndim; // this takes us to the beginning of the 4 channel block

			for (fine = 0; fine<nfine; fine++) {// the fine channel
				for (station=0;station<nstation;station = station + 2) {
					// the station gets to increment by 2 because each exapnsion
					// is 4 samples or r+i from 2 stations

					if (station == 0) {
						current_out_ptr = current_out_ptr + (ndim*pfb*64);
					}
					
					expand_4bit((uint16_t *) in_ptr, (int8_t *) current_out_ptr);

					in_ptr++;
					current_out_ptr = current_out_ptr+4;
				}
				// need to shift over the extra PFBs 
				if (pfb != (NPFB-1)) {
					// this is data from not from the last pfb - therefore an additional jump is needed
					current_out_ptr = current_out_ptr + ((NPFB - pfb - 1)*64*ndim);
					// so if the PFB is the last one there is no additional jump;
					// if the PFB is the penultimate one then we must shift by 1 pfb worth of output  
				}
					
			}
			in_ptr++;// checksum
		}
	}
}
void process_MWA_buffer(void *input, void *output) {

    /* this has to take in the input data and reorder it - stripping off the header in the process */

    /* in the MWA case the input order is:
     *
     * antenna, narrow frequency, time - so we still have to corner turn within the buffers ....
     *
     * there are blocks of 2000 packets and it is probably simpler just to unpack these 2000 packets at once
     * then run multiple unpacks in parallel using OMP.
     *
     * A Wrinkle is that these packets are in 500 packet blocks. Each contains 4 narrow channels. So the "natural" loop
     * element is 500 samples. The number of these 500 packet blocks in an individual correlator input buffer is 4 * <NPB> * <NCAPTURE>
     *
     * NPFB dictates how many tiles we have and NCAPTURE dictates how many channels we have
     *
     * In the mode where we are capturing from multiple PFB (default) then each capture machine is combining input lines for each channel group to get 
     * all the tiles in a single block. The receiver then combines all of the blocks to get all the channel blocks.
     *
     * So the output buffer (correlator input) has to be (4 X NPFB X NCAPTURE) therefore in 32T prototype mode (NPFB = 1, NCAPTURE = 8) assuming we capture 8 data lanes from the PFB
     * The required correlator order is [time][channel][station][pol][complexity]
     *
     * Which means each 500 packet block needs to be offset from the beginning of this buffer
     *
     *   
     * Another option is to unpack a whole input buffer - which already contains all of the NCAPTURE lanes - thereby removing the need to book-keep outside this function. 
     *
     * So we still unpack 500 sample blocks of narrow channel groups - but there is NPFB of them and NCAPTURE of those
     *
     * we can then run the expander on the buffer - perhaps run the expand at the same time as the re-order to
     * cut down on memory bandwidth
     *
     */

#define MWA_PFB_STATIONS 32
#define MWA_FINE_CHAN_GRPS 4    
#define MWA_CHAN_PER_GROUP 4
#define MWA_GROUPS_PER_CAPTURE 4    
#define SAMPS_PER_PACKET 128 // how many dual pol samples in a packet
#define NTIME_PER_BUFFER 500 
#define MWA_NPOL 2
#define MWA_NDIM 2

    int capture=0,pfb=0,narrow_chan =0,time=0,station_index=0;

    int8_t *out_ptr = (int8_t *) output;
    uint16_t *in_ptr = (uint16_t *) input;
    int8_t *current_out_ptr = NULL;

    size_t out_offset = 0;
    size_t pkt_num = 0;

    int fine_ch=0;


    /* The increment between output packets is the size of a 500 packet block
     * and of course which narrow channel group we are in.
     *
     */


    for (capture = 0 ; capture < MAX_NCAPTURE ; capture++ ) { // How many PFB lanes are being captured  


        for (pfb = 0 ; pfb < NPFB ; pfb++) { // how many pfb are there

            pkt_num = 0; // this will reset every 2000 packets

            for (narrow_chan = 0; narrow_chan < MWA_FINE_CHAN_GRPS; narrow_chan++) {

            	for (time = 0; time < NTIME_PER_BUFFER; time++) {

                    /* now have to calculate the output offset - where in the GPU buffer do I drop this time sample */
                    /* time is definitely the slowest variable in the GPU block - how many samples are there for each time */
                    /* == nchannel * nstation * npol
                       then the next slowest variable is channel
                       how do I know which channel I am:
                       So first we have [floor (packet count / 500)]%4 * 4 + 16*line_number

                       0/500    = 0%4 * 4 + (0/2000) * 128 + 16*line_number = 0
                       500/500  = 1%4 * 4 + (500/2000) * 128 +16*line_number = 4
                       1000/500 = 2%4 * 4 + (1000/2000) * 128 + 16*line_number= 8
                       1500/500 = 3%4 * 4 = 12
                       2000/500 = 4%4 * 4 +  = 0

*/

                    in_ptr = in_ptr + HEADER_WORDS; // Chop off header

                    for (fine_ch = 0 ; fine_ch < MWA_CHAN_PER_GROUP; fine_ch++) {


                        /* so the output offset for each time step is */

                        out_offset = time * MWA_GROUPS_PER_CAPTURE * MAX_NCAPTURE * MWA_CHAN_PER_GROUP * NPFB * MWA_PFB_STATIONS * MWA_NPOL * MWA_NDIM; 
                        /* the output offset for this channel group is */
                        out_offset = out_offset + (capture * MWA_GROUPS_PER_CAPTURE * MWA_CHAN_PER_GROUP * NPFB * MWA_PFB_STATIONS * MWA_NPOL * MWA_NDIM);
                        /* the output offset for this fine channel group is */
                        out_offset = out_offset + (narrow_chan * MWA_CHAN_PER_GROUP * NPFB * MWA_PFB_STATIONS * MWA_NPOL * MWA_NDIM);
                        /* the output offset for this particular fine channel is */
                        out_offset = out_offset + (fine_ch * NPFB * MWA_PFB_STATIONS * MWA_NPOL * MWA_NDIM);
                        /* the output offset for this pfb group is */
                        out_offset = out_offset + (pfb * MWA_PFB_STATIONS * MWA_NPOL * MWA_NDIM); 

                        current_out_ptr = out_ptr + out_offset;
                        for (station_index = 0; station_index < MWA_PFB_STATIONS ; station_index = station_index + 1) {

                            /* the fastest is station - the algorithm for this is:
                             * 	pfb*64 + floor(index/4) + index%4*16 
                             * 	but reordering into station[pol] may well be too expensive on the GPU and 
                             * 	the correlator doesn't really care ... perhaps do the re-oredering after the correlatr
                             */

                            /* But unfortunately the PFB does not pack dual pol samples next to each other
                             * so I will have to pull out the antenna in order - however .....
                             * what if we say single pol. and twice the number of antenna .....  this requires some reordering on the output
                             * but will allow us to stream straight through the unpacking */

                        	/* there appears to be a bug in the correlator that doesn;t like single pol  -
                        	 * i'l have to find that - but I could just keep the order as it is - after all
                        	 * everything has to get multiplied with everything else anyway and it is just the untangling
                        	 * that will be painfull.
                        	 */

                            /* so the lookup gives us 4x8bit numbers which is the real and imaginary samples from 2 antennas - lets just keep this the same */

                            /* we get to just increment the input pointer */
                            /* but we have to set the out_ptr */
			   
                            expand_4bit((uint16_t *) input, (int8_t *) current_out_ptr);
                            in_ptr++;
                            current_out_ptr = current_out_ptr+4;

                        }
                    }
                    pkt_num++;
                }
            }
        }
    }
}




/* Note on ops count.
 * So as can be seen this will probably take:
 *
 * NCAPTURE x NPFB * MWA_CHAN_PER_GROUP * NTIME_PER_BUFF * SAMPS_PER_PACKET * EXPAND_OPS
 *
 * 6 x 4 x 4 x 500 * 512 = 24E6 * EXPAND_OPS
 *
 * Assuming 10 OPS per EXPAND - implies 24E7 operations per buffer - each buffer is 0.05 seconds - so we need 20 per second
 *
 * therefore the re-order/expand bit operation will take 5GOPS. Note this assumes each machine only deals with 1/24th of the 
 * frequencies.
 *
 * The lookup is probably a factor of 5 faster than a straight expand and we get 4 for the price of one which gives us a 20x speed up.
 * probably bringing the unpack in under the limit for a single core.
 *
 * 
 *
 */


