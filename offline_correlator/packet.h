#ifndef __PACKET_H
#define __PACKET_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
/* From the PFB ICD */

/*
 Each packet contains all of the sample data for 8 channels for exactly one 0.78125 μs time point.
 The data within the packet is ordered as 8 sequential 1.28 MHz frequency points, for each of the 16 signals.
 In order to efficiently pack the 5 bit fields into the 16 bit data words, the sample data for each 8 signals
 is spread across 5 consecutive words, in a pattern that “walks” through the bits (see Table 1).
 Thus there are 128 complex samples in each packet, and each second there are exactly 1280000 packets transmitted.
 The packet size is 1344 bits (i.e. 84 x 16 bit words), which includes a 5% overhead due to headers, etc.
 Note that the negative frequencies are not transmitted, as they are simply the complex conjugate of the positive
 frequencies that are already being transmitted, and could be reconstructed downstream if needed.
 */

#define NCHAN 8 		/* Number of channels */
#define NPOL 2			/* Number of polarizations */
#define NDIM 2			/* dimension of sampling (2 == COMPLEX) */
#define NTILE 8			/* Tiles per receiver */
#define STAT_ELEMENTS 128
#define HIST_BINS 32
#define NSTATION 32 	/* 32 tiles PROTOTYPE */

/* Some convenient consts all sizes in bits */

extern const size_t rvr_pkt_size; /* packet size */
extern const size_t hdr_size; /* header size -- fixed header (0x0800) denotes packet start */
extern const size_t node_siz; /* node id number -- denotes receiver (0 to 63) */
extern const size_t sig_size; /* id of cable (0,1,2) */
extern const size_t seq_size; /* sequence number 0 to 1279999 offset from a second boundary */
extern const size_t packet_data_size; /* payload size */
extern const size_t info_data_size; /* payload size */
extern const size_t chk_size; /* this looks like a longitudinal parity check */
extern const size_t num_samples; /* number of complex samples */

/* convenience data type */

typedef struct packet {
	uint16_t header;
	uint16_t info_block[2];
	uint16_t data[80];
	uint16_t checksum;
} pkt_t;

typedef struct bit_packet {
	uint16_t word[18]; /* word 0 and 1 is the header - then 8x8x2x2x1 */
} bit1_pkt_t;

typedef struct short_packet { 
	// this is equivelant to the gpu_x packet
	// we get to decide this as we are packing ntime x nchan x nant x npol x ndim
	// should we hold a single timestep x single channel x all ants
	// so 
	// 32 ant x 2 pol x 2 dim (all 8 bit) == 128
	uint8_t sample[NSTATION][2][2];
} bit8_pkt_t;

typedef struct pfb_packet {

	/* PFB PACKET HEADERS
	 *
	 * There are a numnber of counters in the PFB packet header as of 22/3/12 they are:
	 *
	 * uint16_t 0x"0800"
	 * uint16_t "0" & mgt_channel(5) & mgt_group(2) & mgt_id (4?) & "000" & sec_tick(1)
	 * uint16_t pfb_id(2) & mgt_bank(5) & mgt_frame (9)
	 *

       pfb_id      : in std_logic_vector(1 downto 0);
       mgt_frame   : in std_logic_vector(8 downto 0);
       mgt_channel : in std_logic_vector(4 downto 0);
       mgt_group   : in std_logic_vector(1 downto 0);
       mgt_bank    : in std_logic_vector(4 downto 0);


		write pfb_id -- pfb_id 4 x
		fill packet with boa9 -- pfb_enable_mgt_null 1 0 0 0 0 0 0 0 0 0 0 0
	 */

	/* FROM THE PFB ICD
	 * Data are packetized and put in a strictly defined sequence, in f1t[f2a]
	 (frequency-time-frequency-antenna) order. The most rapidly varying index, a, is
	 the antenna number (0, 16, 32, 48, 1, 17, 33, 49, 2,18, …, 15, 31, 47, 63),
	 followed by f2, the fine frequency channel (n..n+3), then time (0..499), and
	 most slowly the fine frequency channel group index (0, 4, 8, 12, 128, 132, 136,
	 140, ,…, 2944, 2948, 2952, 2956 for fibre 0; increment by 16n for fibre n). The
	 brackets indicate the portion of the stream that is contained within a single
	 packet.


	 */
	/* note on antenna order */
	/* antenna number is found as index
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


	 There are 16 channel groups
	 Note on narrow channel order

	 Each packet contains a single time sample from all the tiles of 4 narrow channel groups,
	 which is 16 narrow channels in total

	 How do I know what channel group I am on:

	 I have to keep count of the timesteps (packets)

	 there are 500 time steps (packets) before the narrow channel group increments


	 So first we have [floor (packet count / 500)]%4 * 4 + 16*line_number

	 0/500    = 0%4 * 4 + (0/2000) * 128 + 16*line_number = 0
	 500/500  = 1%4 * 4 + (500/2000) * 128 +16*line_number = 4
	 1000/500 = 2%4 * 4 + (1000/2000) * 128 + 16*line_number= 8
	 1500/500 = 3%4 * 4 = 12
	 2000/500 = 4%4 * 4 +  = 0


	 This means if you can capture the packets in 2000 packet blocks from each
	 data line, then you have all the fine channels for a single 1.28MHz channel.

	 Alignment: However there is no guarantee - in fact it is certain - that the data lines
	 will not be aligned.




	 */

	uint16_t header;
	uint16_t tick;
	uint16_t tick2;

	uint16_t ftc0[32]; /* 4bit real/imag for antenna 0-63; channel 0; time 0 ; freq_grp 0 */
	uint16_t ftc1[32]; /* channel 1 */
	uint16_t ftc2[32]; /* channel 2 */
	uint16_t ftc3[32]; /* channel 3 */
	uint16_t checksum;
	/* there are 500 of each narrow channel group */

} pfb_t;
typedef struct pkt_info {


	/* these are additions from the PFB packet */

	uint16_t mgt_channel;
	uint16_t mgt_group;
	uint16_t mgt_id;
	uint16_t sec_tick;
	uint16_t pfb_id;
	uint16_t mgt_bank;
	uint16_t mgt_frame;

    /* end additions */


	int node_id;

	/* no longer receiver packet specific */

	int nchan;
	int ndim;
	int npol;
	int nant;

	int signal;
	int sequence;
	int *data; /*  8 channels 8 tiles 2 pols complex*/
	int checksum;
} pkt_info_t;

typedef struct pkt_stream {
	void *block; /* the data block (packet stream is inside this block) */

	pkt_t *pkt; /* th packet stream */
	bit1_pkt_t *b1_pkt; /* I should really do this with a Union but I dont have time to do the re-write */
	bit8_pkt_t *b8_pkt;

	uint64_t npackets; /* how many packets in this stream */
	uint64_t counter; /* useful counter */
} pkt_stream_t;

typedef struct bit_stream {
	/* Very similar to the Swinburne BitSeries */
	uint64_t nsamples; /* max number of *time* samples this stream can hold */
	/*
	 I am confusing myself with this. nsamples is the number of time samples
	 NOT the number of voltage samples
	 */

	uint64_t counter; /* how many samples are currently in the stream */
	int nbit; /* nbits per sample */
	int ndim; /* 1 for real 2 for complex */
	int npol; /* number of pols */
	int nchan; /* number of channels */
	int ntile; /* number of tiles */
	size_t hdr_size; /* sizeof packet header */
	uint8_t *data; /* For future compatibility nchan,npol,ndim,nbit */

} bit_stream_t;

/* really convenient union */
typedef union full_packet {
	pkt_t pkt;
	uint16_t word[84];
} pkt_block_t;

typedef struct packet_stats {
	float mean[STAT_ELEMENTS][NTILE][NPOL][NDIM];
	float var[STAT_ELEMENTS][NTILE][NPOL][NDIM];
	int hists[STAT_ELEMENTS][NTILE][NPOL][NDIM][HIST_BINS];
	int index;
} pkt_stats_t;

typedef struct pfb_packet_stats {
	float mean[STAT_ELEMENTS][64][2];
	float var[STAT_ELEMENTS][64][2];
	int hists[STAT_ELEMENTS][64][2][16];
	int index;
} pfb_stats_t;


/* a couple of functions */

/* some derived packets - notably the 1bit packet structure */

/* it has the same 2 words as the full 5bit packet - and the data is structured the same way:
 * 8 channels each containing the 16 complex 1 bit samples for a single timestep.
 */


typedef union full_pfb_packet {
	pfb_t pkt;
	uint16_t word[132]; // including checksum
} pfb_block_t;

typedef struct gpu_sample {
	uint8_t re;
	uint8_t im;

} gpu_samp_t;

typedef struct gpu_antenna {
	gpu_samp_t pol0;
	gpu_samp_t pol1;

} gpu_ant_t;


typedef struct gpu_x_packet {
	/* this packet has to hold a single timestep for all antennas and a single channel */
	/* initially it only holds 8 bit integers - but it only has to do that on the host side */
	/* this datatype is equivelant to the bit8_pkt_t */

	gpu_ant_t station[NSTATION];

} gpu_pkt_t;



#ifdef __cplusplus
extern "C" {
#endif

void pkt_info_init(pkt_info_t *pkt);
void pkt_info_init_n(int n, pkt_info_t *pkt);
void pkt_info_destroy(pkt_info_t *pkt);
void bit_stream_init(bit_stream_t *stream);
void bit_stream_destroy(bit_stream_t *stream);

void pack(pkt_t *ptr, pkt_info_t *info); /* this takes the info and packs it */
void unpack(pkt_info_t *info, pkt_t *ptr, int npackets); /* takes packets and decodes them into info structs */
void unpack_1bit(pkt_info_t *info, bit1_pkt_t *ptr, int npackets);
void unpack_8bit(pkt_info_t *info, bit8_pkt_t *ptr, int npackets);
void unpack_1bit_hdr(pkt_info_t *info, bit1_pkt_t *ptr);
void unpack_hdr(pkt_info_t *info, pkt_t *ptr); /* as above - but just unpack the header */
void unpack_pfb_hdr(pkt_info_t *info, pfb_t *ptr); /* just unpack the pfb header to get the tick values */
void print_packet(pkt_info_t *info); /* prints some packet info to standard out */
int verify(pkt_info_t *current, pkt_info_t *previous); /* verifies checksum and contigous info packets */
pkt_t generate(int signal, int node, int sequence); /* generate a packet */
pfb_t generate_pfb(int tick, int number, int pfb_line, int mode); /* generate a pfb packet */
size_t align(pkt_stream_t *); /* align on a packet and second boundary */
int align_file(FILE **, size_t pkt_size); /* align on a packet <NO TIME ALIGNMENT> */
int align_stream(pkt_stream_t *); /* align on a packet <NO TIME ALIGNMENT> */

void get_mean(int chan_to_use, int nbit, pkt_stream_t *pkts,
		pkt_stats_t *stato);
int init_PGplot_stats();
int PGplot_stats_ticker(int nbit, pkt_stats_t *stato);

int align_1bit_file(FILE **fp);
uint16_t checksum(uint16_t *data, size_t nwords, size_t offset);

#ifdef __cplusplus
}
#endif

#endif
