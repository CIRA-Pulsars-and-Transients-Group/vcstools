#ifndef __BUFFER_SIZES_H
#define __BUFFER_SIZES_H

#define PFB_CAPTURE 1

/* Ideally these should be integral numbers of packets - in fact I may force this.
   The code used to be smart enough to resize these. But I think for simplicity I may force these
 */

/* A receiver packet is 168 bytes */

/* we will hold 128000 packets in a buffer as this in 0.1 seconds worth of receiver data */

/* ten buffers are therefore 1 second of data */
#if RCVR_CAPTURE
#define PACKET 168 // for receiver packet
#define NPACKET 128000 // for 0.1 seconds
#endif
#if PFB_CAPTURE
#define PACKET 264 // for PFB packet
#define NPACKET 48000 // a 500 time steps for all tiles;channels on the PFB lanes
#define HEADER_WORDS 3
#endif
#define VCS_PACKET 4096

#define GPU_PACKET (4*32*2*2)

/* A pfb packet is 262 bytes - we need to send these in 2000 packet blocks as these contain 500 time samples 
   from 16 channels - which is the atomic size for this data set */

/* The entire channel set for an individual data line is 48000 packets which contains 500 time samples for all
   the fine channels allotted to that fiber.
 */

/* A GPU packet is a convenience number which is simply the size in bytes of the data for a single packet by the time
 * it is unpacked to 8bits - these packets are never seen being simply a device for bookkeeping
 */

#define NPFB 4			// How many PFB are being captured
#define NCAPTURE 8		// How many PFB lanes per PFB are being captured
#define MAX_NCAPTURE 8 	// total number of lanes per PFB
#define PFB_LANES 2 	// Number of PFB lanes per EDT card
#define NCOARSE 24		// number of coarse PFB channels
#define NEDT 2			// number of EDT cards per demux
/* EDT buffer sizes */

#define EDT_NBUFS 4
#define EDT_BUFSZ (10*PAGE_SIZE*PACKET)

/* IPC for the shared memory */

#define IPC_NBUFS 24
#define IPC_BUFSZ (4*NPACKET*PACKET)
#define IPC_UTIME 200000 // time is useconds in 1 IPC_BUFFER

/* RINGBUFFER for the sender */
#define ATOM 500*PACKET // 500 time steps for 4 narrow channels
#define UNIT 4*ATOM // 500 time steps for 16 contig. channels
#define NATOMS 96
#define NSAMPS_PER_ATOM 500
#define TARGET_NBUFS 32 
#define TARGET_BUFSZ UNIT /* now just have 2000 packets on a line at a time */

/* RINGBUFFER for the receiver */

/* We now have to buffer at least 1 seconds worth of data. ).5 seconds worth of data is being processed and the next 0.5 seconds is being loaded
 * perhaps another 0.5 seconds of data for security
 *
 * This is broken up into frequencies. Each PACKET contains a single timesample from 4 groups of 4 channels. But 500 timesamples are present before the channel group increments.
 *
 * So the very smallest atomic size in time is a 500 packets (500*10us) - 50ms
 * TARGET_BUFSZ contains 50ms worth of data for 16 (conseq.) channels - or 2000 packets. Each data lane from the PFB has a allocation of 384 channels - therefore there are 96 sets of 500 packets covering the whole
 * allocation for that PFB lane.
 *
 * So for 32T - in order to process the whole bandpass we need to  48000 packets from each of 8 data lanes. Which in size is
 *
 * ATOM * 96 * NCAPTURE == 106.25MB
 *
 * But this only represents 500 time samples - (500 * E-4) = 0.05s (datarate is 2.120 GB/s)
 *
 * In order to get 0.5 s we need to add together 10 of these ~ 1GB
 *
 * And in order to have ~ 4 of them for "safety" we are talking 4GB
 *
 * But only ~1GB of it has to be the pinned memory - we should be fine.
 *
 * This means that the gpu library should be compiled with NTIME = 5000 NTIME_PIPE = 500.
 *
 *
 * I've increased NTIME_PER_PROCESSING BLOCK to 20. As this is the natural size of the pfb output
 *
 * But means that each buffer is 2.12GB of raw data and 4.2GB of "expanded" data. - It could be even larger - but the important thing is
 * that it needs to be pinned .... which could create a problem with the GPU.
 *
 * Note that the gpu only has to load in a 20th of this at a time - but the whole lot needs to be pinned.
 *
 */
#define NGPU 2
#define NTIME_PER_PROCESSING_BLOCK 20
#define NSAMPS_PER_RINGBUF (NSAMPS_PER_ATOM * NATOMS * NCAPTURE * NTIME_PER_PROCESSING_BLOCK)

#define RING_NBUFS 2

// Ring buffer size on receiver now calculated at runtime

// 32T MODE
// #define RING_BUFSZ (3072L*64L*2L*10000L) // NCHAN x NINPUT x NDIM x NTIME - now contains unpacked data
// 64T MODE
// #define RING_BUFSZ (256L*128L*2L*10000L) // NCHAN x NINPUT x NDIM x NTIME - now contains unpacked data
// 128T MODE
#define RING_BUFSZ (128L*256L*2L*10000L) // NCHAN x NINPUT x NDIM x NTIME - now contains unpacked data
/* GPU SIZES */
#define STAGING (2000*NPFB*MAX_NCAPTURE*GPU_PACKET)
#define EDGE 20 // how many channels to drop at the edge
/* UTILITY SIZES
 *
 */
#define TGT_PER_IPC 24 // Number of 2000 packet blocks per IPC_BUFFER (NPACKETS/2000)
#define ATOMIC_SIZE 2000*PACKET

/* for laptop testing */
/*
#define NBUFS 8
#define BUFSZ 4194304 
 */
#endif


