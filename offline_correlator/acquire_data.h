#ifndef __ACQUIRE_DATA_H
#define __ACQUIRE_DATA_H

#include <sys/socket.h> /* for socket() and bind() */
#include <arpa/inet.h>  /* for sockaddr_in and inet_ntoa() */
#include <pthread.h>
#include "twopip_header.h"
#include "ringbuffer.h"
#include "buffer_sizes.h"




pthread_mutex_t acquire_from_EDT_lock;
pthread_cond_t acquire_from_EDT_var;

volatile int occupancy[RING_NBUFS][NPFB][NCAPTURE][NTIME_PER_PROCESSING_BLOCK][NCOARSE]; // this variable is intended to be protected by the above lock
volatile int occupancy_count[RING_NBUFS];

typedef struct test_acq_ctxt	{
	ringbuf_t *buf;
	int n_to_acquire;
} test_acq_t;

typedef struct socket_acq_ctxt  {
	ringbuf_t *buf;

	unsigned short socket;
	struct sockaddr_in ServAddr; /* Local address */
	struct sockaddr_in ClntAddr; /* Client address */
	unsigned int cliAddrLen; 
    char obsheader[TWOPIP_HEADER_SIZE]; 
    int id;

    int pfb_lanes;
    /* for 32T I am only capturing pfb0 but for all channels */
    /* for 128T I should specify a coarse channel on the command line and capture some number of PFB */
    /* this logic should be passed through the socket_context argument - but for 32T will be hardcoded */

    int coarse_channels;
    int capture;


    int connected;

} socket_acq_t;

#ifdef __cplusplus
extern "C" {
#endif
 

void * acquire_from_EDT(void *context);
pthread_t acquire_test_data(ringbuf_t * bufs, int numbufs); // perhaps a function pointer would be nicer here
void * acquire_data_from_socket(void *);
 
#ifdef __cplusplus
}
#endif

#endif
