#ifndef __RINGBUF_H
#define __RINGBUF_H

/* simple ring buffer kit - i've based the logic on 2PiP's shared memory ringbuffers
   but this simple invocation can simply be run as a threaded version in the same process space
   (c) Steve Ord August 2010 (steve.ord@gmail.com)
*/

/* configure_ring_buffers - will set up the buffers and the structs
   aquire_data will fill the buffers
   wait_for_buffer will get the next buffer and also raise an EOD flag if the buffer is 
	partially filled
   get_latest_buffer will get the most recent and not nec. the next buffer
   mark_buffer_filled will increment the write count
   mark_buffer_empty will increment the read count
*/

/* 
   This logic is identical to that of 2pip. I think we have a complication in that
   there are multiple input lines?

   how are they arranged. I will assume all baselines but subset of the channels - so they will have to be combined

   if this is wrong at least it makes sense and can probably be fixed later

   abstract that. and just deal with a single line - I can combine them after the basic structure is written.

*/
#include <stdlib.h>
#include <pthread.h>

typedef struct generic_ringbuffer {

	size_t numbufs;
	size_t bufsize;

	pthread_mutex_t count_lock;	
	volatile int filled_count;
	volatile int current_buf;
	volatile int overrun;
	volatile int EOD;

	char **buffers;

} ringbuf_t;

/* so the logic here is configure ringbufs sets this struct up */
/* then aquire_data launces a thread and fills these buffers up */


#ifdef __cplusplus
extern "C" {
#endif

int assign_ring_buffers(size_t numbufs, size_t bufsize, char **buffers, ringbuf_t * bufs); // assign buffers that have been created elsewhere
int configure_ring_buffers(size_t numbufs, size_t bufsize, ringbuf_t * bufs);
char * get_buffer_to_fill(ringbuf_t *bufs); // buffer to write into
char * try_buffer_to_fill_sync(ringbuf_t *bufs); // return empty buffer to write into or NULL (non-blocking)
char * get_buffer_to_fill_sync(ringbuf_t *bufs); // block until <empty> buffer to write into
char * wait_for_buffer(ringbuf_t * bufs); // buffer to read from
char * get_latest_buffer(ringbuf_t * bufs); // latest not next buffer
char * get_advance_buffer(ringbuf_t *bufs, int nadvance); //get next and future buffers

void mark_buffer_filled(ringbuf_t * bufs); // buffer full
void mark_buffer_empty(ringbuf_t * bufs); // buffer finished with
int free_buffers(ringbuf_t * bufs); // free the memory
void get_buffer_status(ringbuf_t *);
int reset_ring_buffers(ringbuf_t *);
int buffer_EOD(ringbuf_t *);

#ifdef __cplusplus
}
#endif


#endif


