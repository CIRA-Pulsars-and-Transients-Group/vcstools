#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>

#include <unistd.h>

#include <pthread.h>

#include "ringbuffer.h"


int buffer_EOD(ringbuf_t *bufs) {

	int eod = 0;
	if (pthread_mutex_lock(&bufs->count_lock)) {
		perror("wait_for_buffer:On obtaining lock around filled_count\n");
	}

	if ((bufs->filled_count == 0) && (bufs->EOD == 1)) {
		eod = 1;
	}

	if (pthread_mutex_unlock(&bufs->count_lock)) {
			perror("get_latest_buffer:On releasing lock around filled_count\n");
	}

	return eod;


}
/* allocate memory and initialise the type */

void get_buffer_status(ringbuf_t *bufs) {
	if (pthread_mutex_lock(&bufs->count_lock)) {
		perror("get_latest_buffer:On obtaining lock arounf filled_count\n");
	}
	fprintf(stdout,"numbufs %lu, current buf %d, fill count %d EOD %d Overrun %d\n",bufs->numbufs,bufs->current_buf,bufs->filled_count,bufs->EOD,bufs->overrun);
	fflush(stdout);
	if (pthread_mutex_unlock(&bufs->count_lock)) {
		perror("get_latest_buffer:On releasing lock around filled_count\n");
	}
}

int assign_ring_buffers(size_t numbufs, size_t bufsize, char **buffers, ringbuf_t * bufs) {
	unsigned int i = 0;
	/* first set the basics */
	pthread_mutex_init(&bufs->count_lock,NULL);
	bufs->numbufs=numbufs;
	bufs->bufsize=bufsize;
	bufs->filled_count = 0;
	bufs->buffers = (char **) malloc (numbufs*sizeof(char *));

	for (i = 0; i<numbufs; i++) {
		bufs->buffers[i] = buffers[i];
		if (bufs->buffers[i] == NULL) {
			return -1;
		}
	}
	/* the book - keeping */
	bufs->current_buf = 0;
	bufs->EOD = 0;
	bufs->overrun = 0;
	return 1;

}
int configure_ring_buffers(size_t numbufs,size_t bufsize, ringbuf_t * bufs) {

	unsigned int i = 0;
	/* first set the basics */
	pthread_mutex_init(&bufs->count_lock,NULL);
	bufs->numbufs=numbufs;
	bufs->bufsize=bufsize;
	bufs->filled_count = 0;
	bufs->buffers = (char **) malloc (numbufs*sizeof(char *));
	if (bufs->buffers == NULL) {
		return -1;
	}
	/* now malloc the ring buffers this is probably best kept 
	   on the CPU side */

	for (i = 0; i<numbufs; i++) {
		bufs->buffers[i] = NULL;
		posix_memalign (&bufs->buffers[i],(size_t) sysconf(_SC_PAGESIZE),bufsize);
		if (bufs->buffers[i] == NULL) {
			return -1;	
		}
		madvise(bufs->buffers[i],bufsize,MADV_SEQUENTIAL);
	}

	/* the book - keeping */
	bufs->current_buf = 0;
	bufs->EOD = 0;
	bufs->overrun = 0;
	return 1;

}
int reset_ring_buffers( ringbuf_t *bufs) {

	bufs->filled_count = 0;
	/* the book - keeping */
	bufs->current_buf = 0;
	bufs->EOD = 0;
	bufs->overrun = 0;
	return 1;
}

char * get_buffer_to_fill(ringbuf_t *bufs) {
	int index=0;
	if (pthread_mutex_lock(&bufs->count_lock)) {
		perror("wait_for_buffer:On obtaining lock arounf filled_count\n");
	}

	index = bufs->current_buf;

	if (pthread_mutex_unlock(&bufs->count_lock)) {
		perror("wait_for_buffer:On releasing lock around filled_count\n");
	}
	return bufs->buffers[index];
}

char * try_buffer_to_fill_sync(ringbuf_t *bufs) {

	/* the same structure as get buffer - but returns NULL if there is
	 * nothing available.
	 */

	int index=0;
	char *buffer = NULL;

	if (pthread_mutex_lock(&bufs->count_lock)) {
		perror("wait_for_buffer:On obtaining lock around filled_count\n");
	}
	if ((unsigned int)bufs->filled_count < bufs->numbufs)	{
		index = bufs->current_buf;
		buffer = bufs->buffers[index];
	}
	else {
		buffer = NULL;
	}

	if (pthread_mutex_unlock(&bufs->count_lock)) {
		perror("wait_for_buffer:On releasing lock around filled_count\n");
	}

	return buffer;

}
char * get_buffer_to_fill_sync(ringbuf_t *bufs) {

	int index=0;
	char *buffer = NULL;

	while (buffer == NULL) {


		if (pthread_mutex_lock(&bufs->count_lock)) {
			perror("wait_for_buffer:On obtaining lock around filled_count\n");
		}
		if ((unsigned int)bufs->filled_count < bufs->numbufs)	{
			index = bufs->current_buf;
			buffer = bufs->buffers[index];

		}
		if (pthread_mutex_unlock(&bufs->count_lock)) {
			perror("wait_for_buffer:On releasing lock around filled_count\n");
		}

		usleep(100);

	}
	return buffer;
}

/* simple free */

int free_buffers (ringbuf_t * bufs){
	unsigned int i = 0;
	for (i=0;i<bufs->numbufs;i++) {
		free(bufs->buffers[i]);
	}
	free (bufs->buffers);
	return 1;
}
char * wait_for_buffer (ringbuf_t * bufs) 
{
	char *buf = NULL;	
	int index = 0;


	/* wait for a full buffer */
	
	/* there is a worry here as there are any number of possible race conditions
	   so the best thing to do is to have a mutex lock around the filled count
	*/	
	/* if I want to check the status of the filled buffer I should aquire the lock */
	while (buf == NULL && !bufs->overrun) { // worried about a clash accessing the overrun flag .... sloppy coding steve
		if (pthread_mutex_lock(&bufs->count_lock)) {
			perror("wait_for_buffer:On obtaining lock arounf filled_count\n");
		}
		if (bufs->filled_count > 0 ) {
			if ((unsigned int)bufs->filled_count <= bufs->numbufs) {
				index = bufs->current_buf - bufs->filled_count;
				if (index < 0) {
					index = bufs->numbufs+index;
				}
				buf= bufs->buffers[index];
			}
			else {
				buf = NULL;
				bufs->overrun = 1;
			}
		}
		else if (bufs->filled_count == 0 && bufs->EOD == 1) {
				
			if (pthread_mutex_unlock(&bufs->count_lock)) {
				perror("wait_for_buffer:On releasing lock around filled_count\n");
			}
			buf = NULL;
			if (pthread_mutex_unlock(&bufs->count_lock)) {
						perror("wait_for_buffer:On releasing lock around filled_count\n");
			}
			break;
		}
		if (pthread_mutex_unlock(&bufs->count_lock)) {
			perror("wait_for_buffer:On releasing lock around filled_count\n");
		}
		usleep(100); // stop the spinlock ...
	}
	return buf;
}
char * get_advance_buffer (ringbuf_t * bufs, int nadvance)
{
	char *buf = NULL;
	int index = 0;


	/* get a buffer further ahead in the queue -
	 * I mimic this by pretending
	 */

	/* there is a worry here as there are any number of possible race conditions
	   so the best thing to do is to have a mutex lock around the filled count
	*/
	/* if I want to check the status of the filled buffer I should aquire the lock */
	while (buf == NULL && !bufs->overrun) { // worried about a clash accessing the overrun flag .... sloppy coding steve
		if (pthread_mutex_lock(&bufs->count_lock)) {
			perror("get_advance_buffer:On obtaining lock around filled_count\n");
		}
		if (bufs->filled_count >= 0 ) { // this lets youget buffers even if there are no full buffers
			if ((unsigned int)(bufs->filled_count+nadvance) <= bufs->numbufs) {

				// we have space to give the advanced buffer - now which is it
				index = bufs->current_buf - bufs->filled_count;

				if (index < 0) {
					index = bufs->numbufs+index;
				}
				index = index+nadvance;
				index = index%bufs->numbufs;
				buf= bufs->buffers[index];
			}
			else {
				buf = NULL;
			}
		}
		else if (bufs->filled_count == 0 && bufs->EOD == 1) {

			if (pthread_mutex_unlock(&bufs->count_lock)) {
				perror("wait_for_buffer:On releasing lock around filled_count\n");
			}
			buf = NULL;
			if (pthread_mutex_unlock(&bufs->count_lock)) {
						perror("wait_for_buffer:On releasing lock around filled_count\n");
			}
			break;
		}
		if (pthread_mutex_unlock(&bufs->count_lock)) {
			perror("wait_for_buffer:On releasing lock around filled_count\n");
		}
	}
	return buf;
}
char * get_latest_buffer (ringbuf_t * bufs) {

	char *buf = NULL;	
	int index = 0;
	/* get the lates buffer regardless of overrun - useful if you know you aren;t going to keep up */
	
	while (buf == NULL) { 
		if (pthread_mutex_lock(&bufs->count_lock)) {
			perror("get_latest_buffer:On obtaining lock arounf filled_count\n");
		}
		if (bufs->filled_count > 0 ) { // we dont care if we have overrun in this case but set the flags anyway
			if ((unsigned int)bufs->filled_count >= bufs->numbufs) {
				bufs->overrun = 1;
			}
			index = bufs->current_buf -1;
			if (index < 0) {
				index = index + bufs->numbufs;
			}
			buf = bufs->buffers[index];
		}
		if (pthread_mutex_unlock(&bufs->count_lock)) {
			perror("get_latest_buffer:On releasing lock around filled_count\n");
		}
	}
	return buf;

}
void mark_buffer_filled (ringbuf_t * bufs) {
	if (pthread_mutex_lock(&bufs->count_lock)) {
		perror("mark_buffer_filled: On obtaining the lock around filled_count\n");
	}
	/* increment the filled count */
	bufs->filled_count++;
	/* increment the current buffer -- should I do this here ..... */
	bufs->current_buf++;
	if ((unsigned int)bufs->current_buf == bufs->numbufs) {
		bufs->current_buf=0;
	}
	if (pthread_mutex_unlock(&bufs->count_lock)) {
		perror("mark_buffer_filled: On releasing the lock around filled_count\n");
	}
}
void mark_buffer_empty (ringbuf_t * bufs) {
	if (pthread_mutex_lock(&bufs->count_lock)) {
		perror("mark_buffer_empty: On obtaining the lock around filled_count\n");
	}
	/* increment the filled count */
	bufs->filled_count--;
	if (pthread_mutex_unlock(&bufs->count_lock)) {
		perror("mark_buffer_empty: On releasing the lock around filled_count\n");
	}
}
