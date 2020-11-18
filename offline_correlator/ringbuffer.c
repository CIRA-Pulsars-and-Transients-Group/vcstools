#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>

#include <unistd.h>

#include <pthread.h>

#include "ringbuffer.h"

void ringbuffer_pthread_lock(ringbuf_t *buf) {
	if (pthread_mutex_lock(&bufs->count_lock)) {
		perror("get_latest_buffer:On obtaining lock\n");
	}
}

void ringbuffer_pthread_unlock(ringbuf_t *buf) {
	if (pthread_mutex_unlock(&bufs->count_lock)) {
		perror("get_latest_buffer:On releasing lock\n");
	}
}

//apply lock on the ring buffer to determine if
//at the end of the day (assume is necessary since filled_count and EOD are volatile)
int buffer_EOD(ringbuf_t *bufs) {
	int eod = 0;
	ringbuf_pthread_lock(bufs);
	if ((bufs->filled_count == 0) && (bufs->EOD == 1)) {
		eod = 1;
	}
	ringbuf_pthread_unlock(bufs);
	return eod;
}

//apply lock and return EOD
int get_EOD(ringbuf_t *bufs) {
	int eod;
	ringbuf_pthread_lock(bufs);
	eod = bufs->EOD;
	ringbuf_pthread_unlock(bufs);
	return eod;
}

//apply lock and return EOD
int get_overrun(ringbuf_t *bufs) {
	ringbuf_pthread_lock(bufs);
	auto overrun = bufs->overrun;
	ringbuf_pthread_unlock(bufs);
	return overrun;
}


//apply lock and then get status of ring buffer
void get_buffer_status(ringbuf_t *bufs) {
	ringbuf_pthread_lock(bufs);
	fprintf(stdout,"numbufs %lu, current buf %d, fill count %d EOD %d Overrun %d\n",
		bufs->numbufs, bufs->current_buf, bufs->filled_count, bufs->EOD, bufs->overrun);
	fflush(stdout);
	ringbuf_pthread_unlock(bufs);
}

//init the ring buffer and assign buffers.
//return 1 for success, -1, for failure if buffer unfilled
int assign_ring_buffers(size_t numbufs, size_t bufsize, char **buffers, ringbuf_t * bufs) {
	unsigned int i = 0;
	//initialize the count_lock
	pthread_mutex_init(&bufs->count_lock,NULL);
	bufs->numbufs=numbufs;
	bufs->bufsize=bufsize;
	bufs->filled_count = 0;
	bufs->buffers = (char **) malloc (numbufs*sizeof(char *));

	for (i = 0; i<numbufs; i++) {
		bufs->buffers[i] = buffers[i];
		// if buffer is NULL (ie empty, halt and return failure)
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

//setup the ring buffer, and also allocate memory to the char** buffer of ring buf
//return 1 for success, -1, for failure if unable to allocate
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
	// the above statement makes no sense to me. Memory needs to be on the CPU side
	for (i = 0; i<numbufs; i++) {
		bufs->buffers[i] = NULL;
		posix_memalign ((void **)&bufs->buffers[i],(size_t) sysconf(_SC_PAGESIZE),bufsize);
		if (bufs->buffers[i] == NULL) {
			return -1;
		}
	    // set the memory to be in sequential order so pages in
    	// the given range can be aggressively read ahead, and may be
        // freed soon after they are accessed.
		madvise(bufs->buffers[i],bufsize,MADV_SEQUENTIAL);
	}
	/* the book - keeping */
	bufs->current_buf = 0;
	bufs->EOD = 0;
	bufs->overrun = 0;
	return 1;
}

// this rests the buffer to a empty state
// though why there is no lock applied is confusing
// this is not seemily called in critical sections
// and is altering volitile memory so?
int reset_ring_buffers( ringbuf_t *bufs) {

	bufs->filled_count = 0;
	/* the book - keeping */
	bufs->current_buf = 0;
	bufs->EOD = 0;
	bufs->overrun = 0;
	return 1;
}

//gets a buffer to fill based on the current buffer index
//why does index have a default value? Shouldn't it be uninitialized
char * get_buffer_to_fill(ringbuf_t *bufs) {
	int index=0;
	ringbuf_pthread_lock(bufs);
	index = bufs->current_buf;
	ringbuf_pthread_unlock(bufs);
	return bufs->buffers[index];
}

// the same structure as get buffer - but returns NULL if there is
// nothing available.
char * try_buffer_to_fill_sync(ringbuf_t *bufs) {
	int index=0;
	char *buffer = NULL;
	ringbuf_pthread_lock(bufs);
	if ((unsigned int)bufs->filled_count < bufs->numbufs)	{
		index = bufs->current_buf;
		buffer = bufs->buffers[index];
	}
	else {
		buffer = NULL;
	}
	ringbuf_pthread_unlock(bufs);
	return buffer;
}

//get the index of a buffer that can be filled, which occurs when a buffer
//has been emptied and filled_count < numbufs
//if not possible sleep for some time (why 100?) and try again
char * get_buffer_to_fill_sync(ringbuf_t *bufs) {

	int index=0;
	char *buffer = NULL;

	while (buffer == NULL) {
		ringbuf_pthread_lock(bufs);
		if ((unsigned int)bufs->filled_count < bufs->numbufs)	{
			index = bufs->current_buf;
			buffer = bufs->buffers[index];

		}
		ringbuf_pthread_unlock(bufs);
		usleep(100);
	}
	return buffer;
}

/* simple free */
// why does this return a value at all? it is always 1
int free_buffers (ringbuf_t * bufs){
	unsigned int i = 0;
	for (i=0;i<bufs->numbufs;i++) {
		free(bufs->buffers[i]);
	}
	free (bufs->buffers);
	return 1;
}


/* wait for a full buffer */
char * wait_for_buffer (ringbuf_t * bufs)
{
	char *buf = NULL;
	int index, overrun;
	/* there is a worry here as there are any number of possible race conditions
	   so the best thing to do is to have a mutex lock around the filled count
	*/
	/* if I want to check the status of the filled buffer I should aquire the lock */
	// worried about a clash accessing the overrun flag .... sloppy coding steve

	//PJE: yes, just apply a lock, unlock as overrun is volitile
	//seems this code weird locks and unlocks
	ringbuf_pthread_lock(bufs);
	overrun = bufs->overrun;
	ringbuf_pthread_unlock(bufs);
	while (buf == NULL && !overrun) {
		ringbuf_pthread_lock(bufs);
		// if the number of buffers filled is non zero, get the index to
		// empty buffer.
		if (bufs->filled_count > 0 ) {
			// see if buffer is overrun
			if ((unsigned int)bufs->filled_count <= bufs->numbufs) {
				index = bufs->current_buf - bufs->filled_count;
				// wrap index around if current buffer is passed filled count
				if (index < 0) {
					index = bufs->numbufs+index;
				}

				buf = bufs->buffers[index];
			}
			//if buffer has overrun and there are more filled buffers than allowed
			// return null and flag as overrun
			else {
				buf = NULL;
				bufs->overrun = 1;
			}
		}
		else if (bufs->filled_count == 0 && bufs->EOD == 1) {
			// PJE: the old code was unlocking twice? WHY????
			/*
			if (pthread_mutex_unlock(&bufs->count_lock)) {
				perror("wait_for_buffer:On releasing lock around filled_count\n");
			}
			buf = NULL;
			if (pthread_mutex_unlock(&bufs->count_lock)) {
						perror("wait_for_buffer:On releasing lock around filled_count\n");
			}
			*/
			// exit the loop as the buffer is empty and it is the end of the day
			buf = NULL;
			ringbuf_pthread_lock(bufs);
			break;
		}
		ringbuf_pthread_lock(bufs);
		// PJE: why sleep for 100? Is 100 some special number?
		usleep(100); // stop the spinlock ...
	}
	return buf;
}

// get a buffer further ahead in the queue -
// PJE: Old comment was "I mimic this by pretending"
/// no idea what is meant by it
char * get_advance_buffer (ringbuf_t * bufs, int nadvance)
{
	char *buf = NULL;
	int index = 0, overrun;

	/// PJE: Note before that this region had similar comments to above about
	/// data race issues which are valid

	ringbuf_pthread_lock(bufs);
	overrun = bufs->overrun;
	ringbuf_pthread_unlock(bufs);
	while (buf == NULL && !overrun) {
		ringbuf_pthread_lock(bufs);
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
			//note that here as well there were two unlocks
			buf = NULL;
			ringbuf_pthread_unlock(bufs);
			break;
		}
		ringbuf_pthread_unlock(bufs);
	}
	return buf;
}

// get the lates buffer regardless of overrun - useful if you know you aren;t going to keep up
char * get_latest_buffer (ringbuf_t * bufs) {

	char *buf = NULL;
	int index = 0;

	while (buf == NULL) {
		ringbuf_pthread_lock(bufs);
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
		ringbuf_pthread_unlock(bufs);
	}
	return buf;

}

void mark_buffer_filled (ringbuf_t * bufs) {
	ringbuf_pthread_lock(bufs);
	/* increment the filled count */
	bufs->filled_count++;
	/* increment the current buffer -- should I do this here ..... */
	///PJE: the above comment is confusing as I cannot answer this question
	bufs->current_buf++;
	if ((unsigned int)bufs->current_buf == bufs->numbufs) {
		bufs->current_buf=0;
	}
	ringbuf_pthread_unlock(bufs);
}

void mark_buffer_empty (ringbuf_t * bufs) {
	ringbuf_pthread_lock(bufs);
	/* decrement the filled count */
	bufs->filled_count--;
	ringbuf_pthread_unlock(bufs);
}
