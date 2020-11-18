#ifndef __CORR_UTILS_H
#define __CORR_UTILS_H

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<complex.h>
#include<assert.h>
#include<syslog.h>
#include<string.h>
#include<cstdint>
#include<strings.h> //what is strings.h
#include<cstring> //added this as I think strings.h is a typo
#include <complex.h>

// include the local ringbuffer implementation
#include "ringbuffer.h"
// what is the xgpu library?
#include "xgpu.h"


#ifdef __cplusplus
extern "C" {
#endif

typedef struct management {

        int data_line[8];
        int data_line_count;
        ringbuf_t *ring;
        char start_obs_id[64];
        char start_obs_UTC[256];
        int new_obs;
        int offline;
        int integrate; // how many seconds to integrate
        int chan_to_aver;     // how many channels to average
        int marker;
        int shutdown;
        int nfrequency;
        int ntime;
        int nstation;
        int ndim;
        int npol;
        int edge;
        int nbit;
        int ncoarse;
        int coarse_chan;
        int dumps_per_sec;
	int infile;

} manager_t;

void printerror(int);
inline void integrate(char *,uint8_t *, int, int,int);
void buildFITSBuffer(XGPUInfo,Complex *,size_t,void *,time_t,int, manager_t *);

#ifdef __cplusplus
}
#endif

#endif
