/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef IPFB_H
#define IPFB_H

#include <cuda_runtime.h>
#include "mycomplex.h"

struct gpu_ipfb_arrays
{
    int ntaps;
    int in_size;
    int fils_size;
    int out_size;
    float *in_real,     *in_imag;
    float *d_in_real,   *d_in_imag;
    float *fils_real,   *fils_imag;
    float *d_fils_real, *d_fils_imag;
    float *d_out;
};

void cu_invert_pfb_ord( ComplexDouble ****detected_beam, int file_no,
                        int npointing, int nsamples, int nchan, int npol,
                        struct gpu_ipfb_arrays *g, float *data_buffer_uvdif );

void cu_load_filter( ComplexDouble **fils, struct gpu_ipfb_arrays *g,
        int nchan );

void malloc_ipfb( struct gpu_ipfb_arrays *g, int ntaps, int nsamples,
        int nchan, int npol, int fil_size, int npointing );

void free_ipfb( struct gpu_ipfb_arrays *g );

#endif
