#ifndef __GPU_UTILS_H
#define __GPU_UTILS_H

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <cuComplex.h>

typedef float2 Complex;

/* CUDA Wrapper Definitions */
/* data will be in the in_data location and returned in the result location.
 * Assuming nsamples in both the filter and the in_data (padding already done
*/
#ifdef __cplusplus
extern "C" {
#endif

void Filter1D(Complex * in_data, Complex * result, Complex *filter, int nchan_in, int ntap, int sample );
void cuda_invert_pfb(Complex * in_data, Complex * result, float *filter, int nchan_in, int ntap, int sample );

#ifdef __cplusplus
}
#endif

#endif
/* end */



