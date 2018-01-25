#ifndef FILTER_H
#define FILTER_H

//#define FILTERSIZE  1536
#define REAL_COEFFS  0
#define CPLX_COEFFS  1

typedef struct filter_t
{
    complex double *coeffs;
    int             size;
    int             ntaps;
} filter;

#include <complex.h>

void create_filter( filter *fil, int size, int ntaps );
void destroy_filter( filter *fil );

void load_filter( char *filename, int dtype, int ntaps, filter *fil );

void apply_phase_ramp( filter *in, double slope, filter *out );

void apply_mult_phase_ramps( filter *in, int N, filter outs[] );

void fir_filter_1D( filter *fil, complex double *signal, int size,
                    complex double *res );

int test_fir_filter_1D();

void run_all_tests();

#endif
