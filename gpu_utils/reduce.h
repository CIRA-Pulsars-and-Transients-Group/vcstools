#ifndef __REDUCE_H__
#define __REDUCE_H__

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <cuComplex.h>

#include "cuda_error.h"

/*	Macros	*/
// struct alignment
#if defined(__CUDACC__)	// using nvcc
	#define STRUCT_ALIGN(n) __align__(n)
#elif defined(__GNUC__)
	#define STRUCT_ALIGN(n) __attribute__((aligned(n)))
#else
	#error "No struct alignment attribute defined for this compiler"
#endif

/*	Constants	*/
#define NUM_CHANNELS	128
#define NUM_STATIONS	128
#define NUM_PAIRS 		NUM_CHANNELS * NUM_STATIONS
#define NUM_POL			2
#define NUM_STEPS_PER_SECOND	10000

/*	Types	*/

typedef cuFloatComplex sample_t;

typedef struct STRUCT_ALIGN(2) {
	int8_t r; // real component
	int8_t i; // imaginary
} sample_8b_t;

typedef struct STRUCT_ALIGN(16) {
	cuFloatComplex u;
	cuFloatComplex v;
}	sample_pair_t;

typedef uint16_t sample_pair_4b_t;	// type for 4 bit sample pair

typedef int16_t sample_pair_8b_t;
	// the above should perhaps be a int8_t ...
enum beamformer_mode_t {
	COHERENT,
	INCOHERENT
};

/* CUDA Kernel Definitions */
__global__ void reduce( sample_pair_t * in_data, sample_pair_t * result );
__global__ void reduce4b( sample_pair_4b_t * in_data, sample_pair_t * result );
__global__ void reduce4bcoherent( sample_pair_4b_t * in_data, sample_pair_t * in_weights, sample_pair_t * result );
__global__ void reduce8b( sample_8b_t *in_data, uint8_t * result);

/* static helper functions */
__device__ __host__ static __inline__ sample_pair_t zero_sample_pair()
{
	sample_pair_t r;
	r.u = make_cuFloatComplex( 0.0, 0.0 );
	r.v = make_cuFloatComplex( 0.0, 0.0 );
	return r;
}

/*	Function definitions	*/
void print_sample_pair( FILE* file, sample_pair_t x );
#endif
