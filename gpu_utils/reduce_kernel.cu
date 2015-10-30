#include "reduce.h"

	/*	reduce4b reduces each block of 128 4 bit sample pairs
	 *
	 *	block width of 128 corresponds to 1 channel per block which is convenient and gives full occupancy
	 *		1 block corresponds to 1 channel within a time step
	 * 
	 *	takes in Re(u), Im(u), Re(v), Im(v) as 2's complement signed 4 bit ints concatenated together
	 *	outputs u and v as cuFloatComplex, identical in memory to four floats
	 *
	 *	relatively simple reduction - could be a lot faster
	 */
__global__ void reduce4b( sample_pair_4b_t * in_data, sample_pair_t * out )
{
		// declare variables
	__shared__ sample_pair_t  s[128];
	sample_pair_4b_t sample;
	int32_t value;
		
		// read data into shared memory
	unsigned int tid = threadIdx.x;
	unsigned int ii = threadIdx.x + blockIdx.x * blockDim.x;

	sample = in_data[ii];

	// unrolled -> no loop logic
	// last sample, lowest order 4 bits
	value = sample & 0x000F;
	if( value & 0x0008 ) // if the sign bit of the 4 bit number is set
	{	// extend the sign
		s[tid].v.y = __int2float_rz( value | 0xFFFFFFF0 );
	}
	else
	{	// the sign's already extended
		s[tid].v.y = (float) value;
	}
	sample >>= 4;

	// second last sample, second lowest order 4 bits
	value = sample & 0x000F;
	if( value & 0x0008 ) // if the sign bit of the 4 bit number is set
	{	// extend the sign
		s[tid].v.x = __int2float_rz( value | 0xFFFFFFF0 );
	}
	else
	{	// the sign's already extended
		s[tid].v.x = (float) value;
	}
	sample >>= 4;

	// second sample, third lowest order 4 bits
	value = sample & 0x000F;
	if( value & 0x0008 ) // if the sign bit of the 4 bit number is set
	{	// extend the sign
		s[tid].u.y = __int2float_rz( value | 0xFFFFFFF0 );
	}
	else
	{	// the sign's already extended
		s[tid].u.y = __int2float_rz( value );
	}
	sample >>= 4;

	// first sample, highest order 4 bits
	value = sample & 0x000F;
	if( value & 0x0008 ) // if the sign bit of the 4 bit number is set
	{	// extend the sign
		s[tid].u.x = __int2float_rz( value | 0xFFFFFFF0 );
	}
	else
	{	// the sign's already extended
		s[tid].u.x = __int2float_rz( value );
	}

	__syncthreads();

		// do the reduction
	unsigned int jj;
	for( jj = blockDim.x / 2; jj > 0; jj >>= 1 )
	{
		if( tid < jj )
		{
			s[tid].u.x += s[tid + jj].u.x;
			s[tid].u.y += s[tid + jj].u.y;
			s[tid].v.x += s[tid + jj].v.x;
			s[tid].v.y += s[tid + jj].v.y;
		}

		__syncthreads();
	}

	// write out result for this block
	if( tid == 0 ) out[blockIdx.x] = s[0];
}

	// coherent version
	// this basically just has the weighting bolted on but it probably needs to be redesigned with it in mind
__global__ void reduce4bcoherent( sample_pair_4b_t * in_data, sample_pair_t * weights, sample_pair_t * out )
{
		// declare variables
	__shared__ sample_pair_t  s[128];
	__shared__ sample_pair_t  w[128];
	sample_pair_4b_t sample;
	int value;
		
		// read data into shared memory
	unsigned int tid = threadIdx.x;
	unsigned int ii = threadIdx.x + blockIdx.x * blockDim.x;

	w[tid] = weights[ ( ii / 10000 )  ];

	sample = in_data[ii];

	// last sample, lowest order 4 bits
	value = sample & 0x000F;
	if( value & 0x0008 ) // if the sign bit of the 4 bit number is set
	{	// extend the sign
		s[tid].v.y = __int2float_rz( value | 0xFFFFFFF0 );
	}
	else
	{	// the sign's already extended
		s[tid].v.y = __int2float_rz( value );
	}
	sample >>= 4;

	// third sample, second lowest order 4 bits
	value = sample & 0x000F;
	if( value & 0x0008 ) // if the sign bit of the 4 bit number is set
	{	// extend the sign
		s[tid].v.x = __int2float_rz( value | 0xFFFFFFF0 );
	}
	else
	{	// the sign's already extended
		s[tid].v.x = __int2float_rz( value );
	}
	sample >>= 4;

	// second sample, third lowest order 4 bits
	value = sample & 0x000F;
	if( value & 0x0008 ) // if the sign bit of the 4 bit number is set
	{	// extend the sign
		s[tid].u.y = __int2float_rz( value | 0xFFFFFFF0 );
	}
	else
	{	// the sign's already extended
		s[tid].u.y = __int2float_rz( value );
	}
	sample >>= 4;

	// first sample, highest order 4 bits
	value = sample & 0x000F;
	if( value & 0x0008 ) // if the sign bit of the 4 bit number is set
	{	// extend the sign
		s[tid].u.x = __int2float_rz( value | 0xFFFFFFF0 );
	}
	else
	{	// the sign's already extended
		s[tid].u.x = __int2float_rz( value );
	}

	__syncthreads();

		// do the multiplication
	s[tid].u = cuCmulf( s[tid].u, w[tid].u );
	s[tid].v = cuCmulf( s[tid].v, w[tid].v );

		// do the reduction
	unsigned int jj;
	for( jj = blockDim.x / 2; jj > 0; jj >>= 1 )
	{
		if( tid < jj )
		{
			s[tid].u = cuCaddf( s[tid].u, s[tid + jj].u );
			s[tid].v = cuCaddf( s[tid].v, s[tid + jj].v );
		}

		__syncthreads();
	}

	// write out result for this block
	if( tid == 0 ) out[blockIdx.x] = s[0];
}
__global__ void reduce8b( sample_8b_t * in_data, uint8_t * out )
{

	// This simple version is an 8 bit total power reduction operation 
	// It does not care about the polarisation pairs and tries to reduce
	// over all inputs (256) So it will need 256 threads

	// declare variables
	__shared__ unsigned int  s[256];
		
	// read data into shared memory
	
	unsigned int tid = threadIdx.x;
	unsigned int ii = threadIdx.x + blockIdx.x * blockDim.x;
	
	sample_8b_t sample = in_data[ii];

	s[tid] = (sample.r*sample.r + sample.i*sample.i);

	__syncthreads();

		// do the reduction
	unsigned int jj;
	for( jj = blockDim.x / 2; jj > 0; jj >>= 1 )
	{
		if( tid < jj )
		{
			s[tid] += s[tid + jj];

		}

		__syncthreads();
	}

	// write out result for this block (normalised)
	if( tid == 0 ) out[blockIdx.x] = s[0]/blockDim.x;
}


	/*	reduce reduces each block of 128 sample pairs
	 *
	 *	takes pairs of complex numbers represented as 2 cuFloatComplex
	 *	
	 */
__global__ void reduce( sample_pair_t * in_data, sample_pair_t * out ) // not used anymore
{
		// declare shared memory
	__shared__ sample_pair_t  s[128];
		
		// read data into shared memory
	unsigned int tid = threadIdx.x;
	unsigned int ii = threadIdx.x + blockIdx.x * blockDim.x;

	s[tid] = in_data[ii]; //since our block size divides our total data size by definition we don't
							// actually need bounds checking

	__syncthreads();

		// do the reduction
	unsigned int jj;
	for( jj = blockDim.x / 2; jj > 0; jj >>= 1 )
	{
		if( tid < jj )
		{
			s[tid].u.x += s[tid + jj].u.x;
			s[tid].u.y += s[tid + jj].u.y;
			s[tid].v.x += s[tid + jj].v.x;
			s[tid].v.y += s[tid + jj].v.y;
		}

		__syncthreads();
	}

	// write out result for this block
	if( tid == 0 ) out[blockIdx.x] = s[0];
}
