#ifndef _CUDA_ERROR_CHECK_H_
#define _CUDA_ERROR_CHECK_H_

#define CHECK_CUDA_ERROR(call)	{ gpuAssert( (cudaError_t) (call), __FILE__, __LINE__); }

inline void gpuAssert( cudaError_t code, char * file, int line, bool abort=true)
{
	if( code != cudaSuccess)
	{
		fprintf( stderr, "GPUassert: %s %s %d\n", cudaGetErrorString( code ), file, line);
		if( abort )
			exit( code );
	}
}

#endif
