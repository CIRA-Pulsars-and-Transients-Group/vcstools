#ifndef _CUDA_ERROR_CHECK_H_
#define _CUDA_ERROR_CHECK_H_

#include <stdio.h>
#include <string.h>

#define CHECK_CUDA_ERROR(call)	{ gpuAssert( (cudaError_t) (call), __FILE__, __LINE__); }

inline void gpuAssert( const cudaError_t code, const char * file, const int line, bool abort=true)
{
        
        char *error_str = strdup(cudaGetErrorString( code )); 
	if( code != cudaSuccess)
	{
		fprintf( stderr, "GPUassert: %s %s %d\n", error_str, file, line);
		if( abort )
			exit( code );
	}
}
#define getLastCudaError(msg)      __getLastCudaError (msg, __FILE__, __LINE__)

inline void __getLastCudaError(const char *errorMessage, const char *file, const int line)
{
    cudaError_t err = cudaGetLastError();

    if (cudaSuccess != err)
    {
        fprintf(stderr, "%s(%i) : getLastCudaError() CUDA error : %s : (%d) %s.\n",
                file, line, errorMessage, (int)err, cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}
#endif


