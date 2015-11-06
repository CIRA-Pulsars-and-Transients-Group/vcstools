#ifndef _CUDA_ERROR_CHECK_H_
#define _CUDA_ERROR_CHECK_H_

#define CHECK_CUDA_ERROR(call)	{ gpuAssert((call), __FILE__, __LINE__); }

inline void gpuAssert( cudaError_t code, char * file, int line, bool abort=true)
{
	if( code != cudaSuccess)
	{
		fprintf( stderr, "GPUassert: %s %s %d\n", cudaGetErrorString( code ), file, line);
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


