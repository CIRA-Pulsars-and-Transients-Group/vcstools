#include <cuda_runtime.h>
#include <cufft.h>
#include <helper_cuda.h>
#include <helper_functions.h>

#include "gpu_utils.h"

static __device__ __host__ inline Complex ComplexAdd(Complex, Complex);
static __device__ __host__ inline Complex ComplexScale(Complex, float);
static __device__ __host__ inline Complex ComplexMul(Complex, Complex);
static __global__ void ComplexPointwiseMulAndScale(Complex *, const Complex *, int, float);
static __global__ void UpSample(Complex *, Complex *, int, int);
static __global__ void ElementWiseMultiply(Complex *, Complex *, int );
static __global__ void ShiftandReduce(Complex *signal, int offset, int shift, int nelements); 
void Filter1D( Complex  * input, Complex * output, Complex * filter, int nchan_in, int ntap, int nsamples )
{
    int mem_size = sizeof(Complex) * nsamples;

    // Allocate device memory for signal
    Complex *d_signal;
    checkCudaErrors(cudaMalloc((void **)&d_signal, mem_size));
    // Copy host memory to device
    checkCudaErrors(cudaMemcpy(d_signal, input, mem_size,cudaMemcpyHostToDevice));

    // Allocate device memory for filter kernel
    Complex *d_filter_kernel;
    checkCudaErrors(cudaMalloc((void **)&d_filter_kernel, mem_size));

    // Copy host memory to device
    checkCudaErrors(cudaMemcpy(d_filter_kernel, filter, mem_size,cudaMemcpyHostToDevice));

    // CUFFT plan
    cufftHandle plan;
    checkCudaErrors(cufftPlan1d(&plan, nsamples, CUFFT_C2C, 1));

    // Transform signal and kernel
    printf("Transforming signal cufftExecC2C\n");
    checkCudaErrors(cufftExecC2C(plan, (cufftComplex *)d_signal, (cufftComplex *)d_signal, CUFFT_FORWARD));
    checkCudaErrors(cufftExecC2C(plan, (cufftComplex *)d_filter_kernel, (cufftComplex *)d_filter_kernel, CUFFT_FORWARD));

    // Multiply the coefficients together and normalize the result
    printf("Launching ComplexPointwiseMulAndScale<<< >>>\n");
    ComplexPointwiseMulAndScale<<<32, 256>>>(d_signal, d_filter_kernel, nsamples, 1.0f / nsamples);

    // Check if kernel execution generated and error
    getLastCudaError("Kernel execution failed [ ComplexPointwiseMulAndScale ]");

    // Transform signal back
    printf("Transforming signal back cufftExecC2C\n");
    checkCudaErrors(cufftExecC2C(plan, (cufftComplex *)d_signal, (cufftComplex *)d_signal, CUFFT_INVERSE));

    // Copy device memory to host
    Complex *h_convolved_signal = output;
    checkCudaErrors(cudaMemcpy(h_convolved_signal, d_signal, mem_size,cudaMemcpyDeviceToHost));

    // Clean up

    checkCudaErrors(cudaFree(d_signal));
    checkCudaErrors(cudaFree(d_filter_kernel));

}
void cuda_invert_pfb(Complex  * input, Complex * output, float * filter, int nchan_in, int ntap, int nsamples) {
// A more specific version of the filter that tries to invert the PFB a bit more efficiently

    int mem_size = sizeof(Complex) * nsamples * nchan_in * nchan_in;
    int nelements = nsamples * nchan_in;

    static int planned = 0;    
    //Allocate memory for buffer
    Complex *d_buffer;
    checkCudaErrors(cudaMalloc((void **)&d_buffer, nelements*sizeof(Complex)));
    // Move host memeory to device
    checkCudaErrors(cudaMemcpy(d_buffer, input, nelements*sizeof(Complex),cudaMemcpyHostToDevice));

    // Allocate device memory for signal matrix
    Complex *d_signal;
    checkCudaErrors(cudaMalloc((void **)&d_signal, mem_size));
    getLastCudaError("Kernel execution failed [ cudaMalloc ]");
    // default to zero
    checkCudaErrors(cudaMemset((void *) d_signal,0,mem_size));

    getLastCudaError("Kernel execution failed [ cudaMemset ]");
    // Allocate device memory for filter kernel
    Complex *d_filter_kernel;
    checkCudaErrors(cudaMalloc((void **)&d_filter_kernel, nelements*sizeof(Complex)));
    
    getLastCudaError("Kernel execution failed [ cudaMalloc ]");
    // Kernel to copy host signal to upsampled device matrix
    // and simultaneously perform the transpose

    UpSample<<<nchan_in,256>>>((Complex *) d_signal, (Complex *) d_buffer, nchan_in,nsamples);
    getLastCudaError("Kernel execution failed [ UpSample ]");
    // CUFFT plan
    // Batch FFT to do all the FFTs of the rows of the upsampled device matrix
    static cufftHandle plan;
    if (!planned) {
        checkCudaErrors(cufftPlan1d(&plan, nsamples, CUFFT_C2C, nchan_in));
        planned = 1;
    }
    // Transform signal and kernel
    // printf("Transforming signal cufftExecC2C\n");
    checkCudaErrors(cufftExecC2C(plan, (cufftComplex *)d_signal, (cufftComplex *)d_signal, CUFFT_FORWARD));

    // FFT padded filter
    
    checkCudaErrors(cufftExecC2C(plan, (cufftComplex *)d_filter_kernel, (cufftComplex *)d_filter_kernel, CUFFT_FORWARD));
    // multiply each row of the matrix by the filter vector
    ElementWiseMultiply<<<nchan_in,256>>>((Complex *) d_signal, (Complex *) d_filter_kernel, (int) nelements);
    getLastCudaError("Kernel execution failed [ElementWiseMultiply]");
    // Fourier Transform Back 
   
    // printf("Transforming signal back cufftExecC2C\n");
    checkCudaErrors(cufftExecC2C(plan, (cufftComplex *)d_signal, (cufftComplex *)d_signal, CUFFT_INVERSE));
   
    // Permute each row by nsamp entries
    // and Reduce down the columns
//    for (int c = 0; c < nchan_in-1; c++){
 //       int offset = c*nelements;
 //       int shift = c*nsamples;
 //       ShiftandReduce<<<nelements,1>>>((Complex *) d_signal, (int) offset, (int) shift, (int) nelements );
 //       getLastCudaError("Kernel execution failed [ShiftandReduce]");
 //   } 
    // Copy last row to host

    int offset = (nchan_in-1)*nelements;

    // Copy device memory to host
    Complex *h_convolved_signal = output;
    Complex *d_outptr = &d_signal[offset];
    checkCudaErrors(cudaMemcpy(h_convolved_signal, d_outptr, nelements*sizeof(Complex),cudaMemcpyDeviceToHost));

    // Clean up

    checkCudaErrors(cudaFree(d_signal));
    checkCudaErrors(cudaFree(d_buffer));
    checkCudaErrors(cudaFree(d_filter_kernel));


}
static __device__ __host__ inline Complex ComplexAdd(Complex a, Complex b)
{
    Complex c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    return c;
}

// Complex scale
static __device__ __host__ inline Complex ComplexScale(Complex a, float s)
{
    Complex c;
    c.x = s * a.x;
    c.y = s * a.y;
    return c;
}

// Complex multiplication
static __device__ __host__ inline Complex ComplexMul(Complex a, Complex b)
{
    Complex c;
    c.x = a.x * b.x - a.y * b.y;
    c.y = a.x * b.y + a.y * b.x;
    return c;
}

// Complex pointwise multiplication
static __global__ void ComplexPointwiseMulAndScale(Complex *a, const Complex *b, int size, float scale)
{
    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

    for (int i = threadID; i < size; i += numThreads)
    {
        a[i] = ComplexScale(ComplexMul(a[i], b[i]), scale);
    }
}

static __global__ void UpSample(Complex *signal, Complex *input, int nchan, int nsamp) {

    const int input_chan = blockIdx.x ;
    
    int inputID = 0;
    int outputID = 0;
    for (int sampnum = threadIdx.x ; sampnum < nsamp; sampnum=sampnum+blockDim.x) {

        inputID = input_chan + sampnum*nchan;
        outputID = blockIdx.x * (nsamp * nchan) + (sampnum * nchan);
        signal[outputID] = input[inputID];
        sampnum++;
    }
}
static __global__ void ElementWiseMultiply(Complex *signal, Complex *filter, int nelements) {

    const int offset = blockIdx.x * nelements;

    for (int entry = threadIdx.x ; entry<nelements; entry=entry+gridDim.x){
        int inputID = entry + offset;
        int outputID = inputID;
        signal[outputID] = ComplexMul(signal[inputID],filter[entry]);
    }
}
static __global__ void ShiftandReduce(Complex *signal, int offset, int shift, int nelements) {

    const int element = blockDim.x ;
   
    int output_col = element+shift;
    if (output_col > nelements) {
       output_col = output_col - nelements;
    }

    int inputID = offset+element;
    int outputID = offset + output_col + nelements;

    signal[outputID] = ComplexAdd(signal[inputID],signal[outputID]);

}  
