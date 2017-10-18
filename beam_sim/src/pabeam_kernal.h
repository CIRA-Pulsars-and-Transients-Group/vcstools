#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuComplex.h>

#define PI (acos(-1.0)) // Ensures PI is defined on all systems

__global__ void calcArrayFactor(int nel, int ntiles, double a, double *za, double *az, float *xp, float *yp, float *zp, double tkx, double tky, double tkz, cuDoubleComplex *af);
