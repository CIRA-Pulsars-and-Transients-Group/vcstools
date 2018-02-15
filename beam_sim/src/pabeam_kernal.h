#ifndef KERNEL_H
#define KERNEL_H

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuComplex.h>
#include "utils.h"

__global__ void calcArrayFactor(int nel, int ntiles, double a, 
                                double *za, double *az, 
                                float *xp, float *yp, float *zp, 
                                double tkx, double tky, double tkz, 
                                cuDoubleComplex *af);
#endif
