#ifndef IPFB_H
#define IPFB_H

#include <cuda_runtime.h>
#include <cuComplex.h>

void cu_invert_pfb_ord( ComplexFloat ***detected_beam, int file_no,
                        int nsamples, int nchan, int npol,
                        ComplexDouble **fils, int fil_size,
                        float *data_buffer_uvdif );

#endif
