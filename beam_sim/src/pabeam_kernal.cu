#include "pabeam_kernal.h"

__global__ void calcArrayFactor(int nel, int ntiles, double a, double *za, double *az, float *xp, float *yp, float *zp, double tkx, double tky, double tkz, cuDoubleComplex *af)
{
    /* Kernal which takes (az,za) coordinates and tile positions to 
       compute the array factor. This minimises data transfer between 
       host and device.
   

       ##########################
       # Wavenumber computation #
       ##########################
       a     = amplitude factor (2*pi/lambda)
       za    = array of zenith angles
       az    = array of azimuths
       tkx, tky, tkz = target wavenumber components

       NOTE: The standard equations are:
             kx = a * sin(theta) * cos(phi)
             ky = a * sin(theta) * sin(phi)
             kz = a * cos(theta)

             where:
             lambda is the observing wavelength, and
             assuming that (theta,phi) are in the convention from Sutinjo et al. 2015:
             phi = pi/2 - az   AND   theta = za

             The azimuth is measured clockwise from East (standard for antenna theory, offset from astronomy)
       
       ############################
       # Array factor computation #
       ############################
       nel    = total number of elements in az,za,af arrays
       ntiles = number of tiles used to form tied-array beam
       xp     = array of tile x-positions (East)
       yp     = array of tile y-positions (North)
       zp     = array of tile z-positions (above array centre)
       af     = array containing the complex valued array factor

       NOTE: The analytical form for this is:
             
                f(theta,phi;tza,taz) = (1/ntiles) * sum{n=1,n=ntiles}( conj(psi_n(tza,taz)) * psi_n(theta,phi) )

             where:
             (taz,tza) are the target source azimuth and zenith angle
             (theta,phi) are the azimuth and zenith angle pixels over which we evalute, theta=[0,90], phi=[0,360)
             
             psi_n is the planar wave front as detected by the "nth" tile, defined as
                
                psi_n(theta,phi) = exp[ (2*pi*I/lambda) * (x_n*k_x + y_n*k_y + z_n*k_z) ]
             
             where x_n is the x-coordinate of tile n (similarly for y_n, z_n), with k_x, k_y, k_z and lambda as defined above.      
    */

    // set up CUDA thread indexing
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    
    // other intermediate variables
    double ast=0, phi=0;
    double ph=0;
    double kx=0, ky=0, kz=0;
    cuDoubleComplex n = make_cuDoubleComplex(ntiles,0);

    
    for (int i = index; i < nel; i += stride)
    {
        // pre-calculate coefficients/transforms
        ast = a * sin(za[i]);
        phi = PI/2 - az[i];

        // calculate (k - k_target)
        kx = ast * cos(phi) - tkx; 
        ky = ast * sin(phi) - tky;
        kz = a * cos(za[i]) - tkz;

        // initialise this pixel's array factor value
        af[i] = make_cuDoubleComplex(0,0);
        
        // calculate array factor contribution from each tile and sum
        for (int j = 0; j < ntiles; j++)
        {
            ph = (kx * xp[j]) + (ky * yp[j]) + (kz * zp[j]);
            af[i] = cuCadd(af[i], make_cuDoubleComplex(cos(ph), sin(ph)));         
        }

        // normalise the array factor
        af[i] = cuCdiv(af[i], n);
    }
    __syncthreads();
}
