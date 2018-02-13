#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <cuda_runtime.h>

extern "C" {
#include "beam_common.h"
#include "form_beam.h"
#include "mycomplex.h"
}

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    /* Wrapper function for GPU/CUDA error handling. Every CUDA call goes through
       this function. It will return a message giving your the error string,
       file name and line of the error. Aborts on error. */

    if (code != 0)
    {
        fprintf(stderr, "GPUAssert:: %s - %s (%d)\n", cudaGetErrorString(code), file, line);
        if (abort)
        {
            exit(code);
        }
    }
}

// define a macro for accessing gpuAssert
#define gpuErrchk(ans) {gpuAssert((ans), __FILE__, __LINE__);}


__global__ void beamform_kernel( uint8_t *data,
                                 ComplexDouble *W,
                                 ComplexDouble *J,
                                 int nstation,
                                 int invw,
                                 ComplexDouble *Bd,
                                 float *C,
                                 float *I )
/* Layout for input arrays:
 *   data [nsamples] [nchan] [NPFB] [NREC] [NINC] -- see docs
 *   W    [nstation] [nchan] [npol]               -- weights array
 *   J    [nstation] [nchan] [npol] [npol]        -- jones matrix
 * Layout for output arrays:
 *   Bd   [nsamples] [nchan]   [npol]             -- detected beam
 *   C    [nsamples] [nstokes] [nchan]            -- coherent full stokes
 *   I    [nsamples] [nchan]                      -- incoherent
 */
{
    // Translate GPU block/thread numbers into meaningful names
    int sample = blockIdx.x;
    int nchan  = blockDim.x;
    int ch     = threadIdx.x;

    // Assume something about the number of polarisations
    const int nstokes = 4;
    const int npol    = 2;

    // Calculate the indices for the input arrays
    int Di[nstation][npol];
    int Wi[nstation][npol];
    int Ji[nstation][npol][npol];

    int ant, pol, pol2, st;
    int pfb, rec, inc;
    for (ant = 0; ant < nstation; ant++)
    {
        pfb = ant / 32;
        inc = (ant / 8) % 4;
        for (pol = 0; pol < npol; pol++)
        {
            rec = (2*ant+pol) % 16;

            Di[ant][pol] = sample * (NINC*NREC*NPFB*nchan) +
                           ch     * (NINC*NREC*NPFB)       +
                           pfb    * (NINC*NREC)            +
                           rec    * (NINC)                 +
                           inc;

            Wi[ant][pol] = ant * (npol*nchan) +
                           ch  * (npol)       +
                           pol;

            for (pol2 = 0; pol2 < npol; pol2++)
            {
                Ji[ant][pol][pol2] = ant  * (npol*npol*nchan) +
                                     ch   * (npol*npol)       +
                                     pol  * (npol)            +
                                     pol2;
            }
        }
    }

    // Calculate the indices for the output arrays
    int Bdi[npol];
    int Ci[nstokes];
    int Ii;

    for (pol = 0; pol < npol; pol++)
        Bdi[pol] = sample * (npol*nchan) +
                   ch     * (npol)       +
                   pol;

    for (st = 0; st < nstokes; st++)
        Ci[st] = sample * (nchan*nstokes) +
                 st     * (nchan)         +
                 ch;

    Ii = sample*nchan + ch;

    // Calculate the beam and the noise floor
    ComplexDouble B[npol];
    ComplexDouble D[npol];
    ComplexDouble WD[npol];
    ComplexDouble N[npol][npol];

    for (pol = 0; pol < npol; pol++)
    {
        // Initialise beams and noise floor to zero
        Bd[Bdi[pol]] = CMaked( 0.0, 0.0 );
        I[Ii]        = 0.0;
        for (pol2 = 0; pol2 < npol; pol2++)
            N[pol][pol2] = CMaked( 0.0, 0.0 );

        for (ant = 0; ant < nstation; ant++)
        {
            // Calculate the coherent beam (B = J*W*D)
            B[pol]  = CMaked( 0.0, 0.0 );
            D[pol]  = UCMPLX4_TO_CMPLX_FLT(data[Di[ant][pol]]);
            WD[pol] = CMuld( W[Wi[ant][pol]], D[pol] );

            // (... and along the way, calculate the incoherent beam...)
            I[Ii] += DETECT(D[pol]);

            for (pol2 = 0; pol2 < npol; pol2++)
            {
                B[pol] = CAddd( B[pol], CMuld( J[Ji[ant][pol][pol2]],
                                               WD[pol2] ) );
            }

            // Detect the coherent beam
            Bd[Bdi[pol]] += CAddd( Bd[Bdi[pol]], B[pol] );

            // Calculate the noise floor (N = B*B')
            for (pol2 = 0; pol2 < npol; pol2++)
            {
                N[pol][pol2] = CAddd( N[pol][pol2],
                                      CMuld( B[pol], CConjd[B[pol2]] ) );
            }
        }
    }

    // Form the stokes parameters for the coherent beam
    float bnXX = DETECT(Bd[Bdi[0]]) - CReald(N[0][0]);
    float bnYY = DETECT(Bd[Bdi[1]]) - CReald(N[1][1]);
    ComplexDouble bnXY = CSubd( CMuld( Bd[Bdi[0]], CConj( Bd[Bdi[1]] ) ),
                                N[0][1] );

    // Stokes I, Q, U, V:
    C[Ci[0]] = invw*(bnXX + bnYY);
    C[Ci[1]] = invw*(bnXX - bnYY);
    C[Ci[2]] =  2.0*invw*CReald( bnXY );
    C[Ci[3]] = -2.0*invw*CImagd( bnXY );
}

void cu_form_beam( uint8_t *data, struct make_beam_opts *opts, ComplexDouble ***W,
                   ComplexDouble ****J, int file_no, int nstation, int nchan,
                   int npol, int outpol_coh, int outpol_incoh, int invw,
                   ComplexDouble ***detected_beam, float *coh, float *incoh )
/* The CPU version of the beamforming operations, using OpenMP for
 * parallelisation.
 *
 * Inputs:
 *   data    = array of 4bit+4bit complex numbers. For data order, refer to the
 *             documentation.
 *   opts    = passed option parameters, containing meta information about the
 *             obs and the data
 *   W       = complex weights array. [nstation][nchan][npol]
 *   J       = inverse Jones matrix. [nstation][nchan][npol][npol]
 *   file_no = number of file we are processing, starting at 0.
 *   nstation     = 128
 *   nchan        = 128
 *   npol         = 2 (X,Y)
 *   outpol_coh   = 4 (I,Q,U,V)
 *   outpol_incoh = 1 (I)
 *   invw         = the reciprocal of the sum of the antenna weights
 *
 * Outputs:
 *   detected_beam = result of beamforming operation, summed over antennas
 *                   [2*nsamples][nchan][npol]
 *   coh           = result in Stokes parameters (minus noise floor)
 *                   [nsamples][nstokes][nchan]
 *   incoh         = result (just Stokes I)
 *                   [nsamples][nchan]
 *
 * Assumes "coh" and "incoh" contain only zeros.
 */
{
    int coherent_requested = opts->out_coh    ||
                             opts->out_vdif   ||
                             opts->out_uvdif;

    int data_size = opts->sample_rate * nchan * npol * sizeof(uint8_t);
    int W_size    = nstation * nchan * npol * sizeof(ComplexDouble);
    int J_size    = nstation * nchan * npol * npol * sizeof(ComplexDouble);

    // UP TO HERE

    // Because detected_beam is large enough to contain 2 seconds' worth
    // of data, we need an index that keeps track of which "second" we're
    // in, effectively maintaining a circular buffer filled with the last
    // two seconds.
    int db_sample = (file_no % 2)*opts->sample_rate + sample;

    ComplexDouble  beam[nchan][nstation][npol];
    float          incoh_beam[nchan][nstation][npol];
    float          detected_incoh_beam[nchan*outpol_incoh];
    float          spectrum[nchan*outpol_coh];
    ComplexDouble  noise_floor[nchan][npol][npol];
    ComplexDouble  e_true[npol], e_dash[npol];

    // Initialise beam arrays to zero
    if (opts->out_coh)
    {
        // Initialise noise floor to zero
        for (ch    = 0; ch    < nchan; ch++   )
            for (opol1 = 0; opol1 < npol;  opol1++)
                for (opol2 = 0; opol2 < npol;  opol2++)
                    noise_floor[ch][opol1][opol2] = CMaked( 0.0, 0.0 );
    }

    if (coherent_requested)
    {
        // Initialise detected beam to zero
        for (ch  = 0; ch  < nchan; ch++ )
            for (pol = 0; pol < npol ; pol++)
                detected_beam[db_sample][ch][pol] = CMaked( 0.0, 0.0 );
    }

    if (opts->out_incoh)
        for (ch  = 0; ch  < nchan   ; ch++ )
            detected_incoh_beam[ch] = 0.0;

    // Calculate the beam, noise floor
    for (ant = 0; ant < nstation; ant++) {

        // Get the index for the data that corresponds to this
        //   sample, channel, antenna, polarisation
        // Justification for the (bizarre) mapping is found in the docs
        pfb = ant / 32;
        inc = (ant / 8) % 4;
        // rec depends on pol, so is calculated in the inner loop

        for (ch = 0; ch < nchan; ch++ ) {

            // Calculate quantities that depend only on "input" pol
            for (pol = 0; pol < npol; pol++) {

                rec = (2*ant+pol) % 16;

                data_idx = sample * (NINC*NREC*NPFB*nchan) +
                    ch     * (NINC*NREC*NPFB)       +
                    pfb    * (NINC*NREC)            +
                    rec    * (NINC)                 +
                    inc;

                // Form a single complex number
                e_dash[pol] = UCMPLX4_TO_CMPLX_FLT(data[data_idx]);

                // Detect the incoherent beam, if requested
                if (opts->out_incoh)
                    incoh_beam[ch][ant][pol] = DETECT(e_dash[pol]);

                // Apply complex weights
                if (coherent_requested)
                    e_dash[pol] = CMuld( e_dash[pol], W[ant][ch][pol] );

            }

            // Calculate quantities that depend on output polarisation
            // (i.e. apply inv(jones))
            if (coherent_requested)
            {
                for (pol = 0; pol < npol; pol++)
                {
                    e_true[pol] = CMaked( 0.0, 0.0 );

                    for (opol = 0; opol < npol; opol++)
                        e_true[pol] = CAddd( e_true[pol],
                                CMuld( J[ant][ch][pol][opol],
                                    e_dash[opol] ) );

                    if (opts->out_coh)
                        for (opol = 0; opol < npol; opol++)
                            noise_floor[ch][pol][opol] =
                                CAddd( noise_floor[ch][pol][opol],
                                        CMuld( e_true[pol],
                                            CConjd(e_true[opol]) ) );

                    beam[ch][ant][pol] = e_true[pol];
                }
            }
        }
    }

    // Detect the beam = sum over antennas
    for (ant = 0; ant < nstation; ant++)
        for (pol = 0; pol < npol    ; pol++)
            for (ch  = 0; ch  < nchan   ; ch++ )
            {
                // Coherent beam
                if (coherent_requested)
                    detected_beam[db_sample][ch][pol] =
                        CAddd( detected_beam[db_sample][ch][pol],
                                beam[ch][ant][pol] );

                // Incoherent beam
                if (opts->out_incoh)
                    detected_incoh_beam[ch] += incoh_beam[ch][ant][pol];
            }

    if (opts->out_coh)
    {
        // Calculate the Stokes parameters
        form_stokes( detected_beam[db_sample],
                noise_floor,
                nchan,
                invw,
                spectrum );

        int offset_in_coh = sizeof(float)*nchan*outpol_coh*sample;
        memcpy((void *)((char *)coh + offset_in_coh),
                spectrum,
                sizeof(float)*nchan*outpol_coh);
    }

    if (opts->out_incoh)
    {
        int offset_in_incoh = sizeof(float)*nchan*outpol_incoh*sample;
        memcpy((void *)((char *)incoh + offset_in_incoh),
                detected_incoh_beam,
                sizeof(float)*nchan*outpol_incoh);
    }

}

