/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "beam_common.h"
#include "form_beam.h"
#include "mycomplex.h"
#include <omp.h>

void form_beam( uint8_t *data, struct make_beam_opts *opts, ComplexDouble ***W,
                ComplexDouble ****J, int file_no, int nstation, int nchan,
                int npol, int outpol_coh, int outpol_incoh, double invw,
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
    int sample;
#pragma omp parallel for
    for (sample = 0; sample < (int)opts->sample_rate; sample++) {

        int ch, ant, pol, opol, opol1, opol2;
        int pfb, rec, inc;
        int data_idx;

        // Find out if any of the scenarios that need coherent
        // output are wanted.
        int coherent_requested = opts->out_coh    ||
                                 opts->out_vdif   ||
                                 opts->out_uvdif;

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

    } // End OpenMP parallel for

}


void form_stokes( ComplexDouble **detected_beam,
                  ComplexDouble noise_floor[][2][2],
                  int nchan, double invw, float *spectrum )
/* This function forms the Stokes parameters IQUV from the detected beam and
 * the constructed "noise floor". It packs the IQUV parameters into a single
 * 1D array in the following format (assuming nchan = 128):
 *
 *   I1, I2, I3, ..., I128,
 *   Q1, Q2, Q3, ..., Q128,
 *   U1, U2, U3, ..., U128,
 *   V1, V2, V3, ..., V128
 *
 * Where the numbers indicate the fine channel number. The above 4x128 block
 * represents a single time sample.
 *
 * The following array sizes are assumed:
 *
 *   detected_beam[nchan][2]
 *   noise_floor[nchan][2][2]
 *   spectrum[nchan*4]
 *
 * IQUV are all weighted by multiplication by "invw"
 * These calculations are described in .../doc/doc.pdf.
 */
{
    float beam00, beam11;
    float noise0, noise1, noise3;
    ComplexDouble beam01, beam01_n;
    unsigned int stokesIidx, stokesQidx, stokesUidx, stokesVidx;

    int ch;
    for (ch = 0; ch < nchan; ch++)
    {
        beam00 = CReald( CMuld( detected_beam[ch][0], CConjd(detected_beam[ch][0]) ) );
        beam11 = CReald( CMuld( detected_beam[ch][1], CConjd(detected_beam[ch][1]) ) );
        beam01 = CMuld( detected_beam[ch][0], CConjd(detected_beam[ch][1]) );

        noise0 = CReald( noise_floor[ch][0][0] );
        noise1 = CReald( noise_floor[ch][0][1] );
        noise3 = CReald( noise_floor[ch][1][1] );

        stokesIidx = 0*nchan + ch;
        stokesQidx = 1*nchan + ch;
        stokesUidx = 2*nchan + ch;
        stokesVidx = 3*nchan + ch;

        beam01_n = CSubd( beam01, CMaked( noise1, 0.0 ) );

        // Looking at the dspsr loader the expected order is <ntime><npol><nchan>
        // so for a single timestep we do not have to interleave - I could just stack these
        spectrum[stokesIidx] = (beam00 + beam11 - noise0 - noise3) * invw;
        spectrum[stokesQidx] = (beam00 - beam11 - noise0 + noise3) * invw;
        spectrum[stokesUidx] =  2.0 * CReald(beam01_n) * invw;
        spectrum[stokesVidx] = -2.0 * CImagd(beam01_n) * invw;
    }
}

