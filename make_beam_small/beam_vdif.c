#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "beam_vdif.h"
#include "mwa_header.h"
#include "ascii_header.h"

void vdif_write_second( vdifinfo *vf, int8_t *output )
{
    // form the filename
    // there is a standard naming convention
    char  filename[1030];
    sprintf( filename, "%s.vdif", vf->basefilename );

    //fprintf(stderr,"Attempting to open VDIF file for writing: %s\n",filename);
    FILE *fs = fopen( filename, "a" );
    fwrite( output, vf->block_size ,1, fs );
    fclose( fs );

    // write a CPSR2 test header for DSPSR

    char ascii_header[MWA_HEADER_SIZE] = MWA_HEADER_INIT;
    //ascii_header_set( ascii_header, "UTC_START", "%s", vf->date_obs  );
    ascii_header_set( ascii_header, "DATAFILE",   "%s", filename      );
    ascii_header_set( ascii_header, "INSTRUMENT", "%s", "VDIF"        );
    ascii_header_set( ascii_header, "TELESCOPE",  "%s", vf->telescope );
    ascii_header_set( ascii_header, "MODE",       "%s", vf->obs_mode  );
    ascii_header_set( ascii_header, "FREQ",       "%f", vf->fctr      );

    ascii_header_set( ascii_header, "BW",         "%f", vf->BW        );
    ascii_header_set( ascii_header, "RA",         "%s", vf->ra_str    );
    ascii_header_set( ascii_header, "DEC",        "%s", vf->dec_str   );
    ascii_header_set( ascii_header, "SOURCE",     "%s", vf->source    );

    sprintf( filename, "%s.hdr", vf->basefilename );
    fs = fopen( filename,"w" );
    fwrite( ascii_header, MWA_HEADER_SIZE, 1, fs );
    fclose( fs );

}


complex float get_std_dev_complex(complex float *input, int nsamples)
{
    // assume zero mean
    float rtotal = 0;
    float itotal = 0;
    float isigma = 0;
    float rsigma = 0;
    int i;

    for (i=0;i<nsamples;i++){
         rtotal = rtotal+(crealf(input[i])*crealf(input[i]));
         itotal = itotal+(cimagf(input[i])*cimagf(input[i]));

     }
    rsigma = sqrtf((1.0/(nsamples-1))*rtotal);
    isigma = sqrtf((1.0/(nsamples-1))*itotal);

    return rsigma+I*isigma;
}

void set_level_occupancy(complex float *input, int nsamples, float *new_gain)
{
    float percentage = 0.0;
    //float occupancy = 17.0;
    float limit = 0.00001;
    float step = 0.001;
    int i = 0;
    float gain = *new_gain;

    float percentage_clipped = 100;
    while (percentage_clipped > 0 && percentage_clipped > limit) {
        int count = 0;
        int clipped = 0;
        for (i=0;i<nsamples;i++) {
            if (gain*creal(input[i]) >= 0 && gain*creal(input[i]) < 64) {
                count++;
            }
            if (fabs(gain*creal(input[i])) > 127) {
                clipped++;
            }
        }
        percentage_clipped = ((float) clipped/nsamples) * 100;
        if (percentage_clipped < limit) {
            gain = gain + step;
        }
        else {
            gain = gain - step;
        }
        percentage = ((float)count/nsamples)*100.0;
        fprintf(stdout,"Gain set to %f (linear)\n",gain);
        fprintf(stdout,"percentage of samples in the first 64 (+ve) levels - %f percent \n",percentage);
        fprintf(stdout,"percentage clipped %f percent\n",percentage_clipped);
    }
    *new_gain = gain;
}


void get_mean_complex(complex float *input, int nsamples, float *rmean,float *imean, complex float *cmean)
{

    int i=0;
    float rtotal = 0;
    float itotal = 0 ;
    complex float ctotal = 0 + I*0.0;
    for (i=0;i<nsamples;i++){
        rtotal = rtotal+crealf(input[i]);
        itotal = itotal+cimagf(input[i]);
        ctotal = ctotal + input[i];
    }
    *rmean=rtotal/nsamples;
    *imean=itotal/nsamples;
    *cmean=ctotal/nsamples;

}

void normalise_complex(complex float *input, int nsamples, float scale)
{
    int i=0;

    for (i=0;i<nsamples;i++){
        input[i]=input[i]*scale;
    }
}


