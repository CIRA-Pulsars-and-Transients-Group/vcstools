#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <plot.h>
#include <math.h>
#include <complex.h>
#include <inttypes.h>
#include <slamac.h>
#include <fftw3.h>
#include <assert.h>
#include <fcntl.h>
#include <unistd.h>


#ifndef __APPLE__
#include <omp.h>
#endif

#include "mwac_utils.h"
#include "antenna_mapping.h"

/* some externals for the ipfb call */

int nfrequency;
int npol;
int nstation;


map_t corr_mapping[NINPUT][NINPUT];
int pfb_output_to_input[NINPUT];
int miriad_to_mwac[256];
int single_pfb_mapping[64];

void cp2x2(complex double *Min, complex double *Mout) {
    Mout[0] = Min[0];
    Mout[1] = Min[1];
    Mout[2] = Min[2];
    Mout[3] = Min[3];
}

void inv2x2(complex double *Min, complex double *Mout) {
    complex double inv_det = 1.0 / (Min[0] * Min[3] - Min[1] * Min[2]);
    Mout[0] = inv_det * Min[3];
    Mout[1] = -inv_det * Min[1];
    Mout[2] = -inv_det * Min[2];
    Mout[3] = inv_det * Min[0];
}

void mult2x2d(complex double *M1, complex double *M2, complex double *Mout) {
    Mout[0] = M1[0] * M2[0] + M1[1] * M2[2];
    Mout[1] = M1[0] * M2[1] + M1[1] * M2[3];
    Mout[2] = M1[2] * M2[0] + M1[3] * M2[2];
    Mout[3] = M1[2] * M2[1] + M1[3] * M2[3];
}

void mult2x2t(complex double *M1, complex double *M2, complex double *Mout) {
    Mout[0] = M1[0] * M2[0] + M1[1] * M2[1];
    Mout[1] = M1[0] * M2[2] + M1[1] * M2[3];
    Mout[2] = M1[2] * M2[0] + M1[3] * M2[1];
    Mout[3] = M1[2] * M2[2] + M1[3] * M2[3];
}

void mult2x2h(complex double *M1, complex double *M2, complex double *Mout) {
    Mout[0] = M1[0] * conj(M2[0]) + M1[1] * conj(M2[1]);
    Mout[1] = M1[0] * conj(M2[2]) + M1[1] * conj(M2[3]);
    Mout[2] = M1[2] * conj(M2[0]) + M1[3] * conj(M2[1]);
    Mout[3] = M1[2] * conj(M2[2]) + M1[3] * conj(M2[3]);
}

void multaccum2x2dd(complex double *M1, complex double *M2, complex double *Mout) {
    Mout[0] += M1[0] * M2[0] + M1[1] * M2[2];
    Mout[1] += M1[0] * M2[1] + M1[1] * M2[3];
    Mout[2] += M1[2] * M2[0] + M1[3] * M2[2];
    Mout[3] += M1[2] * M2[1] + M1[3] * M2[3];
}

void multaccum2x2dt(complex double *M1, complex double *M2, complex double *Mout) {
    Mout[0] += M1[0] * M2[0] + M1[1] * M2[1];
    Mout[1] += M1[0] * M2[2] + M1[1] * M2[3];
    Mout[2] += M1[2] * M2[0] + M1[3] * M2[1];
    Mout[3] += M1[2] * M2[2] + M1[3] * M2[3];
}

void multaccum2x2dh(complex double *M1, complex double *M2, complex double *Mout) {
    Mout[0] += M1[0] * conj(M2[0]) + M1[1] * conj(M2[1]);
    Mout[1] += M1[0] * conj(M2[2]) + M1[1] * conj(M2[3]);
    Mout[2] += M1[2] * conj(M2[0]) + M1[3] * conj(M2[1]);
    Mout[3] += M1[2] * conj(M2[2]) + M1[3] * conj(M2[3]);
}

void multaccum2x2hd(complex double *M1, complex double *M2, complex double *Mout) {
    Mout[0] += conj(M1[0]) * M2[0] + conj(M1[2]) * M2[2];
    Mout[1] += conj(M1[0]) * M2[1] + conj(M1[2]) * M2[3];
    Mout[2] += conj(M1[1]) * M2[0] + conj(M1[3]) * M2[2];
    Mout[3] += conj(M1[1]) * M2[1] + conj(M1[3]) * M2[3];
}

void mult2x2tlum(complex double *M1, complex double *M2, complex double *M3, complex double *Mout) {
    complex double tmp[4];
    tmp[0] = M1[0] * M2[0] + M1[1] * M2[2];
    tmp[1] = M1[0] * M2[1] + M1[1] * M2[3];
    tmp[2] = M1[2] * M2[0] + M1[3] * M2[2];
    tmp[3] = M1[2] * M2[1] + M1[3] * M2[3];
    Mout[0] = tmp[0] * conj(M3[0]) + tmp[1] * conj(M3[1]);
    Mout[1] = tmp[0] * conj(M3[2]) + tmp[1] * conj(M3[3]);
    Mout[2] = tmp[2] * conj(M3[0]) + tmp[3] * conj(M3[1]);
    Mout[3] = tmp[2] * conj(M3[2]) + tmp[3] * conj(M3[3]);
}

void invert_pfb(complex float *input, complex float *output, int nchan_in, int nchan_out, int npol_in, int npol_out,int mode,int last, void *fcontext) {
    
    int edges = 0;
    int fftmode = 1;
  
    int ntaps = 0;
    int nsamples = 0;
    float *filter = 0;
    struct filter_context *f = (struct filter_context *) fcontext; 
    
    if (f != NULL) {
        ntaps = f->ntaps;
        nsamples = f->nsamples;
        filter = f->filter;
    }

    switch (mode) {
        
        case 1: {
            
            edges = 0;
            fftmode = 1;
            break;
            
        }
        case 2: {
            // with edges 20 each side
            edges = 20;
            fftmode = 1;
            break;
        }
        case 3:
            edges = 0;
            fftmode = 0;        
        default:
            break;
            
    }

    
    int nchan = nchan_in + 2*edges;
    int direction = FFTW_BACKWARD;
    int ch;
    int down_sample; // which sample we are on be


    int upsample_factor = (int) nchan_in/nchan_out;

    static fftwf_complex **in = NULL;
    static fftwf_complex **out = NULL;
    static fftwf_plan *p_forward = NULL;
    static fftwf_complex *filter_dft_in = NULL;
    static fftwf_complex *filter_dft_out = NULL;
    static fftwf_plan p,p_backward;
    static fftwf_complex *upsample_working = NULL;

    if (fftmode == 1) {
    // this simply does a backward FFT 
        if (in == NULL || out == NULL) {
            in = (fftwf_complex **) fftwf_malloc(sizeof(fftwf_complex *));
            out = (fftwf_complex **) fftwf_malloc(sizeof(fftwf_complex *));
            in[0] = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * nchan);
            out[0] = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * nchan);
            p = fftwf_plan_dft_1d(nchan,in[0],out[0],direction,FFTW_ESTIMATE);

        }
        if (in == NULL || out == NULL){
            fprintf(stderr,"Failed to allocate FFTW arrays\n");
            return;
        }


        bzero(in[0],nchan*sizeof(fftwf_complex));
        bzero(out[0],nchan*sizeof(fftwf_complex));

        for (ch=0;ch<nchan_in;ch++) {

            if (ch < nchan_in/2) {
                // these are the negative frequencies
                // pack them in the second half of the array
                // skip over the edges
                in[0][((nchan/2) + edges) + ch] = input[ch];
            }
            if (ch > nchan_in/2) {
                // positive frequencies -- shift them to the first half
                in[0][ch-(nchan_in/2)] = input[ch];
            }
            if (ch == nchan_in/2) {
                // Nyquist bin - give it a zero mean
                in[0][0] = 0.0;
            }


        }
        /*
           FILE *tmp;
           tmp = fopen("fft.txt","a");
           for (ch=0;ch<nchan;ch++){
           fprintf(tmp,"%d %f %f\n",ch,crealf(in[ch]),cimagf(in[ch]));
           }
           fclose(tmp);
           */

        fftwf_execute(p);

        for (ch=0;ch<nchan;ch++) {
            output[ch] = out[0][ch];

        }
        if (last == 1) {
            fftwf_destroy_plan(p);
            fftwf_free(in[0]);
            fftwf_free(in);
            fftwf_free(out[0]);
            fftwf_free(out);
            in = NULL;
            out = NULL;
        }
    }
    else if (fftmode == 0) { // the stitch-up mode
        
         fprintf(stdout,"In synthesis mode\n");    
        /* in order to build a single time step at the new <high> resolution
         * 
         * 1) first take each channel and upsample the correct factor by adding the upsample-1 number of zeros
         * 2) filter the upsampled channel to remove the images 
         * 3) phase rotate the channel to shift the frequency response
         * 4) Sum all to get the final time series
         * *
         */


        int out_sample;

        if (filter == NULL) {
            return;
        }
        if (ntaps == 0) {
            return;
        }
        if (nsamples == 0) {
            return;
        }
        fprintf(stderr,"Processing %d time samples\n",nsamples);
        if (upsample_working==NULL) {
            fprintf(stderr,"Making %d working buffer to hold filter product\n", upsample_factor*nsamples);
            upsample_working = (complex float *) malloc (sizeof(fftwf_complex)*upsample_factor*nsamples);
            fprintf(stderr,"Done\n");
        }
        if (filter_dft_in == NULL || filter_dft_out == NULL) {
            fprintf(stderr,"Building FFT of filter response\n");
            filter_dft_in = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * upsample_factor*nsamples);
            filter_dft_out = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * upsample_factor*nsamples);
     
            bzero(filter_dft_in,upsample_factor*nsamples*sizeof(fftwf_complex));
            int t = 0;
            for (t=0;t<ntaps;t++) {
                filter_dft_in[t] = filter[t] + I*0.0;
            }
            p = fftwf_plan_dft_1d(upsample_factor*nsamples,filter_dft_in,filter_dft_out,FFTW_FORWARD,FFTW_ESTIMATE);
            fftwf_execute(p);
            fprintf(stderr,"Done\n");
        }
        if (in == NULL || out == NULL) {
            // in order to permit some openMP lovAeliness I will make an in and out array for each channel
            in = (fftwf_complex **) fftwf_malloc(sizeof(fftwf_complex *) * nchan_in);
            out = (fftwf_complex **) fftwf_malloc(sizeof(fftwf_complex *) * nchan_in);
            p_forward = (fftwf_plan *) fftwf_malloc(sizeof(fftwf_plan) * nchan_in);

            for (ch = 0; ch<nchan_in; ch++){
                in[ch] = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * upsample_factor * nsamples);
                out[ch] = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * upsample_factor * nsamples);
                p_forward[ch] = fftwf_plan_dft_1d(upsample_factor*nsamples,in[ch],out[ch],FFTW_FORWARD,FFTW_ESTIMATE);
            }
            p_backward = fftwf_plan_dft_1d(upsample_factor*nsamples,upsample_working,out[0],FFTW_BACKWARD,FFTW_ESTIMATE);
        }
        // zero the output
        bzero(upsample_working,upsample_factor*nsamples*sizeof(fftwf_complex));
#ifndef __APPLE__
        fprintf(stderr,"Max threads available %d\n",omp_get_max_threads());
#pragma omp parallel for 
#endif        
        for (ch = 0; ch < nchan_in ; ch++) { // each channel is done separately
            // zero the input
            bzero(in[ch],upsample_factor*nsamples*sizeof(fftwf_complex));
            for (down_sample = 0; down_sample < nsamples; down_sample++) { // we need many samples
                // put in the good samples and allow the intervening zeros to be the upsampling
                in[ch][down_sample*upsample_factor] = input[down_sample*nchan_in + ch];
            }
            // Need to Fourier Transform the input
            //
            fftwf_execute(p_forward[ch]);
            // perform the convolution
         }
         for (ch = 0; ch < nchan_in ; ch++) { // each channel is done separately
            size_t out_index;  
            for (out_sample = 0; out_sample<(upsample_factor*nsamples);out_sample++) { 
                // we are working on the out_sample
                // it is made of the sum of the products of the filter taps
                // but rotate the convolution (ala shift theorem to place the RF in the correct spot
                out_index = out_sample+(ch*nsamples);
                if (out_index > (upsample_factor * nsamples)){
                    out_index = out_index - (upsample_factor * nsamples);
                }
                upsample_working[out_index] = upsample_working[out_index] + filter_dft_out[out_sample]*out[ch][out_sample];
                
            }

        }

        // the RF has been populated - I just have to invert back to the time domain and we are done;       
        fftwf_execute(p_backward);

        // memcpy into the output array
        //
        memcpy(output,out[0],(upsample_factor*nsamples));

        if (last == 1) {
            
            free(upsample_working);

            for (ch=0;ch<nchan_in;ch++){
                fftwf_free(in[ch]);
                fftwf_free(out[ch]);
            }
            fftwf_free(in);
            fftwf_free(out);
            fftwf_free(filter_dft_in);
            fftwf_free(filter_dft_out);

            upsample_working = NULL;
            in = NULL;
            out = NULL;
            filter_dft_in = NULL;
            filter_dft_out = NULL;

        }

    }
}
void dec2hms(char *out, double in, int sflag) {
    int sign = 1;
    char *ptr = out;
    int h, m;
    double s;
    if (in<0.0) { sign = -1; in = fabs(in); }
    h = (int)in; in -= (double)h; in *= 60.0;
    m = (int)in; in -= (double)m; in *= 60.0;
    s = in;
    if (sign==1 && sflag) { *ptr='+'; ptr++; }
    else if (sign==-1) { *ptr='-'; ptr++; }
    sprintf(ptr, "%2.2d:%2.2d:%07.4f", h, m, s);
}
int read_casa_gains_file(char *gains_file,complex double **antenna_gain, int nant, int chan_to_get) {
    FILE *fp = NULL;
    fp = fopen(gains_file,"r");
    if (fp == NULL) {
        fprintf(stderr,"Failed to open %s: quitting\n",gains_file);
        exit(0);
    }
    int linenum=0;
    char line[BUFSIZE];
    char time_str[BUFSIZE];
    char target[BUFSIZE];
    int channel;
    char entry[24][BUFSIZE];
    int nscan;
    int input = 0;

    *antenna_gain = (complex double *) calloc(nant*2,sizeof(complex double)); // X,Y
    if (*antenna_gain == NULL) {
        fprintf(stderr,"Failed to allocate memory for complex gains \n");
        exit(-1);
    }
    input = 0;
    while ((fgets(line, BUFSIZE - 1, fp)) != NULL) {

        if (line[0] == '\n' || line[0] == '#' || line[0] == '\0' || line[0] == '-' || line[0] == 'S' || line[0] =='T' || line[0]==' ' || line[0] == 'L')
            continue; // skip blank/comment lines
        if (line[0] == '/' && line[1] == '/')
            continue; // also a comment (to match other input files using this style)

        nscan = sscanf(line,"%s %s %d|%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",time_str,target,&channel,
                       entry[0],entry[1],entry[2],
                       entry[3],entry[4],entry[5],
                       entry[6],entry[7],entry[8],
                       entry[9],entry[10],entry[11],
                       entry[12],entry[13],entry[14],
                       entry[15],entry[16],entry[17],
                       entry[18],entry[19],entry[20],
                       entry[21],entry[22],entry[23]);


        if (channel == chan_to_get) {
            channel = -99;

            fprintf(stdout,"CASA: %s\n",line);
            int ent_ind = 0;
            while (ent_ind < 24) {
                fprintf(stdout,"CASA:Entry[%d] = %s\n",ent_ind,entry[ent_ind]);
                ent_ind++;
            }
            int ii = 0;
            while (ii < nscan-3 && input < 256) {
                if (strcmp(entry[ii+2],"F") != 0) {
                    (*antenna_gain)[input] = atof(entry[ii])*cexp(I*M_PI*atof(entry[ii+1])/180.0);

                    fprintf(stdout,"CASA: Input %d -- Amp: %f Phase %f \n",input,atof(entry[ii]),M_PI*atof(entry[ii+1])/180.0);
                    ii=ii+2;
                }
                else {
                    (*antenna_gain)[input] = 0.0 + I*0.0;

                    fprintf(stdout,"CASA: Input %d - FLAGGED\n",input);
                    ii = ii+3;
                }
                input++;
            }




        }
        linenum++;

    }
    fclose(fp);
    fprintf(stdout,"Read %d inputs from %s\n",input,gains_file);
    int ii;
    for (ii=0;ii<nant*2;ii++){
        fprintf(stdout,"Complex Gain Correction for input[%d] channel[%d] is %lf %lf\n",ii,chan_to_get,creal((*antenna_gain)[ii]),cimag((*antenna_gain)[ii]));

    }
    return input/2;
}

int read_offringa_gains_file(complex double **antenna_gain, int nant, int coarse_chan, char *gains_file) {
    // Assumes that memory for antenna has already been allocated

    // Open the calibration file for reading
    FILE *fp = NULL;
    fp = fopen(gains_file,"r");
    if (fp == NULL) {
        fprintf(stderr,"Failed to open %s: quitting\n",gains_file);
        exit(1);
    }

    // Read in the necessary information from the header

    uint32_t intervalCount, antennaCount, channelCount, polarizationCount;

    fseek(fp, 16, SEEK_SET);
    fread(&intervalCount,     sizeof(uint32_t), 1, fp);
    fread(&antennaCount,      sizeof(uint32_t), 1, fp);
    fread(&channelCount,      sizeof(uint32_t), 1, fp);
    fread(&polarizationCount, sizeof(uint32_t), 1, fp);

    // Error-checking the info extracted from the header
    if (intervalCount > 1) {
        fprintf(stderr, "Warning: Only the first interval in the calibration ");
        fprintf(stderr, "solution (%s) will be used\n", gains_file);
    }
    if (antennaCount != nant) {
        fprintf(stderr, "Error: Calibration solution (%s) ", gains_file);
        fprintf(stderr, "contains a different number of antennas (%d) ", antennaCount);
        fprintf(stderr, "than specified (%d)\n", nant);
        exit(1);
    }
    if (channelCount != 24) {
        fprintf(stderr, "Warning: Calibration solution (%s) ", gains_file);
        fprintf(stderr, "contains a different number (%d) ", channelCount);
        fprintf(stderr, "than the expected (%d) channels. ", 24);
    }
    if (channelCount <= coarse_chan) {
        fprintf(stderr, "Error: Requested channel number (%d) ", coarse_chan);
        fprintf(stderr, "is more than the number of channels (0-%d) ", channelCount-1);
        fprintf(stderr, "available in the calibration solution (%s)\n", gains_file);
        exit(1);
    }
    int npols = polarizationCount; // This will always = 4

    // Prepare to jump to the first solution to be read in
    int bytes_left_in_header = 16;
    int bytes_to_first_jones = bytes_left_in_header + (npols * coarse_chan * sizeof(complex double));
         //     (See Offringa's specs for details)
         //     Assumes coarse_chan is zero-offset
         //     sizeof(complex double) *must* be 64-bit x 2 = 16-byte
    int bytes_to_next_jones = npols * (channelCount-1) * sizeof(complex double);

    int ant, pol;  // Iterate through antennas and polarisations
    //int ant_idx;   // Used for "re-ordering" the antennas. If not needed (after testing), DELETE ME
    int count = 0; // Keep track of how many solutions have actually been read in
    double re, im;    // Temporary placeholders for the real and imaginary doubles read in

    int first = 1;
    for (ant = 0; ant < nant; ant++) {

        // Jump to next Jones matrix position for this channel
        if (first) {
            fseek(fp, bytes_to_first_jones, SEEK_CUR);
            first = 0;
        }
        else {
            fseek(fp, bytes_to_next_jones, SEEK_CUR);
        }

        // Reorder the antennas. Not necessarily needed. If not, DELETE ME
        //ant_idx = natural_to_mwac[ant*2]/2;

        // Read in the data
        for (pol = 0; pol < npols; pol++) {

            fread(&re, sizeof(double), 1, fp);
            fread(&im, sizeof(double), 1, fp);

            // Check for NaNs
            if (isnan(re) | isnan(im)) {
                antenna_gain[ant][pol] = 0.0 + I*0.0;
                //antenna_gain[ant_idx][pol] = 0.0 + I*0.0; // DELETE ME unless antenna ordering DOES need to be changed
            }
            else {
                antenna_gain[ant][pol] = re  + I*im;
                //antenna_gain[ant_idx][pol] = re  + I*im; // DELETE ME unless antenna ordering DOES need to be changed
            }

            count++;

        }
    }

    // Close the file, print a summary, and return
    fclose(fp);
    fprintf(stdout, "Read %d inputs from %s\n", count, gains_file);

    return count/npols; // should equal the number of antennas
}

int gain_file_id(char *gains_file) {
    FILE *fp = NULL;
    int rvalue=0;
    fp = fopen(gains_file,"r");
    if (fp == NULL) {
        fprintf(stderr,"Failed to open %s: quitting\n",gains_file);
        exit(0);
    }
    char line[BUFSIZE];
    char key[BUFSIZE];

    while ((fgets(line, BUFSIZE - 1, fp)) != NULL) {

        sscanf(line,"%s",key);
        fprintf(stdout,"GET FILE ID: %s\n",key);
        if (strncmp("Spw",key,3) == 0) {
            rvalue = CASA_GAINS_FILE;
            break;
        }
        else {
            rvalue = MIRIAD_GAINS_FILE;
            break;

        }
    }
    fclose(fp);
    return rvalue;
    
}
int read_miriad_gains_file(char *gains_file, complex double **antenna_gain){
    FILE *fp = NULL;
    fp = fopen(gains_file,"r");
    if (fp == NULL) {
        fprintf(stderr,"Failed to open %s: quitting\n",gains_file);
        exit(0);
    }
    char line[BUFSIZE];
    char time_str[BUFSIZE];
    int nscan=0;
    int nant = 0;
    int count = 0;
    int linenum=0;
    int i=0;
    
    float entry[6];
    int index=0;
    
    int real_values=1;
    while ((fgets(line, BUFSIZE - 1, fp)) != NULL) {
        linenum++;
        if (linenum == 4) {
            //parse for antenna number
            nscan = sscanf(line,"%c %s %s %s %d",time_str,time_str,time_str,time_str,&nant);
            if (nscan != 5) {
                fprintf(stderr,"Failed to parse NANT from MIRIAD file %s\n",gains_file);
                exit(0);
            }
            else {
                fprintf(stdout,"Parsed NANT = %d\n",nant);
            }
           
            *antenna_gain = (complex double *) calloc(nant*2,sizeof(complex double)); // X,Y
            if (*antenna_gain == NULL) {
                fprintf(stderr,"Failed to allocate memory for complex gains \n");
                exit(-1);
            }
        }
        if (line[0] == '\n' || line[0] == '#' || line[0] == '\0')
            continue; // skip blank/comment lines
        if (line[0] == '/' && line[1] == '/')
            continue; // also a comment (to match other input files using this style)
        
        if (index == 0) {
            // read the first line this also contains the time etc
            nscan = sscanf(line, "%d %s %f %f %f %f %f %f",&count,time_str,&entry[0],&entry[1],&entry[2],&entry[3],&entry[4],&entry[5]);
            if (nscan != 8) {
                fprintf(stderr,"Failed to parse first data line from MIRIAD file %s\n",gains_file);
                exit(0);
            }
            else {
                fprintf(stdout,"Parsed first data line of miriad file %d %s\n",count,time_str);
                nscan = nscan - 2;
            }
            
        }
        else {
            nscan = sscanf(line, "%f %f %f %f %f %f",&entry[0],&entry[1],&entry[2],&entry[3],&entry[4],&entry[5]);
        }
        
        for (i=0;i<nscan;i++) {
            if (real_values) {
                (*antenna_gain)[index+i] = entry[i];
            }
            else {
                (*antenna_gain)[index+i] = (*antenna_gain)[index+i] + I*entry[i];
            }
        }
        index = index + nscan;
        
        fprintf(stdout,"Parsed index %d (real=%d)\n",index,real_values);
        
        if (index == 2*nant) {
            index = 0;
            real_values = 0;
        }
        else if (index == 4*nant) {
            break;
        }
        
    }
    for (i=0;i<nant*2;i=i+2) {
        fprintf(stdout,"Gain Ant  %d (X) real %f imag %f\n",i/2,crealf((*antenna_gain)[i]),cimagf((*antenna_gain)[i]));
        fprintf(stdout,"Gain Ant  %d (Y) real %f imag %f\n",i/2,crealf((*antenna_gain)[i+1]),cimagf((*antenna_gain)[i+1]));
        
    }
    fclose(fp);
    return nant;
}
int read_rts_file(complex double **G, complex double *Jref, int nant, double *amp,
                      char *fname) {
    
    FILE *fp = NULL;
    if ((fp = fopen(fname, "r")) == NULL) {
        fprintf(stderr, "Error: cannot open gain Jones matrix file: %s\n",
                fname);
        exit(0);
    }
    
    char line[BUFSIZE];
    int index = 0, nscan;
    double re0, im0, re1, im1, re2, im2, re3, im3;
    
    while ((fgets(line, BUFSIZE - 1, fp)) != NULL) {
        
        if (line[0] == '\n' || line[0] == '#' || line[0] == '\0')
            continue; // skip blank/comment lines
        if (line[0] == '/' && line[1] == '/')
            continue; // also a comment (to match other input files using this style)
        
        if (index == 0) {
            
            // read the amplitude and the Alignment Line
            nscan = sscanf(line, "%lf", amp);
            fgets(line, BUFSIZE - 1, fp);
            nscan = sscanf(line, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf", &re0,
                           &im0, &re1, &im1, &re2, &im2, &re3, &im3);
            
            Jref[0] = re0 + im0 * I;
            Jref[1] = re1 + im1 * I;
            Jref[2] = re2 + im2 * I;
            Jref[3] = re3 + im3 * I;
            
        }
        if (index > 0) {
            nscan = sscanf(line, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf", &re0,
                           &im0, &re1, &im1, &re2, &im2, &re3, &im3);
            G[index - 1][0] = re0 + im0 * I;
            G[index - 1][1] = re1 + im1 * I;
            G[index - 1][2] = re2 + im2 * I;
            G[index - 1][3] = re3 + im3 * I;
            
        }
        
        index++;
        
    }
    
    fclose(fp);
    
    return 0;
    
} /* read_cal_file */

int default_read_pfb_call(int in_fd, int out_fd, char *heap) {

// Initialise
//
    int8_t *input_buffer = NULL;
    int8_t *output_buffer = NULL;
    int8_t *binary_buffer = NULL;

    size_t nsteps_read=0;
    
    nfrequency=128;
    npol=2;
    nstation=128;

    input_buffer = malloc((nstation*npol*nfrequency*2*4)/8);
    assert(input_buffer);
    output_buffer = malloc(nstation*npol*nfrequency*2);
    assert(output_buffer);
    binary_buffer = malloc(nstation*npol*nfrequency*2);
    assert(binary_buffer);

    char *heap_ptr = heap;
    
    int index = 0;

    int out_index = 0;

    size_t rtn = 0;
    size_t items_to_read = 0;


    uint8_t original;
    int8_t value;
    static int filled = 0;

    if (filled == 0) {
        fill_mapping_matrix();
        filled = 1;
    }

    size_t gulp = (nstation*npol*nfrequency*2*4)/8;
    int end_of_file = 0;

    while (!end_of_file) {

        items_to_read = gulp;
        rtn = 0;

        while(items_to_read > 0) {
            rtn = read(in_fd,input_buffer+(gulp-items_to_read),items_to_read);
            if (rtn == 0) {
                end_of_file = 1;
                break;
            }
            else if (rtn > 0) {
                
                // fprintf(stderr,"read: got %lu expected %lu\n",rtn,items_to_read);
                items_to_read -= rtn;
                
            }
            else {
                fprintf(stderr,"File conversion error on read\n");
                return -1;
            }
        }
        
        
        
        if (end_of_file && (items_to_read == gulp)){
            fprintf(stderr,"heap reader has read %d steps\n",nsteps_read);
            break;
        }
        else if (end_of_file && (items_to_read != gulp)){
            fprintf(stderr,"EOF on input file with incomplete read <FATAL>\n");
            return -1;
        }
        else {
            nsteps_read++;
        }
        // time then frequency runs slowest in this data block

        int8_t *inp_ptr = NULL;
        int map_index = 0;
        int chan = 0;
        int pol = 0;

        for (chan = 0 ; chan < nfrequency ; chan++) {

            inp_ptr = (int8_t *) input_buffer + ((chan*nstation*npol*2*4)/8);
            for (index = 0; index < nstation; index++) {

                for (pol =0 ;pol < npol; pol++) {
                    // this returns the desired output index (0 to 64)
                    map_index = pfb_output_to_input[index*npol + pol];
                    // output matrix index - this index reorders the
                    // stations AND transposes
                    // the order will now be [time][input][chan][complexity]
                    //
                    out_index = map_index*nfrequency*2 + chan*2;

                    // .. assigns the value to be the current de-ref input_ptr


                    original = (uint8_t) *inp_ptr;
                    original = original & 0xf; // first sample
                    if (original >= 0x8) {
                        value = original - 0x10;
                    }
                    else {
                        value = original;
                    }

                    output_buffer[out_index] = (int8_t) (value & 0xff);


                    out_index++;

                    original = (uint8_t) *inp_ptr ;
                    original = original >> 4;
                    original = original & 0xf; // second sample
                    if (original >= 0x8) {
                        value = original - 0x10;
                    }
                    else {
                        value = original;
                    }
                    output_buffer[out_index] = (int8_t) (value & 0xff);

                    inp_ptr++;


                } // for all pol
            } // for all stations
        } // for all channels
        //}


        int8_t *out_ptr;
        size_t out_size;
        out_ptr = output_buffer;
        out_size = 2*nstation*npol*nfrequency;

        size_t items_to_write = out_size;
        rtn = 0;
        while( items_to_write > 0) {
            if (out_fd > 0) {
                rtn = write(out_fd,out_ptr+(out_size-items_to_write),items_to_write);
                items_to_write -= rtn;
            }
            if (heap != NULL) {
                memcpy(heap_ptr,out_ptr,items_to_write);
                heap_ptr = heap_ptr+items_to_write;
                items_to_write = 0;
            }
        }

    } // next time step


    free(input_buffer);
    free(output_buffer);
    free(binary_buffer);

    return 0;

}
int read_cal_file(complex double **G, int nant, double *amp) {

    FILE *fp = NULL;
    char line[BUFSIZE];
    int index = 0, nscan;
    double re0, im0, re1, im1, re2, im2, re3, im3;
    complex double A[4], invA[4], tmp[4];

    if ((fp = fopen("Gjones.dat", "r")) == NULL) {
        fprintf(stderr,
                "Error: cannot open gain Jones matrix file: Gjones.dat\n");
        exit(0);
    }

    while ((fgets(line, BUFSIZE - 1, fp)) != NULL) {

        if (line[0] == '\n' || line[0] == '#' || line[0] == '\0')
            continue; // skip blank/comment lines
        if (line[0] == '/' && line[1] == '/')
            continue; // also a comment (to match other input files using this style)
        
        if (index == 0) {
            
            // read the amplitude and the Alignment Line
            nscan = sscanf(line, "%lf", amp);
            fgets(line, BUFSIZE - 1, fp);
            nscan = sscanf(line, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf", &re0,
                           &im0, &re1, &im1, &re2, &im2, &re3, &im3);
            
            A[0] = re0 + im0 * I;
            A[1] = re1 + im1 * I;
            A[2] = re2 + im2 * I;
            A[3] = re3 + im3 * I;
            
            inv2x2(A, invA);
            
        }
        if (index > 0) {
            nscan = sscanf(line, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf", &re0,
                           &im0, &re1, &im1, &re2, &im2, &re3, &im3);
            tmp[0] = re0 + im0 * I;
            tmp[1] = re1 + im1 * I;
            tmp[2] = re2 + im2 * I;
            tmp[3] = re3 + im3 * I;
            
            mult2x2d(tmp, invA, G[index - 1]);
            
        }
        
        index++;
        
    }
    
    fclose(fp);
    
    return 0;
    
} /* read_cal_file */


void fill_mapping_matrix() {

	int inp1 = 0, inp2 = 0;
	int pol1 = 0, pol2 = 0;
	int index1 = 0, index2 = 0;

	int p=0,npfb = 4;
	

	
	//	 Output matrix has ordering
	//	 [channel][station][station][polarization][polarization][complexity]

	for (p=0;p<npfb;p++) {
		for (inp1=0;inp1<64;inp1++) {
			pfb_output_to_input[(p*64) + inp1] = single_pfb_mapping[inp1] + (p*64);
            //fprintf(stdout,"input %d - maps to antenna %d\n",((p*64) + inp1),pfb_output_to_input[(p*64) + inp1]);
		}
	}

			 
	for (inp1 = 0; inp1 < nstation; inp1++) {
		for (inp2 = 0; inp2 < nstation; inp2++) {
			for (pol1 = 0; pol1 < npol; pol1++) {
				for (pol2 = 0; pol2 < npol; pol2++) {
					index1 = inp1 * npol + pol1;
					index2 = inp2 * npol + pol2;
					/*
					   fprintf(stdout,
					   "inp1 %d pol1 %d inp2 %d pol2 %d map to index1 %d and index2 %d\n",
					   inp1, pol1, inp2, pol2, index1, index2);
					   fprintf(stdout,
					   "these map to PFB input numbers: %d and %d\n",
					   pfb_output_to_input[index1],
					   pfb_output_to_input[index2]);
					   */
					corr_mapping[pfb_output_to_input[index1]][pfb_output_to_input[index2]].stn1 =
						inp1; // this should give us the pfb input
					corr_mapping[pfb_output_to_input[index1]][pfb_output_to_input[index2]].stn2 =
						inp2;
					corr_mapping[pfb_output_to_input[index1]][pfb_output_to_input[index2]].pol1 =
						pol1;
					corr_mapping[pfb_output_to_input[index1]][pfb_output_to_input[index2]].pol2 =
						pol2;

				}
			}
		}
	}

}

void get_baseline(int st1, int st2, int pol1, int pol2, complex float *data,
		complex float *baseline) {

	int i, j, k, l, m;
	complex float *in, *out;
	extern int npol;
	extern int nstation;
	extern int nfrequency;
	in = data;
	out = baseline;

	for (i = 0; i < nfrequency; i++) {
		for (j = 0; j < nstation; j++) {
			for (k = 0; k < nstation; k++) {
				for (l = 0; l < npol; l++) {
					for (m = 0; m < npol; m++) {
						if (j == st1 && k == st2) {
							if (l == pol1 && m == pol2) {

								*out = *in;
								out++;
								// fprintf(stdout,"%f %f\n",crealf(*in),cimagf(*in));
							}

						}
						in++;
					}

				}
			}
		}
	}
}


void get_baseline_lu(int st1, int st2, int pol1, int pol2, float complex *data,
		float complex *baseline) {

	int i=0;
	float complex *in, *out;	
        
        extern int npol;
	extern int nstation;
	extern int nfrequency;

	size_t in_index=0;

	in = data;
	out = baseline;

	/* direct lookup */

	for (i=0;i<nfrequency;i++) {
		in_index = i*(nstation*nstation*npol*npol) + (st1*nstation*npol*npol) + (st2*npol*npol) + (pol1*npol) + pol2;
		*out = in[in_index];
		out++;
	}

}

void get_baseline_r(int st1, int st2, int pol1, int pol2, complex float *data,
		complex float *reorder,int npol, int nstation, int nfrequency,int true_st1,int true_st2,
		int true_pol1,int true_pol2,int conjugate) {

	int i=0;
	complex float *in, *out;
	size_t out_index =0, in_index=0;;
	in = data;
	out = reorder;
	
/* direct lookup */

	for (i=0;i<nfrequency;i++) {

		in_index = i*(nstation*nstation*npol*npol) + (st1*nstation*npol*npol) + (st2*npol*npol) + (pol1*npol) + pol2;
		out_index = i*(nstation*(nstation+1)*npol*npol/2) + (((true_st1*nstation) - ((true_st1+1)/2)*true_st1) + true_st2)*npol*npol + (pol1*npol) + pol2;
		if (!conjugate) {
			out[out_index] = in[in_index];
		}
		else {
			if (st2>st1) {
				out[out_index] = conj(in[in_index]);
			}
		}
	}

}
// full reorder using the correct mapping - takes the input cube and produces a packed triangular output
// in the correct order

// wacky packed tile order to packed triangular

void full_reorder(complex float *full_matrix_h, complex float *reordered)
{


	int t1=0;
	int t2=0;
	int p1=0;
	int p2=0;


	long long baseline_count = 0;


	for (t1 = 0; t1 < nstation; t1++) {
		for (t2 = t1; t2 < nstation; t2++) {
			for (p1 = 0;p1 < npol;p1++) {
				for (p2 =0; p2 < npol; p2++) {
					baseline_count++;

					int index1 = t1 * npol + p1;
					int index2 = t2 * npol + p2;
					/*
					   fprintf(stdout, "requesting ant1 %d ant 2 %d pol1 %d pol2 %d",
					   antenna1, antenna2, pol1, pol2);
					 */
					map_t the_mapping = corr_mapping[index1][index2];
					int conjugate = 0;
					/*
					   fprintf(stdout,
					   "input ant/pol combination decodes to stn1 %d stn2 %d pol1 %d pol2 %d\n",
					   the_mapping.stn1, the_mapping.stn2, the_mapping.pol1,
					   the_mapping.pol2);
					 */


					if (the_mapping.stn2 > the_mapping.stn1) {
						conjugate = 1;
					}
					else {
						conjugate = 0;
					}

					get_baseline_r(the_mapping.stn1, the_mapping.stn2, the_mapping.pol1,
							the_mapping.pol2, full_matrix_h, reordered,npol,nstation,nfrequency,conjugate,t1,t2,p1,p2);

				}
			}
		}
	}

	// now reoredered should contain a triagular packed array in the correct order
}
// Extracts the full matrix from the packed Hermitian form
void extractMatrix(complex float *matrix, complex float *packed) {
	int f, i, j, pol1, pol2;

	extern int npol;
	extern int nstation;
	extern int nfrequency;

	for (f = 0; f < nfrequency; f++) {
		for (i = 0; i < nstation; i++) {
			for (j = 0; j <= i; j++) {
				int k = f * (nstation + 1) * (nstation / 2) + i * (i + 1) / 2
					+ j;
				for (pol1 = 0; pol1 < npol; pol1++) {
					for (pol2 = 0; pol2 < npol; pol2++) {
						int index = (k * npol + pol1) * npol + pol2;
						matrix[(((f * nstation + i) * nstation + j) * npol
								+ pol1) * npol + pol2] = packed[index];
						matrix[(((f * nstation + j) * nstation + i) * npol
								+ pol2) * npol + pol1] = conjf(packed[index]);
					//	printf("f:%d s1:%d s2:%d %d p1:%d p2:%d %d\n",f,i,j,k,pol1,pol2,index);
					}
				}
			}
		}
	}

}
void extractMatrix_slow(complex float *matrix, complex float *packed) {
	int f, i, j, pol1, pol2;

	extern int npol;
	extern int nstation;
	extern int nfrequency;
	int in_index=0;
	int out_index=0;
	int out_index_conj=0;

	for (f = 0; f < nfrequency; f++) {
		for (i = 0; i < nstation; i++) {
			for (j = 0; j <= i; j++) {
				for (pol1 = 0; pol1 < npol; pol1++) {
					for (pol2 = 0; pol2 < npol; pol2++) {

						out_index = f*(nstation*nstation*npol*npol) + i*(nstation*npol*npol) + j*(npol*npol) + pol1*(npol) + pol2;
						out_index_conj = f*(nstation*nstation*npol*npol) + j*(nstation*npol*npol) + i*(npol*npol) + pol1*(npol) + pol2;
 						matrix[out_index] = packed[in_index];
						matrix[out_index_conj] = conjf(packed[in_index]);
						in_index++;
					}
				}
			}
		}
	}

}
