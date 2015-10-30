/*
 * build_lfiles.c
 *
 *  Created on: Jul 20, 2012
 *      Author: sord
 */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <complex.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include "mwac_utils.h"
#include "antenna_mapping.h"

#include "fitsio.h"

#define MAX_FILES   24
#define MAX_CHAN    3072

/* globals */
FILE *fpd;        //file hand for debug messages. Set this in main()
int npol;
int nstation;
int nfrequency;
int debug=0;
char *input_file[MAX_FILES];

void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];

    if (status)
        fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    fits_get_errstatus(status, status_str);   /* get the error description */
    fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);

    /* get first message; null if stack is empty */
    if ( fits_read_errmsg(errmsg) )
    {
        fprintf(stderr, "\nError message stack:\n");
        fprintf(stderr, " %s\n", errmsg);

        while ( fits_read_errmsg(errmsg) )  /* get remaining messages */
            fprintf(stderr, " %s\n", errmsg);
    }

    exit( status );       /* terminate the program, returning error status */
}

void usage() {

    fprintf(stdout,"build_lfiles:\n");
    fprintf(stdout,"create lfiles (autos and X) from a correlator dump (mode 0) or an NGAS FITS file (mode 1)\n");
    fprintf(stdout,"Options:\n");
    fprintf(stdout,"\t-m\t 0|1|2 [default == 0]: read a raw correlator dump\n\t\t 1: read an NGAS FITS FILE\n\t\t 2: Read an OLD format raw correlator dump\n");
    fprintf(stdout,"\t-v\t <visibility file name> -- this argument can appear multiple times. Each file is concatenated in frequency <NOT TIME>\n");
    fprintf(stdout,"\t-f\t <nfrequency> [default: determine from input filess]: Number of channels <must equal nfiles*nchan_per_file>\n");
    fprintf(stdout,"\t-i\t <number of inputs> [default: determine from input files] : number of input that have been correlated <nstation*npol>\n");
    fprintf(stdout,"\t-o\t <output filename root> [default == \"last_dump\"] root filename of the output\n");
    fprintf(stdout,"\t-s\t <start> -- seconds from beginning of file to start. Default: 0\n");
    fprintf(stdout,"\t-n\t <nseconds> -- how many seconds to process. Default: all\n");
    fprintf(stdout,"\t-T\t <factor> [default == 1] -- by what factor to average in time\n");
    fprintf(stdout,"\t-F\t <factor> [default == 1] -- by what factor to average in frequency (must be power of 2)\n");
    fprintf(stdout,"\t-P\t <type>   Output product type. A: autos, C: cross, B: both. Default: both\n");
    fprintf(stdout,"\t-p\t -- image in primary header -- NOW ASSUMED IF MODE == 0\n");
    fprintf(stdout,"\t-a\t append to output instead of clobber\n");
    fprintf(stdout,"\t-d\t enable debugging messages\n");

    exit(EXIT_FAILURE);
}


/* get info about the number of HDUs and data type in each input file */
int getHDUsInfo(int nfiles, fitsfile *fptr[MAX_FILES], int num_hdus_in_file[MAX_FILES], long dims[2]) {
    int status=0,i,hdutype=0, ncols=0;
    long nrows=0;

    for (i=0; i< nfiles; i++) {
        if (fits_get_num_hdus(fptr[i],num_hdus_in_file+i,&status)) {
            printerror(status);
            return EXIT_FAILURE;
        }
        if (debug) fprintf(fpd,"There are %d HDUs in file %d\n",num_hdus_in_file[i],i);
    }

    /* sanity checks */
    for (i=1; i< nfiles; i++) {
        if ( fabs((num_hdus_in_file[i]-num_hdus_in_file[i-1])/(num_hdus_in_file[i]+1.0)) > 0.1) {
            fprintf(stderr,"ERROR: serious mismatch between number of HDUs in files %s (%d) vs %s (%d). Exiting.\n",
                            input_file[i],num_hdus_in_file[i], input_file[i-1], num_hdus_in_file[i-1]);
            return EXIT_FAILURE;
        }
    }

    /* from the first file, extract the data dimensions.
        In raw output files, this is in the primary image.
        In concatenated archive files, this is in the second and subsequent extensions
        In the old style O-files, this is in a binary table extension */
    if (num_hdus_in_file[0] ==1 ) {
        fits_movabs_hdu(fptr[0], 1, &hdutype, &status);
    }
    else {
        fits_movabs_hdu(fptr[0], 2, &hdutype, &status);
    }
    if (hdutype == BINARY_TBL) {
        fits_get_num_rows(fptr[0], &nrows, &status);
        fits_get_num_cols(fptr[0], &ncols, &status);
        if (status) {
            fits_report_error(stderr,status);
            return status;
        }
        dims[0]=ncols;
        dims[1]=nrows;
        if (debug) fprintf(fpd,"Found %d cols and %ld rows in binary table\n",ncols,nrows);
    }
    else {
        fits_get_img_size(fptr[0], 2, dims, &status);
        if (status) {
            fits_report_error(stderr,status);
            return status;
        }
        if (debug) fprintf(fpd,"Found %ld x %ld image\n",dims[0],dims[1]);
    }

    return EXIT_SUCCESS;
}


void allocate_arrays(size_t matLength, size_t lccspcLength, size_t lacspcLength, size_t fullLength,
                    float complex **cuda_matrix_h, float complex **full_matrix_h,
                    float complex **lccspc_h, float **lacspc_h, float complex **lcc_accum, float **lac_accum){
    *cuda_matrix_h = (float complex *) malloc(matLength * sizeof(float complex));
    *full_matrix_h = (float complex *) malloc(fullLength * sizeof(float complex));
    *lccspc_h = (float complex *) malloc(lccspcLength*sizeof(float complex));
    *lacspc_h = (float *) malloc(lacspcLength*sizeof(float));
    *lcc_accum = (float complex *) calloc(lccspcLength,sizeof(float complex));
    *lac_accum = (float *) calloc(lacspcLength,sizeof(float));
    assert(*cuda_matrix_h!=NULL && *full_matrix_h!=NULL && *lccspc_h!=NULL && *lacspc_h!=NULL);
}




int main(int argc, char **argv) {

    int arg = 0,mem_allocated=0;
    int mode = 0; // mode 0 is single timestep from a correlator dump
    // mode 1 is multiple timesteps from a FITS file
    long datadims[2];

    char *output_file = NULL,prod_type='B';
    int num_hdus_in_file[MAX_FILES];
    int ninput;
    int nfiles = 0;
    int start_sec = 0;
    int nseconds = 99999;
    int fscrunch_factor=1;
    int tscrunch_factor=1;
    int primary = 1; //image NOT in primary header?: see NGAS (1==true)
    int done_extracting=0;
// These are declared in this file and hence are not externals
//    extern int nfrequency;
//    extern int npol;
//    extern int nstation;

//    extern char *optarg;  // This is in getopt.h
    fpd = stdout;

    nstation = 0;
    nfrequency = 0;
    npol = 2;
    ninput = nstation*npol;

    // should we append or overwrite the output file
    int appendtofile=0;
    char lccfilename[1024];        
    char lacfilename[1024];

    float complex *cuda_matrix_h = NULL;
    float complex *full_matrix_h = NULL;
    float complex null = 0x0;
    float complex * lccspc_h = NULL, *lcc_accum=NULL; 
    float complex * lcc_base = NULL;
    float * lacspc_h = NULL, *lac_accum=NULL; 
    float * lac_base = NULL ;

    size_t nvis=0;
    size_t matLength=0, lccspcLength=0, lacspcLength=0, fullLength=0;

    if (argc == 1) {
        usage();
    }

    memset(num_hdus_in_file,0,sizeof(num_hdus_in_file));

    while ((arg = getopt(argc, argv, "m:i:f:F:s:v:n:o:apT:P:d")) != -1) {

        switch (arg) {
            case 'h':
                usage();
                break;
            case 'v':
                input_file[nfiles] = strdup(optarg);
                nfiles++;
                break;
            case 'm':
                mode = atoi(optarg);
                break;
            case 'i':
                ninput = atoi(optarg);
                nstation = ninput/npol;
                break;
            case 'f':
                nfrequency = atoi(optarg);
                break;
            case 'F':
                fscrunch_factor = atoi(optarg);
                break;
            case 's':
                start_sec = atoi(optarg);
                break;
            case 'n':
                nseconds = atoi(optarg);
                break;
            case 'P':
                prod_type = optarg[0];
                break;
            case 'o':
                output_file = strdup(optarg);
                break;
            case 'p':
                primary = 0;
                break;
            case 'a':
                appendtofile = 1;
                break;
            case 'T':
                tscrunch_factor = atoi(optarg);
                break;
            case 'd':
                debug++;
                break;
            default:
                usage();
                break;
        }
    }

    /* sanity checks */
    if (fscrunch_factor !=1 && fscrunch_factor !=2 && fscrunch_factor !=4 && fscrunch_factor !=8 && fscrunch_factor !=16) {
        fprintf(stderr,"freq averaging must be power of 2\n");
        exit(EXIT_FAILURE);
    }

    if (nfiles == 0) {
        mode = 0;
        input_file[0] = "/tmp/last_dump.fits";
        nfiles = 1;

    }
    if (output_file == NULL) {
        output_file = "last_dump";
    }

    if (prod_type!='B' && prod_type !='A' && prod_type!='C') {
        fprintf(stderr,"Unknown product type '%c'\n",prod_type);
        exit(EXIT_FAILURE);
    }

    int ifile,stophdu;
    int ihdu = 1+start_sec;

    if (mode == 1 || mode == 0) {
        int n_averaged=0,iter=0,i;
 
        /* new format vis file: vis are in image, not binary table */
        fitsfile *fptr[MAX_FILES];      /* FITS file pointer, defined in fitsio.h */
        int status = 0;   /*  CFITSIO status value MUST be initialized to zero!  */
        int hdutype, ncols, typecode;
        long nrows,repeat,width_in_bytes;

        /* open files */
        for (ifile = 0; ifile < nfiles; ifile++) {
            if (debug) fprintf(fpd,"Opening file %s \n",input_file[ifile]);
            if (fits_open_file(&(fptr[ifile]), input_file[ifile], READONLY, &status)!=0){
                fprintf(stderr,"Cannot open file %s\n", input_file[ifile]);
                printerror(status);
            }
        }

        if (mode == 0) primary = 0;
        ihdu += primary;    // the data starts in header ihdu.
        stophdu = ihdu + nseconds-1;
        /* find out how many HDUs are in each file if we don't want them all */
        status = getHDUsInfo(nfiles,fptr,num_hdus_in_file,datadims);
        if (status) exit(EXIT_FAILURE);

        /* check that there are the same number of HDUs in each file and only extract the amount of time
        that there is actually data for */
        for (ifile=0; ifile < nfiles; ifile++) {
            if (stophdu > num_hdus_in_file[ifile]) {
                stophdu = num_hdus_in_file[ifile];
            }
        }

        if (debug) fprintf(fpd,"Start HDU: %d, stop HDU: %d\n",ihdu,stophdu);

        /* if we're time averaging, check that an integer number of averages will go into output */
        if ((num_hdus_in_file[0]-1)%tscrunch_factor != 0) {
            fprintf(stderr,"WARNING: There are %d time steps in file and time averaging of %d steps requested.\n",num_hdus_in_file[0]-1,tscrunch_factor);
            fprintf(stderr,"\tThis will truncate data from the output\n");
        }

        /* update globals based on in-file metadata */
        /* The number of correlation products is slightly different to the standard n(n+1)/2
           formula since there is a small fraction of redundant info in the file. If n is the
           number of stations (not stations*pols), and N is the dimension of the data in the file
           then n = (-1 + sqrt(1+N))/2
        */
        nstation = (((int) round(sqrt(datadims[0]+1)))-1)/2;
        ninput = nstation*npol;
        nfrequency = datadims[1]*nfiles;
        nvis = (nstation + 1)*(nstation/2)*npol*npol;
        if (debug) {
            fprintf(fpd,"Calculated there are %d stations in the data.",nstation);
            fprintf(fpd,"\tnpol: %d. nfreq: %d, nvis: %ld\n",npol,nfrequency,(long)nvis);
        }

        matLength = nfrequency * nvis; // cuda_matrix length (autos on the diagonal)
        lccspcLength = nfrequency * (ninput-1)*ninput/2; //Lfile matrix length (no autos)
        lacspcLength = nfrequency * ninput; //Lfile autos length
        fullLength = nfrequency * ninput*ninput;

        /* now we know the dimensions of the data, allocate arrays */
        if (!mem_allocated){
            allocate_arrays(matLength, lccspcLength, lacspcLength, fullLength,
                            &cuda_matrix_h, &full_matrix_h, &lccspc_h, &lacspc_h, &lcc_accum, &lac_accum);
            mem_allocated=1;
        }

        lcc_base = lccspc_h;
        lac_base = lacspc_h;

        sprintf(lccfilename,"%s.LCCSPC",output_file);
        sprintf(lacfilename,"%s.LACSPC",output_file);

        FILE *autos=NULL;
        FILE *cross=NULL;

        if (!appendtofile) {
            if (prod_type=='B' || prod_type=='A') autos = fopen(lacfilename,"w");
            if (prod_type=='B' || prod_type=='C') cross = fopen(lccfilename,"w");
        }
        else {
            if (prod_type=='B' || prod_type=='A') autos = fopen(lacfilename,"a");
            if (prod_type=='B' || prod_type=='C') cross = fopen(lccfilename,"a");
        }

        if ((autos == NULL && (prod_type=='A'||prod_type=='B')) || (cross == NULL && (prod_type=='B'||prod_type=='C'))) {
            fprintf(stderr,"Cannot open %s or %s\n",lacfilename,lccfilename);
            exit(EXIT_FAILURE);
        }

        fill_mapping_matrix();

        printf("Extracting");

        while (!done_extracting) { // (ihdu <= stophdu)

            for (ifile = 0; ifile < nfiles; ifile++) {
                 status=0;

                /* Move to the next hdu and get the HDU type */
                if (fits_movabs_hdu(fptr[ifile], ihdu, &hdutype, &status)) {
                    if (status==107) {  /* end of file */
                        done_extracting=1;
                        break;
                    }
                    printerror(status);
                }

                if (hdutype == BINARY_TBL) {

                    if (debug) fprintf(fpd,"Detected Binary Table HDU\n");

                    fits_get_num_rows(fptr[ifile], &nrows, &status);
                    fits_get_num_cols(fptr[ifile], &ncols, &status);

                    status = 0;
                    fits_get_coltype(fptr[ifile],1,&typecode,&repeat,&width_in_bytes,&status);
                    if (debug) {
                        fprintf(fpd,"table have %ld rows and %d cols\n",nrows,ncols);
                        fprintf(fpd,"each col is type %d with %ld entries and %ld bytes long\n",typecode,repeat,width_in_bytes);
                    }

                    if (nfrequency != (nrows*nfiles)) {
                        fprintf(stderr, "nfiles (%d) * nrows (%ld) not equal to nfrequency (%d) are the FITS files the dimension you expected them to be?\n",\
                                nfiles,nrows,nfrequency);
                        exit(EXIT_FAILURE);
                    }

                    int anynull = 0x0;

                    float complex *ptr = cuda_matrix_h + (ifile * nrows * nvis);

                    status = 0;
                    fits_read_col(fptr[ifile],typecode,1,1,1,(repeat*nrows),&null,ptr,&anynull,&status);
                    if (status != 0 ) {
                        printf("Error reading columns");
                    }
                }
                else {  /* image extension */

                    status=0;
                    long fpixel = 1;
                    float nullval = 0;
                    int anynull = 0x0;
                    long naxes[2];
                    long npixels = 0;

                    if (fits_get_img_size(fptr[ifile],2,naxes,&status)) {
                        printerror(status);
                    }

                    status = 0;
                    npixels = naxes[0] * naxes[1];

                    if (nfrequency != (naxes[1]*nfiles)) {
                        fprintf(stderr, "nfiles (%d) * nrows (%ld) not equal to nfrequency (%d) are the FITS files the dimension you expected them to be?\n",\
                                nfiles,naxes[1],nfrequency);
                        exit(EXIT_FAILURE);
                    }


                    float complex *ptr = cuda_matrix_h + (ifile * naxes[1] * nvis);
                    if (fits_read_img(fptr[ifile],TFLOAT,fpixel,npixels,&nullval,(float *)ptr,&anynull,&status)) {
                        printerror(status);
                    }
                        //
                        // now need to read in the relevant part of the image
                        //
                        //

                }
            }

            if (done_extracting) continue;      // break loop if we didn't read anything

            printf("."); fflush(stdout);

            extractMatrix(full_matrix_h, cuda_matrix_h);
            /* now build lfile data */        
            int input1=0,input2=0;
            for (input1 = 0; input1 < ninput; input1++) {
                for (input2 = input1 ; input2 < ninput; input2++) {
                    map_t the_mapping = corr_mapping[input1][input2];
                    get_baseline_lu(the_mapping.stn1, the_mapping.stn2, the_mapping.pol1,
                            the_mapping.pol2, full_matrix_h, lccspc_h);
                    if (input1 == input2) {
                        /* auto */
                        int i=0;
                        for (i=0;i<nfrequency;i++) {
                            *lacspc_h = crealf(lccspc_h[i]);
                            lacspc_h++;
                        }

                    }
                    else {
                        /* cross */
                        lccspc_h += nfrequency;
                    }
                }
            }
            /* add into average */
            for (i=0; i<lccspcLength; i++) lcc_accum[i] += lcc_base[i]/tscrunch_factor;
            for (i=0; i<lacspcLength; i++) lac_accum[i] += lac_base[i]/tscrunch_factor;
            n_averaged++;

            /* if we have accumulated enough into a time average, write it out */
            if ((n_averaged % tscrunch_factor) == 0) {

                /* average over frequency, if necessary */
                if (fscrunch_factor > 1) {
                    float ac_av;
                    float complex cc_av;
                    long i,j;
                    for (i=0; i<lacspcLength/fscrunch_factor; i++) {
                        ac_av=0.0;
                        for (j=0; j<fscrunch_factor; j++) ac_av += lac_accum[i*fscrunch_factor+j];
                        lac_accum[i] = ac_av/fscrunch_factor;
                    }
                    for (i=0; i<lccspcLength/fscrunch_factor; i++) {
                        cc_av=0.0;
                        for (j=0; j<fscrunch_factor; j++) cc_av += lcc_accum[i*fscrunch_factor+j];
                        lcc_accum[i] = cc_av/fscrunch_factor;
                    }
                }

                if (autos !=NULL) fwrite(lac_accum,sizeof(float),lacspcLength/fscrunch_factor,autos);
                if (cross !=NULL) fwrite(lcc_accum,sizeof(float complex),lccspcLength/fscrunch_factor,cross);
                memset(lcc_accum,0,lccspcLength*sizeof(float complex));
                memset(lac_accum,0,lacspcLength*sizeof(float));
                if (debug) {
                    fprintf(fpd,"iter %d. Writing average\n",iter);
                }
                n_averaged=0;
            }
            lacspc_h = lac_base;
            lccspc_h = lcc_base;

            ihdu++;
            iter++;
            if (ihdu > stophdu) done_extracting=1;
        }
        status=0;
        for (ifile = 0; ifile < nfiles; ifile++) {
            fits_close_file(fptr[ifile], &status);
        }

        printf(" Done.\n");

        if (autos!=NULL) fclose(autos);
        if (cross!=NULL) fclose(cross);

    }
    else if (mode == 2) {
        /* old format vis file */
        FILE *matrix;
        int rtn;
        float complex *cuda_matrix_h = NULL;
        float complex *full_matrix_h = NULL;

        float complex * lccspc_h = NULL; 
        float complex * lcc_base = NULL;
        float * lacspc_h = NULL; 
        float * lac_base = NULL ;

        size_t matLength = nfrequency * (ninput+npol)*ninput/2; // cuda_matrix length (autos on the diagonal)
        size_t lccspcLength = nfrequency * (ninput-1)*ninput/2; //Lfile matrix length (no autos)
        size_t lacspcLength = nfrequency * ninput; //Lfile autos length
        size_t fullLength = nfrequency *ninput*ninput;

        cuda_matrix_h = (float complex *) malloc(matLength * sizeof(float complex));
        full_matrix_h = (float complex *) malloc(fullLength * sizeof(float complex));
        lccspc_h = (float complex *) malloc(lccspcLength*sizeof(float complex));
        lacspc_h = (float *) malloc(lacspcLength*sizeof(float));

        lcc_base = lccspc_h;
        lac_base = lacspc_h;

        matrix = fopen(input_file[0], "r");

        if (matrix == NULL) {
            perror("On opening file");
            exit(EXIT_FAILURE);
        }

        rtn = fread(cuda_matrix_h, sizeof(float complex), matLength, matrix);

        if (rtn != matLength) {
            if (!feof(matrix)) {
                perror("On reading file");
            } else {

                fprintf(stderr,
                        "EOF before full matrix read (partial read of %d elemetns)\n",
                        rtn);
                exit(EXIT_FAILURE);
            }
        } else {


            fprintf(stdout,
                    "Full matrix read in for a single time step  -- will reorder to triangular\n");
            // wacky packed tile order to packed triangular
            // xgpuReorderMatrix((Complex *) cuda_matrix_h);
            // convert from packed triangular to full matrix
            extractMatrix_slow(full_matrix_h, cuda_matrix_h);
            // get the mapping
            fill_mapping_matrix();
        }

        /* now build lfile data */
        int input1=0,input2=0;
        for (input1 = 0; input1 < ninput; input1++) {
            for (input2 = input1 ; input2 < ninput; input2++) {
                map_t the_mapping = corr_mapping[input1][input2];
                get_baseline_lu(the_mapping.stn1, the_mapping.stn2, the_mapping.pol1,
                        the_mapping.pol2, full_matrix_h, lccspc_h);
                if (input1 == input2) {
                    /* auto */
                    int i=0;
                    for (i=0;i<nfrequency;i++) {
                        *lacspc_h = crealf(lccspc_h[i]);
                        lacspc_h++;
                    }

                }
                else {
                    /* cross */
                    lccspc_h += nfrequency;
                }
            }

        }
        FILE *autos=NULL;
        FILE *cross=NULL;
        autos = fopen("last_dump.LACSPC","w");
        fwrite(lac_base,sizeof(float),lacspcLength,autos);
        fclose(autos);
        cross = fopen ("last_dump.LCCSPC","w");
        fwrite(lcc_base,sizeof(float complex),lccspcLength,cross);
        fclose(cross);

        fclose(matrix);
        exit(EXIT_SUCCESS);

    }

    free(lcc_base);
    free(lac_base);

    return 0;
}
