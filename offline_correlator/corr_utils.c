#include "fitsio.h"
#include "corr_utils.h"

//PJE removed common includes to header 

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

inline void integrate(char * output,uint8_t *input, int factor, int nint,int nfreq) {

    /*      integrates by adding samples together -- assumes single pol
     * 	for multiple pol it dependes whether we are packing the pols together - so this would affect both the input and output strides
     */

    int integration = 0;
    uint16_t *output_samp = (uint16_t *) output;

    for (integration=0;integration<nint;integration++) {

        for (int chan=0;chan<nfreq;chan++){
            uint8_t *inp_ptr = input + (integration*factor*nfreq);
            for (int sample=0;sample<factor;sample++) {
                if (sample == 0) {
                    output_samp[integration*nfreq+chan] = *inp_ptr;
                }
                else {
                    output_samp[integration*nfreq+chan] += inp_ptr[sample*nfreq];
                }
            }
        }
    }

}


void buildFITSBuffer( XGPUInfo xgpu_info,
                      Complex *full_matrix_h,
                      size_t blockSize,
                      void *out,
                      time_t current_time_t,
                      int dumps_per_second,
                      manager_t *manager )

/* The purpose of this is to re-order the output array and build the binary
 * FITS table. Somethings only have to be done once (like the mapping) other
 * things have to be done every time.
 *
 * The nice thing about this is we will have a full self describing output
 * file that can be read by any fits reader.
 *
 * The format is a FITS image extension nbaselines (cols) by nchan (rows)
 * This ordering allows easy appending of multi-channel data.
 * The FITS file is created in memory and not written out.
 */

{
    int chans_to_aver = manager->chan_to_aver;

    int marker = manager->marker;

	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
	int status = 0;

	int nfrequency = xgpu_info.nfrequency/chans_to_aver;
	size_t newsize = 2*blockSize;
	int frac_time = (1000/dumps_per_second);
	float int_time = 1.0/dumps_per_second;

	int milli_time = 0;
	int num_baselines = xgpu_info.nbaseline*4; // num baselines x npol x npol
	int num_floats_per_chan = num_baselines * 2; // complex
	void *memblock = NULL;


	/* we need to define a block of memory that is big enough to hold the whole FITS file + headers and padding
	 * this also has to be precisely the size of the memory block we hand over the NGASS so it knows where and how to trim
	 *
	 * I'm going to allocate a buffer here as it seems that the FITS API wants to mess with this
	 */
	memblock = calloc(2*blockSize,1);
	assert(memblock);
	/* create a FITS file that will reside in memory */
	if (fits_create_memfile(&fptr, (void **) &memblock, &newsize,2880,realloc,&status)) {
		printerror(status);
	}


	/* need to create a new binary table and the associated HDU */

	char extname[128];

	sprintf(extname,"%ld",current_time_t);


	int bitpix = LONG_IMG;
	long naxis = 2; // 2 dimensions
	long naxes[2];

	naxes[0] = num_floats_per_chan; // pixels wide (number of baselines x 2 as we are writing as FLOATS not COMPLEX
	naxes[1] = nfrequency;
	int dump = 0;

	float *average = (float *) calloc(num_floats_per_chan,sizeof(float)); // i'm going to average into here
	if (average == NULL) {
		syslog(LOG_CRIT,"FITSBuilder::Failure to allocate memory\n");
		exit(EXIT_FAILURE);
	}



    while (dump < dumps_per_second) {



        status = 0;

        /* append a new empty FITS IMAGE onto the FITS file */
        if (fits_create_img( fptr, bitpix, naxis, naxes,
                            &status) ) {
            printerror( status );
            syslog(LOG_CRIT,"Failure to initialize FITS image\n");
            exit(EXIT_FAILURE);
        }
        /* Now we need to re-order the gpu output as will be required by the
         * astronomers
         */
        // wacky packed tile order to packed triangular
        Complex *ptr = full_matrix_h + (dump * xgpu_info.matLength);


        xgpuReorderMatrix((Complex *) ptr);

        float *buff = (float *) ptr;

        long nelements=0,fpixel[2];


        nelements = naxes[0]*naxes[1];

        fpixel[0] = 1;
        fpixel[1] = 1;

        int f=0;
        int index=0;

        do {

            for (f=0;f<chans_to_aver;f++) {

                for (index=0;index<num_floats_per_chan;index++) {

                    average[index] = average[index] + (*buff);
                    buff++;

                }

            }

            float *tofits = (float *) average;
            status = 0;
            nelements = num_floats_per_chan;

            if (fits_write_pix(fptr,TFLOAT,fpixel,nelements,tofits,&status)) {
                printerror(status);
                syslog(LOG_CRIT,"Failure to write FITS image\n");
                exit(EXIT_FAILURE);

            }

            fpixel[1]++; // next channel (I hope)

            bzero((void *)average,num_floats_per_chan*sizeof(float));

        } while (fpixel[1] <= naxes[1]);


        status = 0;
        if (fits_update_key(fptr,TLONG,(char *) "TIME",&current_time_t,(char *) "Unix time (seconds)",&status)) {
            printerror(status);
            syslog(LOG_CRIT,"Failed to add time_t to the primary fits header\n");
            exit(EXIT_FAILURE);
        }


        int_time = 1.0/dumps_per_second;
        milli_time = dump*frac_time;


        if (fits_update_key(fptr,TINT,(char *) "MILLITIM",&milli_time,(char *) "Milliseconds since TIME",&status)) {
            printerror(status);
            syslog(LOG_CRIT,"Failed to add milliseconds to the primary fits header\n");
            exit(EXIT_FAILURE);
        }


        if (fits_update_key(fptr,TFLOAT,(char *) "INTTIME",&int_time,(char *) "Integration time (s)",&status)) {
            printerror(status);
            syslog(LOG_CRIT,"Failed to add integer seconds to the primary fits header\n");
            exit(EXIT_FAILURE);
        }

        if (fits_update_key(fptr,TINT,(char *) "MARKER",&marker,(char *) "Data offset marker (all channels should match)",&status)) {
            printerror(status);
            syslog(LOG_CRIT,"Failed to add alignment marker to fits header\n");
            exit(EXIT_FAILURE);
        }

        if (manager->coarse_chan >=0)
            if (fits_update_key(fptr,TINT,(char *) "COARSE_CHAN",&manager->coarse_chan,(char *) "Receiver Coarse Channel Number (only used in offline mode)",&status)) {
                printerror(status);
                syslog(LOG_CRIT,"Failed to add coarse channel number to header\n");
                exit(EXIT_FAILURE);
            }

        float bscale_parameter = 1.0/(chans_to_aver*int_time);

        if (fits_update_key(fptr,TFLOAT,(char *) "BSCALE",&bscale_parameter,(char *) "Incorporates channel and time averaging",&status)) {
            printerror(status);
            syslog(LOG_CRIT,"Failed to add time_t to the primary fits header\n");
            exit(EXIT_FAILURE);
        }

        dump++;
    }
	status = 0;
	if (fits_close_file(fptr,&status)) {
		printerror(status);
		syslog(LOG_CRIT,"Failed to close FITS file in memory\n");
		exit(EXIT_FAILURE);
	}

	/* now I am going to copy to file over to the buffer
	 *
	 */
	if (newsize != 2*blockSize) {
		fprintf(stderr,"FITS API has resized the memory area to %lu from %lu",newsize, 2*blockSize);
	}
	memcpy(out,memblock,blockSize);

	free(memblock);
	free(average);



	/* done */

}
