#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "psrfits.h"

void usage()
{
    printf( "splice_psrfits from VCSTools %s\n", VERSION_BEAMFORMER );
    printf( "usage: splice_psrfits <file 1> <file 2> ..... <basefilename>\n\n" );
    printf( "Utility to splice psrfits files together providing they have the same start time\n" );
    printf( "The output will get written to \"[basefilename]_0001.fits\". This file must not already exist.\n" );
}

int prep_data(void *to_load)
{
    struct psrfits *pf = (struct psrfits *) to_load;
    int rv = 0;
    pf->sub.dat_freqs = (float *)malloc(sizeof(float) * pf->hdr.nchan);
    pf->sub.dat_weights = (float *)malloc(sizeof(float) * pf->hdr.nchan);
    pf->sub.dat_offsets = (float *)malloc(sizeof(float) * pf->hdr.nchan * pf->hdr.npol);
    pf->sub.dat_scales  = (float *)malloc(sizeof(float) * pf->hdr.nchan * pf->hdr.npol);
    pf->sub.data  = (unsigned char *)malloc(pf->sub.bytes_per_subint);
    pf->sub.rawdata  = (unsigned char *)malloc(pf->sub.bytes_per_subint);
    rv = 1;
    return rv;
}

int cleanup(void *to_clean)
{
    struct psrfits *pf = (struct psrfits *) to_clean;
    free(pf->sub.dat_freqs);
    free(pf->sub.dat_weights);
    free(pf->sub.dat_offsets);
    free(pf->sub.dat_scales);
    free(pf->sub.data);
    free(pf->sub.rawdata);
    return 1;
}

float round_three(float var)
{
    // round to the nearest three decimal places
    float value = (int)(var * 1000 + .5);
    return (float)value / 1000;
}

int append(char **to_append, void *total, int n){

    // this would perhaps be better as C++
    // this method appends a number of psrfits files to the total

    // sets the total pointer to the input argument
    // this has already been populated by a open/close pair

    struct psrfits *pf_total = (struct psrfits *) total;


    int ii=0,jj=0,f=0,rv=0,sub=0;

    /* first check the channels */
    /* last channel of the total */

    struct psrfits pf1,pf2;
    struct psrfits *list;
    struct psrfits *pf = NULL;

    // list of psrfits structures
    list = (struct psrfits *) calloc(n,sizeof(struct psrfits));

    // zero the fitst struct
    bzero(&pf1,sizeof(struct psrfits));

    // set to the first filename
    pf1.filenames = (char **) malloc (sizeof(char *));
    pf1.filenames[0] = to_append[0];

    pf1.status=0;
    pf1.filenum=0;
    pf1.numfiles=1; // using explicit filenames

    rv = psrfits_open(&pf1);
    if (rv) { fits_report_error(stderr, rv); exit(-1);}

    // the above open command fills the struct - so we now know something about the file

    // allocate memory for 1 based upon the contents of the struct
    prep_data(&pf1);
    // read a subint
    psrfits_read_subint(&pf1);
    // close thge file
    rv = psrfits_close(&pf1);
    // zero the 2nd struct
    bzero(&pf2,sizeof(struct psrfits));

    pf2.filenames = (char **) malloc (sizeof(char *));
    //set to the name of file 2

    pf2.filenames[0] = to_append[1];

    pf2.status=0;
    pf2.filenum=0;
    pf2.numfiles=1; // using explicit filenames

    // this opne populates the struct
    psrfits_open(&pf2);
    // allocate memory
    prep_data(&pf2);
    psrfits_read_subint(&pf2);
    // read a subint
    // close the file
    rv = psrfits_close(&pf2);

    if (rv) { fits_report_error(stderr, rv); exit(-1);}

    // now do some bookkeeping to ensure we can actually append these files

    // what is the last channel of the first file in our list
    float last_chan_of_total = pf1.sub.dat_freqs[pf1.hdr.nchan-1];
    // what is the last channel of the second
    float first_chan_of_to_append = pf2.sub.dat_freqs[0];
    // what is deltaf
    float df = pf1.sub.dat_freqs[1] - pf1.sub.dat_freqs[0];
    // what is the freq diff between the files
    float diff_freq = first_chan_of_to_append - last_chan_of_total;
    int chans_to_pad = 0;

    if (round_three(diff_freq) != round_three(df)){
        fprintf( stderr, "warning: channels are not contiguous; attempting "
                         "to pad\n" );
        fprintf( stderr, "         trying to append channel %f to channel %f.\n",
                 first_chan_of_to_append,
                 last_chan_of_total );

        chans_to_pad = rint((diff_freq-df)/df);

        if (chans_to_pad < 0)
        {
            fprintf( stderr, "error: %d pads required -- channel order "
                             "switched - check file order\n", chans_to_pad );
            return(-1);
        }
        else
        {
            fprintf( stderr, "%d pads required\n", chans_to_pad );
        }
    }

    fprintf(stdout,"Done checking channels\n");

    pf_total->rownum=1; // set to the beginning
    pf_total->filenum=0;
     // set the nchan of the total - remember this has already been partially filled by
    // an open / close pair
    // total number of channels is n* the number we already have plus the pads
    pf_total->hdr.nchan = (pf1.hdr.nchan + chans_to_pad)*n; // dont forget the two half pads at the edges
    // centre freq - this is the original base frequency + half the channels
    // this seems to assuming that fctr is the base of central channel - dont forget that we are going to pad
    pf_total->hdr.fctr = (pf1.hdr.fctr - pf1.hdr.BW/2.0) - (chans_to_pad*0.5)*pf_total->hdr.df + (pf_total->hdr.nchan*0.5)*pf_total->hdr.df;
    // BW number of channels time their width
    pf_total->hdr.BW = pf_total->hdr.nchan * pf_total->hdr.df;
    // bytes per subint
    pf_total->sub.bytes_per_subint = (pf_total->hdr.nbits * pf_total->hdr.nchan * pf_total->hdr.npol * pf_total->hdr.nsblk) / 8;


    /* empty the arrays from the earlier channel checks */

    cleanup(&pf2);
    cleanup(&pf1);

    // zero the first struct
    bzero(&pf1,sizeof(struct psrfits));
    // zero the second struct
    bzero(&pf2,sizeof(struct psrfits));


    fprintf(stdout,"Allocating buffers\n");


    pf_total->sub.dat_freqs = (float *)malloc(sizeof(float) * pf_total->hdr.nchan);
    pf_total->sub.dat_weights = (float *)malloc(sizeof(float) * pf_total->hdr.nchan);
    // this changes the frequency labling for the output
    // this is consistent with fctr being the base of the central channel and
    // each channel label being the centre of the channel
    float dtmp = pf_total->hdr.fctr - 0.5 * pf_total->hdr.BW + 0.5 * pf_total->hdr.df;
    for (ii = 0 ; ii < pf_total->hdr.nchan ; ii++) {
        pf_total->sub.dat_freqs[ii] = dtmp + ii * pf_total->hdr.df;
        pf_total->sub.dat_weights[ii] = 1.0;
    }
    pf_total->sub.dat_offsets = (float *)malloc(sizeof(float) * pf_total->hdr.nchan * pf_total->hdr.npol); // definitely needed for 8 bit numbers
    pf_total->sub.dat_scales = (float *)malloc(sizeof(float) * pf_total->hdr.nchan * pf_total->hdr.npol);

    pf_total->sub.data = (unsigned char *)malloc(pf_total->sub.bytes_per_subint);
    pf_total->sub.rawdata = pf_total->sub.data;


    printf( "Build combined dat_weights array\n" );
    printf( "Build combined scales and offsets\n" );

    for (ii = 0 ; ii < pf_total->hdr.nchan * pf_total->hdr.npol ; ii++) {

        pf_total->sub.dat_offsets[ii] = 0.0;
        pf_total->sub.dat_scales[ii] = 1.0;

    }


    pf_total->filenames = (char **)malloc(sizeof(char *));
    pf_total->filenames[0] = "total.fits";

    for (f=0;f<n;f++) { // for each file in the list

        pf = &list[f]; // get the struct

        pf->filenames = (char **) malloc (sizeof(char *));
        pf->filenames[0] = to_append[f]; // get the filename

        pf->status=0;
        pf->filenum=0;
        pf->numfiles=1; // using explicit filenames

        //fprintf(stdout,"Opening file (%s) %d of %d\n",pf->filenames[0],f+1,n);
        // populate the struct
        psrfits_open(pf);
        // allocate the arrays
        prep_data(pf);
    }
    // now we have to transpose a bit
    // we need to step through every constituent file for every subint

    int nsub = 1000;
    for (sub = 0; sub < nsub; sub++) { // foreach subint

        size_t start_offset = 0;
        // i have to interlace the channels unfortunately
        // we need full total_nchan of each polarisation - not
        // nchan*npol blocks
        // loop over npol

        printf("Read subint %d (file 0-%d, pol 1-4, row %d/%d)\n", sub,
                    n, pf->rownum-1, pf->rows_per_file);
        for (f=0;f<n;f++) {

            pf = &list[f]; // get the psrfits struct - it is already open

            rv = psrfits_read_subint(pf); // read a subint into it

            if (rv) { fits_report_error(stderr, rv); return rv;}
            int p=0; // the current Polarisation
            for (p=0;p<pf_total->hdr.npol;p++) {
                // ok so I now increment the start offset by the file number
                start_offset = p*pf_total->hdr.nchan + f*(pf->hdr.nchan + chans_to_pad);

                // where do i put this file in the list of offsets and scales
                // start at the freq offset which is npad*filenum*npol
                // then add half the pad
                // keep going for the nchan samples for this pol
                // We start at half the pad size in:
                // start_offset + edge_pad*0.5*npol
                // which means we stop at
                //      start_offset + edge_pad*0.5*npol + nchan*npol
                // The orders are also interleaved - nchan of each pol - which
                // means we have to split out the offsets - unless they are the same for all
                // pols.
                //
                int ii_start = start_offset+(chans_to_pad/2);
                int ii_stop = ii_start + pf->hdr.nchan;

                for (ii = ii_start, jj=0 ; ii < ii_stop ; ii++, jj++) {
                // our index is
                    // 0 to nchan
                    pf_total->sub.dat_offsets[ii] = pf->sub.dat_offsets[jj];
                    pf_total->sub.dat_scales[ii] = pf->sub.dat_scales[jj];

                    //fprintf(stdout,"output scales start index: %d\n",ii);
                }

                /* now we need to append the arrays */

                size_t bytes_transferred = 0;
                size_t bytes_to_copy = 0;
                size_t bytes_per_subint_per_pol = pf->sub.bytes_per_subint/pf->hdr.npol;

                unsigned char *output_pos_init = pf_total->sub.rawdata + start_offset;
                unsigned char *pf_input_pos_init = pf->sub.rawdata + p*pf->hdr.nchan;


                int timestep = 0;
                size_t output_step_between_row = pf_total->hdr.nchan*(pf_total->hdr.npol);
                size_t input_step_between_row = pf->hdr.nchan*(pf->hdr.npol);

                while (bytes_transferred < bytes_per_subint_per_pol){
                //    printf("step %d bytes: %lu\r",timestep,bytes_transferred);
                    unsigned char *output_pos = output_pos_init+(timestep*output_step_between_row);
                    unsigned char *pf_input_pos = pf_input_pos_init+(timestep*input_step_between_row);
                // memcpy the pad and set to zero

                    bytes_to_copy = (chans_to_pad/2);
                    bzero(output_pos,bytes_to_copy);
                    output_pos += bytes_to_copy;

                //memcpy the channels

                    bytes_to_copy = pf->hdr.nchan;

                    memcpy(output_pos,pf_input_pos,bytes_to_copy);
                    bytes_transferred += bytes_to_copy;
                    output_pos += bytes_to_copy;

                // memcpy the pad

                    bytes_to_copy = chans_to_pad/2;
                    bzero(output_pos,bytes_to_copy);
                    output_pos += bytes_to_copy;

                    timestep++;


                }
            }
        }
	    // Now we need to get the MetaData for this subint
	    pf_total->sub.feed_ang = pf->sub.feed_ang; // Feed angle at subint centre (deg)
	    pf_total->sub.tsubint  = pf->sub.tsubint;  // Length of subintegration (sec)
	    pf_total->sub.offs     = pf->sub.offs;     // Offset from Start of subint centre (sec)
	    pf_total->sub.lst      = pf->sub.lst;      // LST at subint centre (sec)
	    pf_total->sub.ra       = pf->sub.ra;       // RA (J2000) at subint centre (deg)
	    pf_total->sub.dec      = pf->sub.dec;      // Dec (J2000) at subint centre (deg)
	    pf_total->sub.glon     = pf->sub.glon;     // Gal longitude at subint centre (deg)
	    pf_total->sub.glat     = pf->sub.glat;     // Gal latitude at subint centre (deg)
	    pf_total->sub.par_ang  = pf->sub.par_ang;  // Parallactic angle at subint centre (deg)
	    pf_total->sub.tel_az   = pf->sub.tel_az;   // Telescope azimuth at subint centre (deg)
	    pf_total->sub.pos_ang  = pf->sub.pos_ang;  // Position angle of feed at subint centre (deg)
	    pf_total->sub.tel_zen  = pf->sub.tel_zen;  // Telescope zenith angle at subint centre (deg)

        psrfits_write_subint(pf_total);
    }
    return 1;
}



int main(int argc, char *argv[]) {

    int i=0;

    // If "-V" appears in any of the arguments, give the version number
    // and exit
    int a;
    for (a = 0; a < argc; a++)
    {
        if (strcmp( argv[a], "-V" ) == 0)
        {
            printf( "splice_psrfits from VCS Tools %s\n",
                    VERSION_BEAMFORMER );
            exit(EXIT_SUCCESS);
        }
    }

    // Otherwise there must be at least 2 arguments
    if (argc < 3)
    {
        usage();
        exit(EXIT_FAILURE);
    }
    else
    {
        printf( "There are %d files to append\n", argc-2 );
    }

    // output file psrfits structure
    struct psrfits pf_total;

    // space for the input filenames and output filename
    char ** filenames = (char **) malloc ((argc-2) * sizeof(char *));

    // copy the file name strings

    for (i=1;i < argc-1;i++) {
        filenames[i-1] = strdup(argv[i]);
    }

    // Copy the output filename. If exists, abort.
    char * basefilename = strdup(argv[argc-1]);
    char outfile[200];
    sprintf(outfile, "%s_0001.fits", basefilename);
    FILE * foutfile;

    if ((foutfile = fopen(outfile, "r")))
    {
        fclose( foutfile );
        fprintf( stderr, "warning: File %s already exists. Overwriting.\n", basefilename);
        remove( outfile );
    }

    // this sets up the filename of the output - note initially this sets the total to be
    // the first file in the list
    pf_total.filenames = (char **) malloc (sizeof(char *));
    pf_total.filenames[0] = filenames[0];

    pf_total.status=0;
    pf_total.filenum=0;
    pf_total.numfiles=1; // using explicit filenames

    int rv = 0;

    //opening the file (the first one in the list)
    // this populates the psrfits structure

    rv = psrfits_open(&pf_total);
    if (rv) { fits_report_error(stderr, rv); exit(-1);}
    rv = psrfits_close(&pf_total);
    if (rv) { fits_report_error(stderr, rv); exit(-1);}

    sprintf(pf_total.basefilename, "%s", basefilename);
    printf( "Will append %d files to %s\n", argc-2, outfile );


    rv = append(filenames,(void *)&pf_total,argc-2);
    // append reopens the file so we must close it again
    if (rv)
    {
        psrfits_close(&pf_total);
    }

    printf( "Done\n" );

}
