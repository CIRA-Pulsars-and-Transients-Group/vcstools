#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <assert.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include "psrfits.h"
#include "star/pal.h"
#include "star/palmac.h"

size_t readn( int inSock, void *outBuf, size_t inLen ) {
    size_t  nleft;
    int nread;
    char *ptr;

    assert( inSock >= 0 );
    assert( outBuf != NULL );
    assert( inLen > 0 );

    ptr   = (char*) outBuf;
    nleft = inLen;

    while ( nleft > 0 ) {
        nread = read( inSock, ptr, nleft );
        if ( nread < 0 ) {
            if ( errno == EINTR )
                nread = 0;  /* interupted, call read again */
            else
                return -1;  /* error */
        } else if ( nread == 0 )
            break;        /* EOF */

        nleft -= nread;
        ptr   += nread;
    }

    return(inLen - nleft);
} /* end readn */

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

int main() {
    int ii;
    double dtmp;
    struct psrfits pf;
    char user_input[256];
    char *rval = 0;
    char input_pipe[256];
    int nbit_input=1;
    // Only set the basefilename and not "filename"
    // Also, fptr will be set by psrfits_create_searchmode()
    
    pf.filenum = 0;             // This is the crucial one to set to initialize things
    pf.rows_per_file = 200;     // I assume this is a max subint issue
				// Need to set this based on PSRFITS_MAXFILELEN

    fprintf(stdout,"==================================================================================\n");
    fprintf(stdout,"|| MWA-PSRFITS Generator (v0.1) - based on the Ransom and Demorest psrfits_utils ||\n");
    fprintf(stdout,"==================================================================================\n");

    strcpy(input_pipe, "/tmp/mk_psrfits_in");

    fprintf(stdout,"Name of input data pipe [%s]:",input_pipe);

    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strlen(rval) != 1) {
	    strcpy(input_pipe, strtok(rval,"\n"));
    }

    if(mkfifo(input_pipe,0666) != 0) {
	    perror("Failed to create  data pipe");
	    if (errno == EEXIST) {
		    fprintf(stderr,"pipe %s already exists ... choose a unique name or delete this pipe (if you are sure noone else is using it).",input_pipe);
		    exit(0);
	    }	
    }
    else {
	    fprintf(stdout,"Streaming input ... will use the pipe (%s) as input. Please cat the input data to this pipe\n",input_pipe);

    }
    fprintf(stdout,"Input format (1=8b) [%d]:",nbit_input);  
  
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    nbit_input = atoi(rval);

    }

    if (nbit_input != 1) {
	fprintf(stderr,"Invalid nbit input\n");
	exit(0);
    }   
    
    strcpy(pf.basefilename, "test_psrfits");

    fprintf(stdout,"Basename of output file [%s]:",pf.basefilename);

    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strlen(rval) != 1) {
	    strcpy(pf.basefilename, strtok(rval,"\n"));
    }

    // Now set values for our hdrinfo structure
    strcpy(pf.hdr.obs_mode,"SEARCH");
    fprintf(stdout,"OBS_MODE [%s]:",pf.hdr.obs_mode);

    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strlen(rval) != 1) {
	    strcpy(pf.hdr.obs_mode, strtok(rval,"\n"));
    }


    pf.hdr.scanlen = 5.0; // in sec

    fprintf(stdout,"Length of scan in this file [%lf]:",pf.hdr.scanlen);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.scanlen = atof(rval);
    }

   
    strcpy(pf.hdr.observer, "MWA User");

    fprintf(stdout,"Observer [%s]:",pf.hdr.observer);

    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strlen(rval) != 1) {
	    strcpy(pf.hdr.observer, strtok(rval,"\n"));
    }
 
    strcpy(pf.hdr.telescope, "MWA");

    fprintf(stdout,"Telescope [%s]:",pf.hdr.telescope);

    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strlen(rval) != 1) {
	    strcpy(pf.hdr.telescope, strtok(rval,"\n"));
    }


    strcpy(pf.hdr.source, "<vcs_obs>");

    fprintf(stdout,"Source [%s]:",pf.hdr.source);

    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strlen(rval) != 1) {
	    strcpy(pf.hdr.source, strtok(rval,"\n"));
    }


    strcpy(pf.hdr.frontend, "MWA-RECVR");
    fprintf(stdout,"Front End [%s]:",pf.hdr.frontend);

    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strlen(rval) != 1) {
	    strcpy(pf.hdr.frontend, strtok(rval,"\n"));
    }

    strcpy(pf.hdr.backend, "MWA-VCS");
    fprintf(stdout,"Back End [%s]:",pf.hdr.backend);

    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strlen(rval) != 1) {
	    strcpy(pf.hdr.backend, strtok(rval,"\n"));
    }


    strcpy(pf.hdr.project_id, "MWA-D0004");
    fprintf(stdout,"project_id [%s]:",pf.hdr.project_id);

    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strlen(rval) != 1) {
	    strcpy(pf.hdr.project_id, strtok(rval,"\n"));
    }


    strcpy(pf.hdr.date_obs, "2010-01-01T05:15:30.000");
    fprintf(stdout,"Date Obs [%s]:",pf.hdr.date_obs);

    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strlen(rval) != 1) {
	    strcpy(pf.hdr.date_obs, strtok(rval,"\n"));
    }


    strcpy(pf.hdr.poln_type, "LIN");
    strcpy(pf.hdr.track_mode, "TRACK");
    strcpy(pf.hdr.cal_mode, "OFF");
    strcpy(pf.hdr.feed_mode, "FA");

    pf.hdr.dt = 0.00010;			// sample rate (s)

    fprintf(stdout,"Sample Rate (s) [%lf]:",pf.hdr.dt);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.dt = atof(rval);
    }
   	
    pf.hdr.fctr = 150.0;			// frequency (MHz)

    fprintf(stdout,"Centre Frequency (MHz) [%lf]:",pf.hdr.fctr);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.fctr = atof(rval);
    }

    pf.hdr.BW = 1.28;

    fprintf(stdout,"Bandwidth (MHz) [%lf]:",pf.hdr.BW);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.BW = atof(rval);
    }


    pf.hdr.ra2000 = 302.0876876;

    fprintf(stdout,"RA (2000) (deg) [%lf]:",pf.hdr.ra2000);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.ra2000 = atof(rval);
    }

    dec2hms(pf.hdr.ra_str, pf.hdr.ra2000/15.0, 0);
    
    pf.hdr.dec2000 = -3.456987698;
 
    fprintf(stdout,"Dec (2000) (deg) [%lf]:",pf.hdr.dec2000);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.dec2000 = atof(rval);
    }

    dec2hms(pf.hdr.dec_str, pf.hdr.dec2000, 1);
    
    pf.hdr.azimuth = 123.123;
  
    fprintf(stdout,"Azimuth (deg) [%lf]:",pf.hdr.azimuth);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.azimuth = atof(rval);
    }   

    pf.hdr.zenith_ang = 23.0;
  
    fprintf(stdout,"Zenith Angle (deg) [%lf]:",pf.hdr.zenith_ang);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.zenith_ang = atof(rval);
    }   


    pf.hdr.beam_FWHM = 0.25;

    fprintf(stdout,"Beam FWHM (deg) [%lf]:",pf.hdr.beam_FWHM);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.zenith_ang = atof(rval);
    }   



    pf.hdr.start_lst = 10000.0;
    

    fprintf(stdout,"Seconds past 00h LST [%lf]:",pf.hdr.start_lst);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.start_lst = atof(rval);
    }   

    pf.hdr.start_sec = 25000.82736876;
 
    fprintf(stdout,"Seconds past 00h UTC [%lf]:",pf.hdr.start_sec);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.start_sec = atof(rval);
    }   


    pf.hdr.start_day = 55000; 
 
    fprintf(stdout,"Start MJD (whole day) [%d]:",pf.hdr.start_day);
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.start_day = atoi(rval);
    }   


    pf.hdr.scan_number = 3;
     fprintf(stdout,"Scan Number [%d]:",pf.hdr.scan_number);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.scan_number = atoi(rval);
    }   

    pf.hdr.rcvr_polns = 2;
    fprintf(stdout,"Receiver Polarisation [%d]:",pf.hdr.rcvr_polns);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.rcvr_polns = atoi(rval);
    }   



    pf.hdr.summed_polns = 1;
    fprintf(stdout,"Summed Polarisations? [1/0]  [%d]:",pf.hdr.summed_polns);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.summed_polns = atoi(rval);
    }   


    pf.hdr.offset_subint = 0;
    fprintf(stdout,"Offset Subint   [%d]:",pf.hdr.offset_subint);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.summed_polns = atoi(rval);
    }   

    pf.hdr.nchan = (pf.hdr.BW*100);
    fprintf(stdout,"Number of Channels  [%d]:",pf.hdr.nchan);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.nchan = atoi(rval);
    }   


    pf.hdr.orig_nchan = pf.hdr.nchan;
    pf.hdr.orig_df = pf.hdr.df = pf.hdr.BW / pf.hdr.nchan;


    pf.hdr.nbits = 8;

    fprintf(stdout,"Number of bits per sample  [%d]:",pf.hdr.nbits);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.nbits = atoi(rval);
    }   

    if (pf.hdr.summed_polns) {
	    pf.hdr.npol = 1;
    }
    else {
	    pf.hdr.npol = pf.hdr.rcvr_polns;
    }

    pf.hdr.nsblk = 10000;

    fprintf(stdout,"Number of spectra per row  [%d]:",pf.hdr.nsblk);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    pf.hdr.nsblk = atoi(rval);
    }   

   
    pf.hdr.MJD_epoch = pf.hdr.start_day + (pf.hdr.start_sec/86400.0);  
 
    fprintf(stdout,"Start MJD  [%Lf]:",pf.hdr.MJD_epoch);
 
    rval = fgets(user_input,256,stdin);

    if (rval == NULL) {
	perror("Failed to parse");
	exit(0);
    } 

    if (strncmp(rval,"\n",1) != 0) {
	    double temp = atof(rval);
	    pf.hdr.MJD_epoch = (long double) temp;
    }   

    pf.hdr.ds_freq_fact = 1;
    pf.hdr.ds_time_fact = 1;
    

    // some things that we are unlikely to change
    //
    pf.hdr.fd_hand = 0;
    pf.hdr.fd_sang = 0.0;
    pf.hdr.fd_xyph = 0.0;
    pf.hdr.be_phase = 0.0;
    pf.hdr.chan_dm = 0.0;

 
    // Now set values for our subint structure
    pf.tot_rows = 0;
    pf.sub.tsubint = pf.hdr.nsblk * pf.hdr.dt;
    pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;
    pf.sub.lst = pf.hdr.start_lst;
    pf.sub.ra = pf.hdr.ra2000;
    pf.sub.dec = pf.hdr.dec2000;
    palEqgal(pf.hdr.ra2000*PAL__DD2R, pf.hdr.dec2000*PAL__DD2R, 
             &pf.sub.glon, &pf.sub.glat);
    pf.sub.glon *= PAL__DR2D;
    pf.sub.glat *= PAL__DR2D;
    pf.sub.feed_ang = 0.0;
    pf.sub.pos_ang = 0.0;
    pf.sub.par_ang = 0.0;
    pf.sub.tel_az = pf.hdr.azimuth;
    pf.sub.tel_zen = pf.hdr.zenith_ang;
    pf.sub.bytes_per_subint = (pf.hdr.nbits * pf.hdr.nchan * 
                               pf.hdr.npol * pf.hdr.nsblk) / 8;
    pf.sub.FITS_typecode = TBYTE;  // 11 = byte

    // Create and initialize the subint arrays
    pf.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
    pf.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
    dtmp = pf.hdr.fctr - 0.5 * pf.hdr.BW + 0.5 * pf.hdr.df;
    for (ii = 0 ; ii < pf.hdr.nchan ; ii++) {
        pf.sub.dat_freqs[ii] = dtmp + ii * pf.hdr.df;
        pf.sub.dat_weights[ii] = 1.0;
    }
    pf.sub.dat_offsets = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
    pf.sub.dat_scales = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
    for (ii = 0 ; ii < pf.hdr.nchan * pf.hdr.npol ; ii++) {
        pf.sub.dat_offsets[ii] = 0.0;
        pf.sub.dat_scales[ii] = 1.0;
    }
 

    pf.sub.data = (unsigned char *)malloc(pf.sub.bytes_per_subint);
    pf.sub.rawdata = pf.sub.data;
    // This is what you would update for each time sample (likely just
    // adjusting the pointer to point to your data)
    //
    // We are using a pipe open it
    int fp = open(input_pipe, O_RDONLY);
 
    if(fp == -1) {
        perror("Could not open the pipe\n");
	exit(0);
    }
    pf.status=0;
    // Here is the real data-writing loop
    do {
        // Update the pf.sub entries here for each subint
        // as well as the pf.sub.data pointer
        //
	    size_t readval = readn(fp, pf.sub.data, pf.sub.bytes_per_subint);
	    if ((int)readval == pf.sub.bytes_per_subint) {
		    fprintf(stderr,"got %zu of %d\n",readval,pf.sub.bytes_per_subint);
		    psrfits_write_subint(&pf);
	    }
	    else {
		fprintf(stderr,"got %zu of %d\n",readval,pf.sub.bytes_per_subint);
		perror("Failed to read full allocation from pipe");
		break;
            }   
	    pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;
	    pf.sub.lst += pf.sub.tsubint;;

    } while (pf.T < pf.hdr.scanlen && !pf.status);

    //close the pipe
    close(fp);

    // Close the last file and cleanup
    fits_close_file(pf.fptr, &(pf.status));
    free(pf.sub.dat_freqs);
    free(pf.sub.dat_weights);
    free(pf.sub.dat_offsets);
    free(pf.sub.dat_scales);
    free(pf.sub.data);

    printf("Done.  Wrote %d subints (%f sec) in %d files.  status = %d\n", 
           pf.tot_rows, pf.T, pf.filenum, pf.status);

    exit(0);
}
