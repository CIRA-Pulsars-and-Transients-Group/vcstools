#include <stdio.h>
#include <fitsio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <glob.h>

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
/*--------------------------------------------------------------------------*/


void usage() {

	fprintf(stdout,"read_mwac -s -e -+ <list of files> \n");
	fprintf(stdout,"Simple utility to read the dump files from the correlator (or NGAS) " \
			" and dump the contents to stdout, you can choose to average and only read a subset of channels or both:\n");
    fprintf(stdout,"\t-c\t\t -- conjugate <all> visibilities\n");
	fprintf(stdout,"\t-s\t\t -- start channel\n");
	fprintf(stdout,"\t-e\t\t -- end channel\n");
	fprintf(stdout,"\t-+\t\t -- fscrunch factor\n");
	fprintf(stdout,"\t-a <root> \t\t -- glob the directory for files and reads them in order of creation (oldest to most recent). It will continue for until no new files appear\n");
	fprintf(stdout,"NEW RAW MODE\n");
	fprintf(stdout,"\t-r \"<nbaselines> <nchan>\" \t -- dimensions of raw cube expected on stdin\n");
}

int  checkFile(time_t index, char * fileName, time_t gap) {
	struct stat *s;
	//time_t accessTime;
	time_t modTime;
	//time_t changeTime;
	s = malloc(sizeof(struct stat));
	if(stat(fileName,s) == 0) {
		//accessTime = s->st_atime;
		modTime = s->st_mtime;      
		//changeTime = s->st_ctime;
		if (index-modTime < gap) {
			return 1;
		}
		else {
			return 0;
		}
	}
	free(s);
        return 0;
}
 
int main(int argc, char **argv) {

	int c;
	int rval = 0;

	/* --- parse arguments -----:
	 *
	 * input and output channels if the input channel number is smaller than the output
	 * then we assume that the input pipe is frequency ordered in the order that you want to combine them
	 *
	 * If the output is smaller than the input the output will be freuqnecy ordered.
	 *
	 * Therefore you should be able to put this back to back and split a big file into small ones then pipe 
	 * the small ones into a big one
	 *
	 * ---------------------------:
	 */  
	int startchan=0;
	int endchan=0;
	int fscrunch=1;
	int status=0;
	int nbaselines=0;
	int nchan=0;
	int raw=0;
	char *pattern = NULL;
    int conjugate = 0;
	int always_on = 0;

	while ((c = getopt(argc, argv, "Cr:s:e:+:a:")) != -1) {
		switch(c) {
			case 'a':
				pattern = strdup(optarg);
				always_on = 1;

				break;
            case 'c':
                conjugate = 1;
                break;
			case 'r':
				if ((sscanf(optarg,"%d %d",&nbaselines,&nchan)) != 2) {
					fprintf(stderr,"Using expert RAW mode must specify nbaselines and nchan in the output cube\n");
					usage();
				}
				raw = 1;
				break;
			case 'h':
				usage();
				exit(EXIT_SUCCESS);
			case 's':
				startchan=atoi(optarg);
				break;
			case 'e':
				endchan=atoi(optarg);
				break;
			case '+':
				fscrunch=atoi(optarg);
				break;
			default:
				usage();
				break;
		}
	}	
	if (raw==1) {

		fprintf(stderr,"Raw mode %d visibilites and %d channels\n",nbaselines,nchan);

		/* just stream data through usualy used for averaging */
		size_t nbuffer=0;
		size_t fpixel=0;
		int ii=0,f=0;
		float *buffer = NULL;
		int index = 0;

		nbuffer = nbaselines*nchan*2; // num floats per channel;

		fpixel = startchan * nbaselines;

		buffer = calloc (nbuffer,sizeof(float));

		while (fread(buffer,nbuffer*sizeof(float),1,stdin) != 0) {

			if (fscrunch > 1) {

				float *average = calloc(nbaselines*2,sizeof(float));

				float *buff = buffer+fpixel;
				ii=fpixel;
				do {

					for (f=0;f<fscrunch;f++) {

						for (index=0;index<nbaselines*2;index++) {

							average[index] = average[index] + (*buff/fscrunch);
							buff++;
							ii++;
						}

					}

					fwrite(average,sizeof(float),nbaselines*2,stdout);
					bzero((void *)average,nbaselines*2*sizeof(float));
				} while (ii<(int)nbuffer);
				free(average);

			}
			else {
				// write them all
				fwrite(buffer,sizeof(float),nbuffer,stdout);
			}
		}

		free(buffer);
	}

	else { 
		if (!always_on) {
			while (optind < argc) {

				status=0;
				fitsfile *fptr = NULL;
				char *filename = argv[optind];
				rval = fits_open_file(&fptr,filename,READONLY,&status);
				if (rval) {
					printerror(status);
				} 
				int anynull=0, hdupos=0;

				fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */
				status = 0;

				while(!status)  /* Main loop through each extension */
				{
					long naxes[2], fpixel=0, nbuffer=0, ii=0;
					int index=0,f=0;	
					float nullval  = 0;                /* don't check for null values in the image */
					int numchan=0;
					float *buffer = NULL;
					int hdutype = 0;

					fits_get_hdu_type(fptr, &hdutype, &status);

					/* read the NAXIS1 and NAXIS2 keyword to get image size */
					//if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
					naxes[0] = naxes[1] = 0;
					if (fits_get_img_size(fptr,2,naxes,&status))                 
						printerror( status );

					fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */

					if (naxes[0] != 0) {

						if (startchan != 0 || endchan != 0) {
							numchan = endchan - startchan + 1;
						}
						else {
							numchan = naxes[1];
						}

						nbuffer = naxes[0]*numchan; // num floats per channel;

						fpixel = startchan * naxes[0]  + 1;

						buffer = calloc (nbuffer,sizeof(float));


						/* Note that even though the FITS images contains unsigned integer */
						/* pixel values (or more accurately, signed integer pixels with    */
						/* a bias of 32768),  this routine is reading the values into a    */
						/* float array.   Cfitsio automatically performs the datatype      */
						/* conversion in cases like this.                                  */

						if ( fits_read_img(fptr, TFLOAT, fpixel, nbuffer, &nullval,
									buffer, &anynull, &status) )
							printerror( status );

						if (fscrunch > 1) {

							float *average = calloc(naxes[0],sizeof(float));

							float *buff = buffer;
							ii=0;
							do {

								for (f=0;f<fscrunch;f++) {

									for (index=0;index<naxes[0];index++) {

										average[index] = average[index] + (*buff/fscrunch);
										buff++;
										ii++;
									}

								}

								fwrite(average,sizeof(float),naxes[0],stdout);
								bzero((void *)average,naxes[0]*sizeof(float));
							} while (ii<nbuffer);
							free(average);

						}
                        else if (conjugate){
                            for (index=0;index<nbuffer;index=index+2) {
                                buffer[index] = -1*buffer[index];
                            }
                            fwrite(buffer,sizeof(float),nbuffer,stdout);

                        }
						else {
							// write them all
							fwrite(buffer,sizeof(float),nbuffer,stdout);
						}

						free(buffer);
					}
					fits_movrel_hdu(fptr, 1, NULL, &status);  /* try to move to next HDU */
				}
				status=0;
				optind++;
				fits_close_file(fptr,&status);	

			}
		} // always on
		else if (always_on ==1) {
			time_t now = time(NULL);
			glob_t globbuf;
			glob(pattern, 0, NULL, &globbuf);
			if(globbuf.gl_pathc == 0)
				printf("there were no matching files");
			else {
				int i;
				for(i = 0; i < (int)globbuf.gl_pathc; i++) {
					if (checkFile(now,globbuf.gl_pathv[i],10)) {
						printf("%s \n",globbuf.gl_pathv[i]);
					}
				}
			}
		}
	}
	return 0;
}
