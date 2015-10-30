#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <unistd.h>
#include <strings.h>
#include <math.h>

static int compare_floats(const void *a,const void *b) {
	
	float compare = *( const float*) a - *(const float*)b;
	if (compare > 0) {
		return 1;
	}
	else if (compare == 0) {
		return 0;
	}
	else {
		return -1;
	}
}
void usage() {
	fprintf(stderr,"get_weights -n <nchan> [128] -t <nsteps> [-1], takes input from stdin and dumps to stdout\n");
	fprintf(stderr,"options: -r <reference ant> [0]\n");
	fprintf(stderr,"\t -a nstation [128]\n");
	fprintf(stderr,"\t -n nchan [128]\n");
	fprintf(stderr,"\t -p npol [128]\n");
	fprintf(stderr,"\t -T tolerance [0.1] -- fraction of difference from median to reject\n");
	fprintf(stderr,"\t -t nsteps [-1]: number of timesteps (10000 per second) -1 is ALL\n");
        fprintf(stderr,"Expects ascii <real> <imaginary> EXCEPT if the -s flag is given\n");
	fprintf(stderr,"-s\t sigproc mode: read a binary single precision float total intensity filterbank. Caculates a shift and scale for the conversion to 8 bit unsigned chars that sigproc uses");
	exit(-1);
}

int main(int argc, char **argv) {

	int nchan =128;
	int nstation=32;
	int c = 0;
	int ch=0;
	int nsteps=-1;
	int npol = 2;
	int i = 0;
	int reference = 0;
        int sigproc = 0;
        float tolerance = 0.1;
	if (argc > 1) {

		while ((c = getopt(argc, argv, "a:hST:t:n:p:r:")) != -1) {
			switch(c) {
				case 'a':
					nstation = atoi(optarg);
					break;
				case 'h':
					usage();
					exit(-1);
					break;
				case 'n':
					  nchan = atoi(optarg);
					  break;
				case 'S':
					  sigproc = 1;
					  break;
				case 'T':
					  tolerance = atof(optarg);
					  break;
				case 't':
					  nsteps = atoi(optarg);
					  break;
				case 'p':
					  npol = atoi(optarg);
					  break;
				case 'r':
					 reference = atoi(optarg);
					 break;
				default:
					  usage();
					  break;
			}
	    }
	}



	float **spectra= (float **) calloc(nstation*npol,sizeof(float *));
	float *total_spec_power = (float *) calloc (nstation*npol,sizeof(float));
	float *weight = (float *) calloc (nstation*npol,sizeof(float));
	float *tosort = (float *) calloc (nstation*npol,sizeof(float));

	for (i=0;i<nstation*npol;i++) {
		spectra[i] = calloc(nchan,sizeof(float));
		bzero(spectra[i],nchan*sizeof(float));
	}
	if (sigproc) {
		int samp=0;
		float var=0;
		float mean = 0;
		float min = 0;
		float max = 0;
		int refchan = reference;
		float scale=0;
                float shift = 0;
		for (samp=0;samp<nsteps;samp++) {
			fread(spectra[0],nchan,sizeof(float),stdin);
			mean += spectra[0][refchan]/nsteps;
			if (spectra[0][refchan] < min) {
				min = spectra[0][refchan];
			}
			if (spectra[0][refchan] > max) {
				max = spectra[0][refchan];
			}

		}
		for (samp=0;samp<nsteps;samp++) {
			fread(spectra[0],nchan,sizeof(float),stdin);
			var += (spectra[0][refchan]-mean)*(spectra[0][refchan]-mean)/nsteps;
		}
				
		scale = 6.0*sqrt(var)/255.0;
				
		shift = -1.0*mean + (1.5*sqrt(var));


		fprintf(stderr,"MAX: %f MIN: %f MEAN: %f SIGMA: %f SCALE: %f SHIFT: %f\n",max,min,mean,sqrt(var),scale,shift);
		exit(EXIT_SUCCESS);

	}
	float real=0,cmp=0;
	int step = 0;
	// build up all the spectra
	while(!feof(stdin) && (nsteps < 0 || step < nsteps)) {

		for (i=0;i<nstation*npol;i++) {
			fprintf(stderr,"                                                       \r");
			fprintf(stderr,"Input step %d (input %d of %d)\r",step,i,nstation*npol);
			float complex spec = 0 + 0*I;
			for (ch=0;ch<nchan;ch++) {
				if ((fscanf(stdin,"%f %f\n",&real,&cmp)) == 2) {
					spec = real + I*cmp;
					spectra[i][ch] += cabsf(spec);
				}
			}
		}
		step++;
	}
	// finished building a total spectra for nsteps
	bzero(total_spec_power,(nstation*npol*sizeof(float)));
	for (i=0;i<nstation*npol;i++) {
		for (ch = 0; ch < nchan; ch++) {
			total_spec_power[i] += spectra[i][ch]/step;
		}
	}
	// above loop just averaged
	bzero(weight,(nstation*npol*sizeof(float)));
	bzero(tosort,(nstation*npol*sizeof(float)));
	for (i=0;i<nstation*npol;i++) {
		tosort[i] = 1.0/total_spec_power[i];
		weight[i] = tosort[i];
		
	}

	qsort(tosort,nstation*npol,sizeof(float),compare_floats);
	// get median
	
	float normalise = tosort[nstation*npol/2];

	for (i=0;i<nstation*npol;i++) {
		weight[i] = weight[i]/normalise;
		if ((weight[i] < (1-tolerance)) || (weight[i] > (1+tolerance))) { 
			weight[i] = 0.0;
		}
		fprintf(stdout,"%f \n",weight[i]);
	}


}
	
