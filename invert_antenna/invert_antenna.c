#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <complex.h>
#include <fftw3.h>

void usage () {
	fprintf(stderr,"invert_antenna -d [+1|-1] apply an FFT to stdin\n");
}
int main(int argc, char **argv) {

	int nchan =128;
	int direction = FFTW_FORWARD;
	int c;
	fftwf_complex *in,*out;
	fftwf_plan p;
	int real_in = 0;
	if (argc > 1) {

		while ((c = getopt(argc, argv, "hrn:d:")) != -1) {
			switch(c) {
				case 'h':
					usage();
					exit(-1);
					break;
				case 'd': {
						  int dir = atoi(optarg);
						  if (dir == -1) {
							  direction = FFTW_BACKWARD;
						  }
						  else if (dir == 1) {
							  direction = FFTW_FORWARD;
						  }
					  }
					  break;
				case 'r':
					real_in = 1;
					break;
				case 'n':
					  nchan = atoi(optarg);
					  break;
				default:
					  usage();
					  break;
			}
	    }
    }

	
		
	in = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * nchan);
	out = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * nchan);
	if (real_in)
		p = fftwf_plan_dft_r2c_1d(nchan,(float *) in,out,FFTW_ESTIMATE);
	else 
		p = fftwf_plan_dft_1d(nchan,in,out,direction,FFTW_ESTIMATE);

	int ch=0;
	float real,cmp;
	while(!feof(stdin)) {

		for (ch=0;ch<nchan;ch++) {
			fscanf(stdin,"%f %f\n",&real,&cmp);
			in[ch] = real + I*cmp;
		}
		if (feof(stdin))
			break;
		fftwf_execute(p);

		for (ch = 0; ch < nchan; ch++) {
			fprintf(stdout,"%f %f\n",crealf(out[ch]),cimagf(out[ch]));
		}

	}
}
	
