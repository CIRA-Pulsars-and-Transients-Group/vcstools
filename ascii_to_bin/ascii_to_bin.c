#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

void usage() {

	fprintf(stderr,"ascii_to_bin -s <scale> -t <type> [1]\n");
	fprintf(stderr,"Read an ascii file from stdin and parse it as a list of \
types the types are:\n");
	fprintf(stderr,"\t1--32 bit float\n");
	fprintf(stderr,"\t2--32 bit int\n");
	fprintf(stderr,"\t3--16 bit int\n");
	exit(-1);
}
int main(int argc, char **argv) {
	int type=1;
	int scale=1;
	char buffer[64];

	int c;

	if (argc > 1) {

	    while ((c = getopt(argc, argv, "hs:t:")) != -1) {
		    switch(c) {
			case 'h':
				usage();
				break;
			case 's':
				scale = atoi(optarg);
				break;
			case 't':
				type = atoi(optarg);
				break;
			default:
				usage();
				break;
			}
	    }
	}

	
	
	while (!feof(stdin)) {
		fscanf(stdin,"%s",buffer);
		if (feof(stdin))
			break;
		if (type == 1) {
			float number = atof(buffer);
			number *= scale;
			fwrite(&number,sizeof(float),1,stdout);
		}
		else if (type == 2) {
			float temp = atof(buffer);
			temp *= scale;
			int32_t number = (int32_t) temp;;
			fwrite(&number,sizeof(int32_t),1,stdout);
		}
		else if (type == 3) {
			float temp = atof(buffer);

			temp *= scale;
			int16_t number = (int16_t) temp;
			fwrite(&number,sizeof(int16_t),1,stdout);
		}
	}
}
		
