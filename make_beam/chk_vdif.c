#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<assert.h>
#include<string.h>

#include<vdifio.h>

void usage() {
	fprintf(stderr,"chk_vdif - a simple utility to check the time stamps of a vdif file and test for discontinuites\n");
	exit(-1);
}

int main(int argc, char **argv) {

    int c=0;
    int check = 0;
    int fix = 0;
    int verbose = 0;
    int test = 0;
    char *infile = NULL;
    FILE *input = stdin;
    FILE *output = stdout;
    int rate = 10000;
    
	if (argc > 1) {

		while ((c = getopt(argc, argv, "i:ctfr:v")) != -1) {
			switch(c) {
				case 'c':
					check = 1;
					break;
                    
				case 'f':
                    fix = 1;
                    break;
                case 'i':
                    infile = strdup(optarg);
                    input = fopen(infile,"r");
                    break;
                case 'r':
                    rate = atoi(optarg);
                    break;
                case 't':
                    test = 1;
                    break;
                case 'v':
                    verbose = 1;
                    break;
                default:
                    usage();
                    break;
                    
            }
        }
	}
	else {
		usage();
	}
	
    vdif_header vhdr;
    vdif_header last_vhdr;
    
	unsigned int frame = 0;
	int finished = 0;
	unsigned int rtn = 0;
    int got_first = 0;
    char *buffer = NULL;
    char *pattern = NULL;
    int second = 0;
    int last_second = 0;
    int last_frame=0;
    double MJD = 0.0;
    unsigned long frame_count = 0;
    int failed = 0;
    int diff_frame=0;
    int diff_seconds=0;
    
	if (check || fix || test) {
        
		while (!finished) {
            
			rtn = fread(&vhdr,32,1,input);
            
			if (feof(input) || rtn != 1) {
                if (verbose)
				fprintf(stderr,"chk_vdif finished:\n");
				finished = 1;
				continue;
			}
            else {
                frame_count = frame_count + 1;
            }
            if (got_first == 0) {
                if (check && verbose)
                    fprintf(stdout,"## Frame length: %d",vhdr.framelength8*8);
    
                buffer = malloc(vhdr.framelength8*8);
                assert(buffer);
                if (fix) {
                    // initialise to a mean zero 2bit offset binary pattern
                    pattern = malloc(vhdr.framelength8*8);
                    assert(pattern);
                    uint16_t mean_zero = 0xAA55;
                    size_t count = 0;
                    while (count < vhdr.framelength8*8) {
                        memcpy(pattern+count,&mean_zero,2);
                        count = count+2;
                    }
                    
                }
                got_first = 1;
            }
            else {
                
                frame = getVDIFFrameNumber(&vhdr);
                if (check && verbose)
                    fprintf(stdout,"Frame: %d\n",frame);
                second = getVDIFFrameSecond(&vhdr);
                if (check && verbose)
                    fprintf(stdout,"Second: %d\n", second) ;
                MJD = getVDIFDMJD(&vhdr,10000);
                if (check && verbose)
                    fprintf(stdout,"MJD %f\n",MJD);
                
                if (frame_count > 2) {
                    diff_frame = frame-last_frame;
                    if (diff_frame > 1) {
                        if (check)
                            fprintf(stderr,"Skipped a frame between %d and %d\n",last_frame,frame);
                        
                        failed = 1;
                    }
                    
                    diff_seconds = second-last_second;
                    if (diff_seconds > 1) {
                            if (check)
                                fprintf(stderr,"Skipped a second between %d and %d\n",last_second,second);
                        if (fix) {
                            while (diff_seconds > 1) {
                                nextVDIFHeader(&last_vhdr,rate);
                                fwrite(&last_vhdr,VDIF_HEADER_BYTES,1,stdout);
                                // repeat the last buffer
                                fwrite(&buffer,(vhdr.framelength8*8-VDIF_HEADER_BYTES),1,output);
                                last_second = getVDIFFrameSecond(&last_vhdr);
                                diff_seconds = second-last_second;
                            }
                        }
                        if (fix == 0)
                            failed = 1;
                    }
                   
                    
                }
                
            }
            
            
            if (failed == 1) {
                fprintf(stdout,"failed\n");
                finished = 1;
                continue;
            }
            rtn = fread(buffer,(vhdr.framelength8*8-VDIF_HEADER_BYTES),1,input);
            if (feof(input) || rtn != 1) {
                if (verbose)
                fprintf(stderr,"chk_vdif finished:\n");
                finished = 1;
                continue;
            }
            
            last_frame = frame;
            last_second = second;
            last_vhdr = vhdr;
            
            if (fix) {
                fwrite(&vhdr,VDIF_HEADER_BYTES,1,output);
                fwrite(buffer,(vhdr.framelength8*8-VDIF_HEADER_BYTES),1,output);
            }
			
		}
	}
    free(buffer);
    if (test && failed == 0) {
        fprintf(stdout,"passed\n");
    }
    if (infile) {
        fclose(input);
    }
}
