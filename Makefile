pabeam: pabeam.c
	gcc -Wall -O3 -o $@ $^ -L/usr/lib -lm -lsla -lcfitsio


pabeam_gpu: pabeam.cu
	nvcc -g -G -O3 -o $@ $^ -L/usr/lib -lm -lsla -lcfitsio 


saxpy: saxpy.cu
	nvcc -o saxpy $^


