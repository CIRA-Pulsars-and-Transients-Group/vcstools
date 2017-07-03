pabeam: pabeam.c
	gcc -Wall -O3 -o $@ $^ -L/usr/lib -lm -lsla -lcfitsio


pabeam_gpu: pabeam.cu
	nvcc -O3 -o $@ $^ -arch=sm_35 -L/usr/lib -lm -lsla -lcfitsio

