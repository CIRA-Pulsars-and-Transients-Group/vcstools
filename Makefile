pabeam: pabeam.c
	gcc -Wall -O3 -o $@ $^ -L/usr/lib -lm -lsla -lcfitsio


pabeam_gpu: pabeam.cu
	nvcc -O3 -o $@ $^ --gpu-architecture=compute_35 -L/usr/lib -lm -lsla -lcfitsio -lcublas


saxpy: saxpy.cu
	nvcc -o saxpy $^


