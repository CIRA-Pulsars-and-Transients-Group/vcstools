

CFLAGS=-Wall -O3


pabeam: pabeam.c
	gcc $(CFLAGS) -o $@ $^ -L/usr/lib -lm -lsla -lcfitsio


saxpy: saxpy.cu
	nvcc -o saxpy $^


