

CFLAGS=-Wall -O3


wavenum: calcWaveNumbers.c
	gcc $(CFLAGS) -o $@ $^ -L/usr/lib -lm -lsla 


saxpy: saxpy.cu
	nvcc -o saxpy $^


