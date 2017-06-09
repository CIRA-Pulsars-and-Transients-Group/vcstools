

CFLAGS=-Wall


saxpy: saxpy.cu
	nvcc -o saxpy $^

wavenum: calcWaveNumbers.c
	gcc $(CFLAGS) -o $@ $^ -L/usr/lib -lm -lsla 
