all: pabeam_gpu

pabeam_gpu: pabeam.cu
	nvcc -O3 -o $@ $^ -arch=sm_35 -L/usr/lib -lm -lsla -lcfitsio

test: pabeam_gpu
	./pabeam_gpu -f 184.96e6 -r "05:34:31.97" -d "+22:00:52.06" -t "2014-11-07T16:53:20" -m 1099414416_metafits_ppds.fits -b flagged_tiles.txt -x 0.01 -y 0.05

profile: pabeam_gpu
	nvprof ./pabeam_gpu -f 184.96e6 -r "05:34:31.97" -d "+22:00:52.06" -t "2014-11-07T16:53:20" -m 1099414416_metafits_ppds.fits -b flagged_tiles.txt -x 0.1 -y 0.1


pabeam: pabeam.c
	gcc -Wall -O3 -o $@ $^ -L/usr/lib -lm -lsla -lcfitsio



