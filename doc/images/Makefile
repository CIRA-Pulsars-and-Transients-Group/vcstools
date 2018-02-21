CFLAGS = -Wall -Wextra
LDLIBS = -lm

TARGETS = weights_abs.png \
          weights_arg.png \
		  filter_coeffs.png

all: $(TARGETS)

weights_abs.png weights_arg.png: weights.gpi weights_abs.txt weights_arg.txt
	gnuplot $<

weights_abs.txt weights_arg.txt: weights filter_coeffs.txt
	./$<

weights: weights.o filter.o

filter_coeffs.png: filter_coeffs.gpi
	gnuplot $<

filter.h filter.c:
	ln -s ../../make_beam_small/$@ .
