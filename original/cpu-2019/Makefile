SHELL = /bin/sh

CC = gcc 

src = $(wildcard *.c)
obj = $(src:.c=.o)

LFLAGS = -fopenmp -L/usr/local/gfftw/lib -lfftw3_omp -lfftw3 -lm 

CFLAGS = -O3 

spin1D: $(obj)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS) 

$(obj) : binary.h

.PHONY: clean
clean:
	rm -f $(obj) spin1D









# -L tells the path to find libfftw3.a    
#  this library is double-precision FFTW library
#  -I tells the path for fftw3.h
#LIBS : -L/usr/local/lib -lfftw3 -I/usr/local/include

