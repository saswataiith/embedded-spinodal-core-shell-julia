SHELL = /bin/sh
NVCC = nvcc
CUDAFLAGS = -Xptxas -O3
ARCHFLAG = -arch=sm_70   
LIBS = -lcufft -lcurand -lcudart
INCDIRS = -I/usr/local/cuda-11.2/include
LIBDIRS = -L/usr/local/cuda-11.2/lib64
HEADER = binary.h

SRC  = src_3d
EXTR = ext
all: build

build: main.cu $(HEADER) 
	$(NVCC) $(CUDAFLAGS) $(ARCHFLAG) $(LIBDIRS) $(INCDIRS) $(LIBS) main.cu -o sp3d.out

#work_dir:
#	cp $(EXTR)/bin1ary ../
#	cp $(EXTR)/InputParams ../

clean:
	\rm -f profile.in
	\rm -f prof_gp.*
	\rm sp2d.out
