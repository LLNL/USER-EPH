
NVCC=/usr/local/cuda/bin/nvc
NVCCFLAGS= -O2 --std=c++11 --compiler-bindir /usr/bin/gcc
INCFLAGS= -I/usr/local/cuda/include

all: eph.cpp eph.cu eph.h
	$(NVCC) -c eph.cu eph.cpp -o eph.o
	ar rcs libeph.a eph.o

