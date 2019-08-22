
NVCC=/usr/local/cuda/bin/nvcc
CXX=g++
CXXFLAGS= -O2 -std=c++11 
NVCCFLAGS= -O2 --std=c++11 --compiler-bindir /usr/bin/gcc
INCFLAGS= -I/usr/local/cuda/include
LDFLAGS= -L/usr/local/cuda/lib64

all: eph.cpp eph.cu eph.h
	$(NVCC) -c eph.cu  -o eph_cu.o $(INCFLAGS) $(NVCCFLAGS)
	$(NVCC) -c eph.cpp -o eph_cpp.o $(INCFLAGS) $(NVCCFLAGS)
	ar rcs libeph.a eph_cu.o eph_cpp.o
	$(CXX) main.cpp libeph.a -o main -lcuda -lcudart $(LDFLAGS)

clean: 
	rm libeph.a eph_cu.o eph_cpp.o main

