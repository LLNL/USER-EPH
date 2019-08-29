NVCC=/usr/local/cuda/bin/nvcc
CXX=mpicxx 
CXXFLAGS= -O2 -std=c++11 -DFIX_EPH_GPU
NVCCFLAGS= -ccbin=g++ -O2 --std=c++11 
#                                            --compiler-bindir /usr/bin/gcc
INCFLAGS= -I/usr/local/cuda/include
LDFLAGS= -L/usr/local/cuda/lib64

all: eph_gpu.cpp eph_gpu.cu eph_spline_gpu.h eph_beta_gpu.h eph_gpu.h main.cpp
	$(NVCC) -c eph_gpu.cu  -o eph_gpu_cu.o $(INCFLAGS) $(NVCCFLAGS)
	$(NVCC) -c eph_gpu.cpp -o eph_gpu_cpp.o $(INCFLAGS) $(NVCCFLAGS)
	ar rcs libeph_gpu.a eph_gpu_cu.o eph_gpu_cpp.o
	$(CXX) $(CXXFLAGS) main.cpp libeph_gpu.a -o main -lcuda -lcudart $(INCFLAGS) $(LDFLAGS)

clean: 
	rm -f libeph_gpu.a eph_gpu_cu.o eph_gpu_cpp.o main

