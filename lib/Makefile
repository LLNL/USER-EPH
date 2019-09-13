#NVCC=/usr/local/cuda/bin/nvcc
ifndef CUDA_HOME
	CUDA_HOME=/usr/local/cuda
endif

CUDA_ARCH=compute_70
CUDA_CODE=sm_70

NVCC=$(CUDA_HOME)/bin/nvcc
CXX=mpicxx 
CXXFLAGS= -g -O2 -std=c++11 -DFIX_EPH_GPU
NVCCFLAGS= -ccbin=g++ -g -O2 --std=c++11 -arch=$(CUDA_ARCH) -code $(CUDA_CODE) 

INCFLAGS= -I$(CUDA_HOME)/include -I../
LDFLAGS= -L$(CUDA_HOME)/lib64

all: eph_gpu.cpp eph_gpu.cu eph_spline_gpu.h eph_beta_gpu.h eph_gpu.h main.cpp
	$(NVCC) -c eph_gpu.cu  -o eph_gpu_cu.o $(INCFLAGS) $(NVCCFLAGS)
	$(NVCC) -c eph_gpu.cpp -o eph_gpu_cpp.o $(INCFLAGS) $(NVCCFLAGS)
	ar rcs libeph_gpu.a eph_gpu_cu.o eph_gpu_cpp.o
	$(CXX) $(CXXFLAGS) main.cpp libeph_gpu.a -o main -lcuda -lcudart $(INCFLAGS) $(LDFLAGS)

clean: 
	rm -f libeph_gpu.a eph_gpu_cu.o eph_gpu_cpp.o main

