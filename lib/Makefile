#NVCC=/usr/local/cuda/bin/nvcc
ifndef CUDA_HOME
	CUDA_HOME=/usr/local/cuda
endif

CUDA_ARCH=compute_60
CUDA_CODE=sm_60

NVCC=$(CUDA_HOME)/bin/nvcc
MPICXX=mpicxx
CXX=$(MPICXX)

CXXFLAGS= -g -O3 -std=c++11 -DFIX_EPH_GPU
NVCCFLAGS= -g -O3 -ccbin=$(CXX) --std=c++11 -arch=$(CUDA_ARCH) -code $(CUDA_CODE) 

INCFLAGS= -I$(CUDA_HOME)/include -I../
LDFLAGS= -L$(CUDA_HOME)/lib64

all: eph_gpu.cpp eph_gpu.cu eph_spline_gpu.h eph_beta_gpu.h eph_gpu.h main.cpp
	$(NVCC) -c eph_gpu.cu  -o eph_gpu_cu.o $(INCFLAGS) $(NVCCFLAGS)
	$(NVCC) -c eph_gpu.cpp -o eph_gpu_cpp.o $(INCFLAGS) $(NVCCFLAGS)
	ar rcs libeph_gpu.a eph_gpu_cu.o eph_gpu_cpp.o
	$(CXX) $(CXXFLAGS) main.cpp libeph_gpu.a -o main -lcuda -lcudart $(INCFLAGS) $(LDFLAGS)

profile: eph_gpu.cpp eph_gpu.cu eph_spline_gpu.h eph_beta_gpu.h eph_gpu.h
	$(NVCC) -pg -c eph_gpu.cu  -o eph_gpu_cu.o $(INCFLAGS) $(NVCCFLAGS)
	$(NVCC) -pg -c eph_gpu.cpp -o eph_gpu_cpp.o $(INCFLAGS) $(NVCCFLAGS)
	ar rcs libeph_gpu.a eph_gpu_cu.o eph_gpu_cpp.o

clean: 
	rm -f libeph_gpu.a eph_gpu_cu.o eph_gpu_cpp.o main

