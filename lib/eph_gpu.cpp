
#include <cuda.h>
#include <cuda_runtime_api.h>

#include <iostream>

#define FIX_EPH_GPU
#include "eph_gpu.h"

void calculate_environment_gpu();

EPH_GPU allocate_EPH_GPU(size_t n)
{
  EPH_GPU eph_gpu;
  eph_gpu.nlocal = n;
  return eph_gpu;
}

void deallocate_EPH_GPU(EPH_GPU& eph_gpu)
{
  if(eph_gpu.x_gpu != nullptr) cudaFree(eph_gpu.x_gpu);
}

