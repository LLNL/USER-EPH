
#include <cuda.h>
#include <cuda_runtime_api.h>

#define FIX_EPH_GPU
#include "eph_gpu.h"

void calculate_environment_gpu();

// TODO: this file does not need much

EPH_Spline_GPU::EPH_Spline_GPU(double dx, std::vector<double> &y) :
  Spline(dx, y)
{
  // TODO: error checks
  cudaMalloc((void**) &c_gpu, c.size() * sizeof(Coefficients));
  cudaMemcpy(c_gpu, c.data(), c.size() * sizeof(Coefficients), cudaMemcpyHostToDevice);
  // cleanup unused memory
  c.resize(0);
}

EPH_Spline_GPU::EPH_Spline_GPU(Spline& spl) :
  Spline(spl)
{
  cudaMalloc((void**) &c_gpu, c.size() * sizeof(Coefficients));
  cudaMemcpy(c_gpu, c.data(), c.size() * sizeof(Coefficients), cudaMemcpyHostToDevice);
  // cleanup unused memory
  c.resize(0);
}

EPH_Spline_GPU::~EPH_Spline_GPU()
{
  cudaFree(c_gpu);
}


