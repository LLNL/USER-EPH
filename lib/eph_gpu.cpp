
#include <cuda.h>
#include <cuda_runtime_api.h>

#define FIX_EPH_GPU
#include "eph_gpu.h"

void run_dummy_test_cu(int, int);

void call_dummy(int myID, int nrPS) 
{
  double *f;
  int n = 1024;

  cudaMalloc((void**) &f, n*sizeof(double));
  run_dummy_test_cu(myID, nrPS);
  cudaFree(f);
}

void calculate_environment_gpu();

EPH_Spline_GPU::EPH_Spline_GPU(double dx, std::vector<double> &y) :
  Spline(dx, y)
{
  // TODO: error checks
  cudaMalloc((void**) &c_gpu, y.size() * sizeof(Coefficients));
  
  cudaMemcpy(c_gpu, c.data(), c.size() * sizeof(Coefficients), cudaMemcpyHostToDevice);
}

EPH_Spline_GPU::EPH_Spline_GPU(Spline& spl)
{
  
}

EPH_Spline_GPU::~EPH_Spline_GPU()
{
  cudaFree(c_gpu);
}


