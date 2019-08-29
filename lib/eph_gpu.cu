
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>

#define FIX_EPH_GPU

#include "eph_spline_gpu.h"
#include "eph_beta_gpu.h"
#include "eph_gpu.h"

__global__ 
void test_interpolator_cu(EPH_Spline_GPU spl)
{
  //int thread_index = threadIdx.x;
  //int block_dimension = blockDim.x;
  //int grid_dimension = gridDim.x;
  
  printf("GPU: %.2f %.2f %.2f %.2f\n", 
    spl(0.0), spl(1.0), spl(1.5), spl(2.0));
}

void test_interpolator_gpu(EPH_Spline_GPU &spl) {
  test_interpolator_cu<<<1, 1>>>(spl);
  cudaDeviceSynchronize();
}

__global__
void test_beta_rho_cu(EPH_Beta_GPU beta)
{
  printf("GPU: %.3f %.3f %.3f %.3f\n", 
    beta.get_rho(0, 1.0), beta.get_rho_r_sq(0, 1.0), 
    beta.get_alpha(0, 0.1), beta.get_beta(0, 0.1));
}

void test_beta_rho_gpu(EPH_Beta_GPU &beta)
{
  test_beta_rho_cu<<<1, 1>>>(beta);
  cudaDeviceSynchronize();
}

__global__
void calculate_environment_cu(EPH_GPU eph_gpu)
{
  printf("calculating environment\n");
}

void calculate_environment_gpu(EPH_GPU& eph_gpu)
{
  calculate_environment_cu<<<1, 1>>>(eph_gpu);
}



