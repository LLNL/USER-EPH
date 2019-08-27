
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>

#define FIX_EPH_GPU

#include "eph_gpu.h"

__global__
void dummy_test(int myID, int nrPS)
{
  int thread_index = threadIdx.x;
  int block_dimension = blockDim.x;
  int grid_dimension = gridDim.x;

  printf("GPU %d %d %d %d %d\n", myID, nrPS, thread_index, block_dimension, grid_dimension);
}

void run_dummy_test_cu(int myID, int nrPS) 
{
  dummy_test<<<1, 1>>>(myID, nrPS);
}

__global__ 
void test_interpolator_cu(EPH_Spline_GPU spl, double *values, int n_values)
{
  //int thread_index = threadIdx.x;
  //int block_dimension = blockDim.x;
  //int grid_dimension = gridDim.x;

  printf("GPU: %.2f %.2f %.2f %.2f\n", spl(0.0), spl(1.0), spl(1.5), spl(2.0));
}

void test_interpolator_gpu(EPH_Spline_GPU spl, double *values, int n_values) {
  test_interpolator_cu<<<1, 1>>>(spl, values, n_values);
}

__global__
void calculate_environment_gpu_cu()
{
  
}

void calculate_environment_gpu()
{
  
}



