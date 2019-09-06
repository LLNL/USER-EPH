
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>

#define FIX_EPH_GPU

#include "eph_spline_gpu.h"
#include "eph_beta_gpu.h"
#include "eph_gpu.h"

__device__
double get_scalar(double3d x, double3d y) 
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

__device__
double get_norm(double3d x) 
{
  return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}

__device__
double get_distance_sq(double3d x, double3d y) 
{
  double3d dxy;
  dxy[0] = x[0] - y[0];
  dxy[1] = x[1] - y[1];
  dxy[2] = x[2] - y[2]; 
  
  return get_norm(dxy);
}

__device__
double get_difference_sq(const double3d& x, const double3d& y, double3d& z) 
{
  z[0] = x[0] - y[0];
  z[1] = x[1] - y[1];
  z[2] = x[2] - y[2];
  
  return get_norm(z);
}

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
  int thread_index = threadIdx.x;
  int block_index = blockIdx.x;
  int block_dimension = blockDim.x;
  int grid_dimension = gridDim.x;
  
  EPH_Beta_GPU& beta = *((EPH_Beta_GPU*) eph_gpu.beta_gpu);
  
  int nlocal = eph_gpu.nlocal;
  int nghost = eph_gpu.nghost;
  
  double3d* x = eph_gpu.x_gpu;
  
  int* type = eph_gpu.type_gpu;
  int* mask = eph_gpu.mask_gpu;
  int groupbit = eph_gpu.groupbit;
  
  int* type_map = eph_gpu.type_map_gpu;
  
  double* rho_i_gpu = eph_gpu.rho_i_gpu;
  
  double r_cutoff_sq = eph_gpu.r_cutoff_sq;
  
  int index = block_index * block_dimension + thread_index;
  int stride = block_dimension * grid_dimension;
  
  int* number_neigh = eph_gpu.number_neigh_gpu;
  int* index_neigh = eph_gpu.index_neigh_gpu;
  int* neighs = eph_gpu.neighs_gpu;
  
  for(size_t i = index; i < nlocal; i += stride)
  {
    if(!(mask[i] & groupbit)) continue;
    
    rho_i_gpu[i] = 0;
    
    // neighbour list version
    for(int j = 0; j < number_neigh[i]; ++j)
    {
      int jj = neighs[index_neigh[i] + j];
      
      double r_sq = get_distance_sq(x[i], x[jj]);
      
      if(r_sq < r_cutoff_sq)
      {
        rho_i_gpu[i] += beta.get_rho_r_sq(type_map[type[jj] - 1], r_sq);
      }
    }
    
    // brute force scan
    /*
    for(size_t j = 0; j < (nlocal + nghost); ++j)
    {
      if(i == j) continue;
      
      double r_sq = get_distance_sq(x[i], x[j]);
      
      if(r_sq < r_cutoff_sq)
      {
        rho_i_gpu[i] += beta.get_rho_r_sq(type_map[type[j] - 1], r_sq);
      }
    }
    */
  }
}

void calculate_environment_gpu(EPH_GPU& eph_gpu)
{
  calculate_environment_cu<<<256, 256>>>(eph_gpu);
  cudaDeviceSynchronize();
}



