
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
void calculate_environment_cu(EPH_GPU eph_gpu)
{
  int thread_index = threadIdx.x;
  int block_index = blockIdx.x;
  int block_dimension = blockDim.x;
  int grid_dimension = gridDim.x;
  
  EPH_Beta_GPU& beta = *((EPH_Beta_GPU*) eph_gpu.beta_gpu);
  
  int nlocal = eph_gpu.nlocal;
  //int nghost = eph_gpu.nghost;
  
  double3d* x = eph_gpu.x_gpu;
  
  int* type = eph_gpu.type_gpu;
  int* mask = eph_gpu.mask_gpu;
  int groupbit = eph_gpu.groupbit;
  
  int* type_map = eph_gpu.type_map_gpu;
  
  double* rho_i = eph_gpu.rho_i_gpu;
  
  double r_cutoff_sq = eph_gpu.r_cutoff_sq;
  
  int index = block_index * block_dimension + thread_index;
  int stride = block_dimension * grid_dimension;
  
  int* number_neigh = eph_gpu.number_neigh_gpu;
  int* index_neigh = eph_gpu.index_neigh_gpu;
  int* neighs = eph_gpu.neighs_gpu;
  
  for(size_t i = index; i < nlocal; i += stride)
  {
    if(not(mask[i] & groupbit)) continue;
    
    rho_i[i] = 0;
    
    // neighbour list version
    for(int j = 0; j != number_neigh[i]; ++j)
    {
      int jj = neighs[index_neigh[i] + j];
      jj &= eph_gpu.neighmask;
      
      double r_sq = get_distance_sq(x[i], x[jj]);
      
      if(r_sq < r_cutoff_sq)
      {
        rho_i[i] += beta.get_rho_r_sq(type_map[type[jj] - 1], r_sq);
      }
    }
  }
}

void calculate_environment_gpu(EPH_GPU& eph_gpu)
{
  calculate_environment_cu<<<256, 256>>>(eph_gpu);
  cudaDeviceSynchronize();
}

__global__
void force_prl_stage1_cu(EPH_GPU eph_gpu)
{
  int thread_index = threadIdx.x;
  int block_index = blockIdx.x;
  int block_dimension = blockDim.x;
  int grid_dimension = gridDim.x;
  
  EPH_Beta_GPU& beta = *((EPH_Beta_GPU*) eph_gpu.beta_gpu);
  
  int nlocal = eph_gpu.nlocal;
  int nghost = eph_gpu.nghost;
  
  double3d* x = eph_gpu.x_gpu;
  double3d* v = eph_gpu.v_gpu;
  double3d* w_i = eph_gpu.w_i_gpu;
  
  int* type = eph_gpu.type_gpu;
  int* mask = eph_gpu.mask_gpu;
  int groupbit = eph_gpu.groupbit;
  
  int* type_map = eph_gpu.type_map_gpu;
  
  double* rho_i = eph_gpu.rho_i_gpu;
  
  double r_cutoff_sq = eph_gpu.r_cutoff_sq;
  
  int index = block_index * block_dimension + thread_index;
  int stride = block_dimension * grid_dimension;
  
  int* number_neigh = eph_gpu.number_neigh_gpu;
  int* index_neigh = eph_gpu.index_neigh_gpu;
  int* neighs = eph_gpu.neighs_gpu;
  
  for(size_t i = index; i < nlocal; i += stride)
  {
    if(!(mask[i] & groupbit)) continue;
    
    w_i[i][0] = 0;
    w_i[i][1] = 0;
    w_i[i][2] = 0;
    
    if(not(rho_i[i] > 0)) continue;
    
    double alpha_i = beta.get_alpha(type_map[type[i] - 1], rho_i[i]);
    
    for(int j = 0; j != number_neigh[i]; ++j)
    {
      int jj = neighs[index_neigh[i] + j];
      jj &= eph_gpu.neighmask;
      
      double3d e_ij;
      double e_r_sq = get_difference_sq(x[jj], x[i], e_ij);
      
      if(e_r_sq > r_cutoff_sq) continue;
      
      double v_rho_ji = beta.get_rho_r_sq(type_map[type[jj] - 1], e_r_sq);
      double prescaler = alpha_i * v_rho_ji / (rho_i[i] * e_r_sq);
      
      double e_v_v1 = get_scalar(e_ij, v[i]);
      double var1 = prescaler * e_v_v1;
      
      double e_v_v2 = get_scalar(e_ij, v[jj]);
      double var2 = prescaler * e_v_v2;
      
      double dvar = var1 - var2;
      w_i[i][0] += dvar * e_ij[0];
      w_i[i][1] += dvar * e_ij[1];
      w_i[i][2] += dvar * e_ij[2];
    }
  }
}

__global__
void force_prl_stage2_cu(EPH_GPU eph_gpu)
{
  int thread_index = threadIdx.x;
  int block_index = blockIdx.x;
  int block_dimension = blockDim.x;
  int grid_dimension = gridDim.x;
  
  EPH_Beta_GPU& beta = *((EPH_Beta_GPU*) eph_gpu.beta_gpu);
  
  int nlocal = eph_gpu.nlocal;
  int nghost = eph_gpu.nghost;
  
  double3d* x = eph_gpu.x_gpu;
  double3d* v = eph_gpu.v_gpu;
  
  double3d* f_EPH = eph_gpu.f_EPH_gpu;
  double3d* w_i = eph_gpu.w_i_gpu;
  
  int* type = eph_gpu.type_gpu;
  int* mask = eph_gpu.mask_gpu;
  int groupbit = eph_gpu.groupbit;
  
  int* type_map = eph_gpu.type_map_gpu;
  
  double* rho_i = eph_gpu.rho_i_gpu;
  
  double r_cutoff_sq = eph_gpu.r_cutoff_sq;
  
  int index = block_index * block_dimension + thread_index;
  int stride = block_dimension * grid_dimension;
  
  int* number_neigh = eph_gpu.number_neigh_gpu;
  int* index_neigh = eph_gpu.index_neigh_gpu;
  int* neighs = eph_gpu.neighs_gpu;
  
  for(size_t i = index; i < nlocal; i += stride)
  {
    if(!(mask[i] & groupbit)) continue;
    
    f_EPH[i][0] = 0;
    f_EPH[i][1] = 0;
    f_EPH[i][2] = 0;
    
    if(not(rho_i[i] > 0)) continue;
    
    double alpha_i = beta.get_alpha(type_map[type[i] - 1], rho_i[i]);
    
    for(int j = 0; j != number_neigh[i]; ++j)
    {
      int jj = neighs[index_neigh[i] + j];
      jj &= eph_gpu.neighmask;
      
      double3d e_ij;
      double e_r_sq = get_difference_sq(x[jj], x[i], e_ij);
      
      if(e_r_sq > r_cutoff_sq or not(rho_i[jj] > 0)) continue;
      
      double alpha_j = beta.get_alpha(type_map[type[jj] - 1], rho_i[jj]);
      
      double v_rho_ji = beta.get_rho_r_sq(type_map[type[jj] - 1], e_r_sq);
      double e_v_v1 = get_scalar(e_ij, w_i[i]);
      double var1 = alpha_i * v_rho_ji * e_v_v1 / (rho_i[i] * e_r_sq);
      
      double v_rho_ij = beta.get_rho_r_sq(type_map[type[i] - 1], e_r_sq);
      double e_v_v2 = get_scalar(e_ij, w_i[jj]);
      double var2 = alpha_j * v_rho_ij * e_v_v2 / (rho_i[jj] * e_r_sq);
      
      double dvar = var1 - var2;
      // friction is negative!
      f_EPH[i][0] -= dvar * e_ij[0];
      f_EPH[i][1] -= dvar * e_ij[1];
      f_EPH[i][2] -= dvar * e_ij[2];
    }
  }
}

__global__
void force_prl_stage3_cu(EPH_GPU eph_gpu)
{
  int thread_index = threadIdx.x;
  int block_index = blockIdx.x;
  int block_dimension = blockDim.x;
  int grid_dimension = gridDim.x;
  
  EPH_Beta_GPU& beta = *((EPH_Beta_GPU*) eph_gpu.beta_gpu);
  
  int nlocal = eph_gpu.nlocal;
  int nghost = eph_gpu.nghost;
  
  double3d* x = eph_gpu.x_gpu;
  double3d* xi_i = eph_gpu.xi_i_gpu;
  
  double3d* f_RNG = eph_gpu.f_RNG_gpu;
  
  int* type = eph_gpu.type_gpu;
  int* mask = eph_gpu.mask_gpu;
  int groupbit = eph_gpu.groupbit;
  
  int* type_map = eph_gpu.type_map_gpu;
  
  double* rho_i = eph_gpu.rho_i_gpu;
  
  double r_cutoff_sq = eph_gpu.r_cutoff_sq;
  
  int index = block_index * block_dimension + thread_index;
  int stride = block_dimension * grid_dimension;
  
  int* number_neigh = eph_gpu.number_neigh_gpu;
  int* index_neigh = eph_gpu.index_neigh_gpu;
  int* neighs = eph_gpu.neighs_gpu;
  
  for(size_t i = index; i < nlocal; i += stride)
  {
    if(!(mask[i] & groupbit)) continue;
    
    f_RNG[i][0] = 0;
    f_RNG[i][1] = 0;
    f_RNG[i][2] = 0;
    
    if(not(rho_i[i] > 0)) continue;
  }
}

// W * v
void force_prl_stage1_gpu(EPH_GPU& eph_gpu)
{
  force_prl_stage1_cu<<<256, 256>>>(eph_gpu);
  cudaDeviceSynchronize();
}

// W * wi
void force_prl_stage2_gpu(EPH_GPU& eph_gpu)
{
  force_prl_stage2_cu<<<256, 256>>>(eph_gpu);
  cudaDeviceSynchronize();
}

// W * xi
void force_prl_stage3_gpu(EPH_GPU& eph_gpu)
{
  force_prl_stage3_cu<<<256, 256>>>(eph_gpu);
  cudaDeviceSynchronize();
}

