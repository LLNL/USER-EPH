
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

/** ZERO DATA ON DEVICE **/
__global__
void zero_data_cu(EPH_GPU eph_gpu)
{
  int thread_index = threadIdx.x;
  int block_index = blockIdx.x;
  int block_dimension = blockDim.x;
  int grid_dimension = gridDim.x;
  
  double3d *f_EPH = eph_gpu.f_EPH_gpu;
  double3d *f_RNG = eph_gpu.f_RNG_gpu;
  
  double3d *xi_i = eph_gpu.xi_i_gpu;
  double3d *w_i = eph_gpu.w_i_gpu;
  
  double *rho_i = eph_gpu.rho_i_gpu;
  double *T_e_i = eph_gpu.T_e_i_gpu;
  
  int index = block_index * block_dimension + thread_index;
  int stride = block_dimension * grid_dimension;
  
  int nlocal = eph_gpu.nlocal;
  int nghost = eph_gpu.nghost;
  int ntotal = nlocal + nghost;
  
  for(size_t i = index; i < ntotal; i += stride)
  {
    f_EPH[i][0] = 0; f_EPH[i][1] = 0; f_EPH[i][2] = 0;
    f_RNG[i][0] = 0; f_RNG[i][1] = 0; f_RNG[i][2] = 0;
    
    xi_i[i][0] = 0; xi_i[i][1] = 0;  xi_i[i][2] = 0;
    w_i[i][0] = 0; w_i[i][1] = 0; w_i[i][2] = 0;
    
    rho_i[i] = 0;
    T_e_i[i] = 0;
  }
}

void zero_data_gpu(EPH_GPU& eph_gpu)
{
  int threads = 1;
  int blocks = (eph_gpu.nlocal + eph_gpu.nghost + threads - 1) / threads;
  
  zero_data_cu<<<blocks, threads>>>(eph_gpu);
  cudaDeviceSynchronize();
}

/** END ZERO DATA ON DEVICE **/


/** CALCULATE ENVIRONMENT **/

__global__
void calculate_environment_cu(EPH_GPU eph_gpu)
{
  int thread_index = threadIdx.x;
  int block_index = blockIdx.x;
  int block_dimension = blockDim.x;
  
  EPH_Beta_GPU& beta = *((EPH_Beta_GPU*) eph_gpu.beta_gpu);
  
  int nlocal = eph_gpu.nlocal;
  
  double3d* x = eph_gpu.x_gpu;
  
  // TODO: maybe create a shared coordinate pool
  
  int* type = eph_gpu.type_gpu;
  int* mask = eph_gpu.mask_gpu;
  int groupbit = eph_gpu.groupbit;
  
  int* type_map = eph_gpu.type_map_gpu;
  
  double* rho_i = eph_gpu.rho_i_gpu;
  
  double r_cutoff_sq = eph_gpu.r_cutoff_sq;
  
  int* number_neigh = eph_gpu.number_neigh_gpu;
  int* index_neigh = eph_gpu.index_neigh_gpu;
  int* neighs = eph_gpu.neighs_gpu;
  
  int i = block_index;
  double rho_thread = 0;
  
  if(i < nlocal && (mask[i] & groupbit))
  {
    for(int j = thread_index; j < number_neigh[i]; j += block_dimension)
    {
      int jj = neighs[index_neigh[i] + j];
      jj &= eph_gpu.neighmask;
      
      double r_sq = get_distance_sq(x[i], x[jj]);
      
      if(r_sq < r_cutoff_sq)
      { 
        rho_thread += beta.get_rho_r_sq(type_map[type[jj] - 1], r_sq);
      }
    }
  }
  
  atomicAdd(&(rho_i[i]), rho_thread);
}

void calculate_environment_gpu(EPH_GPU& eph_gpu)
{
  int threads = 32;
  int blocks = eph_gpu.nlocal;
  
  calculate_environment_cu<<<blocks, threads>>>(eph_gpu);
  cudaDeviceSynchronize();
}

/** END CALCULATE ENVIRONMENT **/

// TODO: parallelise over x,y,z coordinates
// TODO: use atomicAdd_block() instead of +=

/** FORCE PRL STAGE1 **/

__global__
void force_prl_stage1_cu(EPH_GPU eph_gpu)
{
  int thread_index = threadIdx.x;
  int block_index = blockIdx.x;
  int block_dimension = blockDim.x;
  
  EPH_Beta_GPU& beta = *((EPH_Beta_GPU*) eph_gpu.beta_gpu);
  
  int nlocal = eph_gpu.nlocal;
  
  double3d* x = eph_gpu.x_gpu;
  double3d* v = eph_gpu.v_gpu;
  double3d* w_i = eph_gpu.w_i_gpu;
  
  int* type = eph_gpu.type_gpu;
  int* mask = eph_gpu.mask_gpu;
  int groupbit = eph_gpu.groupbit;
  
  int* type_map = eph_gpu.type_map_gpu;
  
  double* rho_i = eph_gpu.rho_i_gpu;
  
  double r_cutoff_sq = eph_gpu.r_cutoff_sq;
  
  int* number_neigh = eph_gpu.number_neigh_gpu;
  int* index_neigh = eph_gpu.index_neigh_gpu;
  int* neighs = eph_gpu.neighs_gpu;
  
  int i = block_index;
  
  if(i < nlocal && (mask[i] & groupbit) && rho_i[i] > 0)
  {
    double3d w_i_thread = {0, 0, 0};
    
    double alpha_i = beta.get_alpha(type_map[type[i] - 1], rho_i[i]);
    
    for(int j = thread_index; j < number_neigh[i]; j += block_dimension)
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
      w_i_thread[0] += dvar * e_ij[0];
      w_i_thread[1] += dvar * e_ij[1];
      w_i_thread[2] += dvar * e_ij[2];
    }
    
    atomicAdd(&(w_i[i][0]), w_i_thread[0]);
    atomicAdd(&(w_i[i][1]), w_i_thread[1]);
    atomicAdd(&(w_i[i][2]), w_i_thread[2]);
  }
}

// W * v
void force_prl_stage1_gpu(EPH_GPU& eph_gpu)
{
  int threads = 32;
  int blocks = eph_gpu.nlocal;
  
  force_prl_stage1_cu<<<blocks, threads>>>(eph_gpu);
  cudaDeviceSynchronize();
}

/** END FORCE PRL STAGE1 **/

/** FORCE PRL STAGE2 **/

__global__
void force_prl_stage2_cu(EPH_GPU eph_gpu)
{
  int thread_index = threadIdx.x;
  int block_index = blockIdx.x;
  int block_dimension = blockDim.x;
  
  EPH_Beta_GPU& beta = *((EPH_Beta_GPU*) eph_gpu.beta_gpu);
  
  int nlocal = eph_gpu.nlocal;
  
  double3d* x = eph_gpu.x_gpu;
  
  double3d* f_EPH = eph_gpu.f_EPH_gpu;
  double3d* w_i = eph_gpu.w_i_gpu;
  
  int* type = eph_gpu.type_gpu;
  int* mask = eph_gpu.mask_gpu;
  int groupbit = eph_gpu.groupbit;
  
  int* type_map = eph_gpu.type_map_gpu;
  
  double* rho_i = eph_gpu.rho_i_gpu;
  
  double r_cutoff_sq = eph_gpu.r_cutoff_sq;
  
  int* number_neigh = eph_gpu.number_neigh_gpu;
  int* index_neigh = eph_gpu.index_neigh_gpu;
  int* neighs = eph_gpu.neighs_gpu;
  
  int i = block_index;
  
  if(i < nlocal && (mask[i] & groupbit) && rho_i[i] > 0) {
    double alpha_i = beta.get_alpha(type_map[type[i] - 1], rho_i[i]);
    double3d f_EPH_thread = {0, 0, 0};
    
    for(int j = thread_index; j < number_neigh[i]; j += block_dimension)
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
      f_EPH_thread[0] -= dvar * e_ij[0];
      f_EPH_thread[1] -= dvar * e_ij[1];
      f_EPH_thread[2] -= dvar * e_ij[2];
    }
    
    atomicAdd(&(f_EPH[i][0]), f_EPH_thread[0]);
    atomicAdd(&(f_EPH[i][1]), f_EPH_thread[1]);
    atomicAdd(&(f_EPH[i][2]), f_EPH_thread[2]);
  }
}

// W * wi
void force_prl_stage2_gpu(EPH_GPU& eph_gpu)
{
  int threads = 32;
  int blocks = eph_gpu.nlocal;
  
  force_prl_stage2_cu<<<blocks, threads>>>(eph_gpu);
  cudaDeviceSynchronize();
}

/** END FORCE PRL STAGE2 **/

/** FORCE PRL STAGE3 **/

__global__
void force_prl_stage3_cu(EPH_GPU eph_gpu)
{
  int thread_index = threadIdx.x;
  int block_index = blockIdx.x;
  int block_dimension = blockDim.x;
  
  EPH_Beta_GPU& beta = *((EPH_Beta_GPU*) eph_gpu.beta_gpu);
  
  int nlocal = eph_gpu.nlocal;
  
  double3d* x = eph_gpu.x_gpu;
  double3d* xi_i = eph_gpu.xi_i_gpu;
  
  double3d* f_RNG = eph_gpu.f_RNG_gpu;
  
  double* T_e_i = eph_gpu.T_e_i_gpu;
  
  int* type = eph_gpu.type_gpu;
  int* mask = eph_gpu.mask_gpu;
  int groupbit = eph_gpu.groupbit;
  
  int* type_map = eph_gpu.type_map_gpu;
  
  double* rho_i = eph_gpu.rho_i_gpu;
  
  double r_cutoff_sq = eph_gpu.r_cutoff_sq;
  
  int* number_neigh = eph_gpu.number_neigh_gpu;
  int* index_neigh = eph_gpu.index_neigh_gpu;
  int* neighs = eph_gpu.neighs_gpu;
  
  int i = block_index;
  
  if(i < nlocal && (mask[i] & groupbit) && rho_i[i] > 0)
  {
    double alpha_i = beta.get_alpha(type_map[type[i] - 1], rho_i[i]);
    
    double3d f_RNG_thread = {0, 0, 0};
    double var = eph_gpu.eta_factor * sqrt(T_e_i[i]);
    
    for(size_t j = thread_index; j < number_neigh[i]; j += block_dimension) 
    {
      int jj = neighs[index_neigh[i] + j];
      jj &= eph_gpu.neighmask;
      
      // calculate the e_ij vector
      double3d e_ij;
      double e_r_sq = get_difference_sq(x[jj], x[i], e_ij);
      
      if((e_r_sq > r_cutoff_sq) or not(rho_i[jj] > 0)) continue;
      
      double alpha_j = beta.get_alpha(type_map[type[jj] - 1], rho_i[jj]);
      
      double v_rho_ji = beta.get_rho_r_sq(type_map[type[jj] - 1], e_r_sq);
      double e_v_xi1 = get_scalar(e_ij, xi_i[i]);
      double var1 = alpha_i * v_rho_ji * e_v_xi1 / (rho_i[i] * e_r_sq);
      
      double v_rho_ij = beta.get_rho_r_sq(type_map[type[i] - 1], e_r_sq);
      double e_v_xi2 = get_scalar(e_ij, xi_i[jj]);
      double var2 = alpha_j * v_rho_ij * e_v_xi2 / (rho_i[jj] * e_r_sq);
      
      double dvar = var1 - var2;
      f_RNG_thread[0] += dvar * e_ij[0] * var;
      f_RNG_thread[1] += dvar * e_ij[1] * var;
      f_RNG_thread[2] += dvar * e_ij[2] * var;
    }
    
    atomicAdd(&(f_RNG[i][0]), f_RNG_thread[0]);
    atomicAdd(&(f_RNG[i][1]), f_RNG_thread[1]);
    atomicAdd(&(f_RNG[i][2]), f_RNG_thread[2]);
  }
}

// W * xi
void force_prl_stage3_gpu(EPH_GPU& eph_gpu)
{
  int threads = 32;
  int blocks = eph_gpu.nlocal;
  
  force_prl_stage3_cu<<<blocks, threads>>>(eph_gpu);
  cudaDeviceSynchronize();
}

/** END FORCE PRL STAGE3 **/
