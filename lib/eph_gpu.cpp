
#include <cuda.h>
#include <cuda_runtime_api.h>

#include <iostream>

#define FIX_EPH_GPU
#include "eph_gpu.h"
#include "eph_beta_gpu.h"

void calculate_environment_gpu(EPH_GPU& eph_gpu);

// EPH_GPU functions

void EPH_GPU::grow(size_t ngrow)
{
  std::cout << "GROWING\n";
  if(n)
  {
    cudaFree(type_gpu);
    cudaFree(mask_gpu);
    
    cudaFree(x_gpu);
    cudaFree(v_gpu);
    cudaFree(f_gpu);
    
    cudaFree(f_EPH_gpu);
    cudaFree(f_RNG_gpu);
    
    cudaFree(xi_i_gpu);
    cudaFree(w_i_gpu);
    
    cudaFree(rho_i_gpu);
  }
  
  n = ngrow;
  
  cudaMalloc((void**) &type_gpu, n * sizeof(int));
  cudaMalloc((void**) &mask_gpu, n * sizeof(int));
  
  cudaMalloc((void**) &x_gpu, n * sizeof(double3d));
  cudaMalloc((void**) &v_gpu, n * sizeof(double3d));
  cudaMalloc((void**) &f_gpu, n * sizeof(double3d));
  
  cudaMalloc((void**) &f_EPH_gpu, n * sizeof(double3d));
  cudaMalloc((void**) &f_RNG_gpu, n * sizeof(double3d));
  
  cudaMalloc((void**) &xi_i_gpu, n * sizeof(double3d));
  cudaMalloc((void**) &w_i_gpu, n * sizeof(double3d));
  
  cudaMalloc((void**) &rho_i_gpu, n * sizeof(double));
}

void cpu_to_device_EPH_GPU(void* dst, void* src, size_t n);
void device_to_cpu_EPH_GPU(void* dst, void* src, size_t n);

void synchronise_EPH_GPU();


EPH_GPU allocate_EPH_GPU(Beta& beta_in, int types, int* type_map)
{
  EPH_GPU eph_gpu;
  
  EPH_Beta_GPU* beta = new EPH_Beta_GPU(beta_in);
  eph_gpu.beta = (void*) beta;
  
  cudaMalloc((void**) &(eph_gpu.beta_gpu), sizeof(EPH_Beta_GPU));
  cudaMemcpy(eph_gpu.beta_gpu, beta, sizeof(EPH_Beta_GPU), cudaMemcpyHostToDevice);
  
  eph_gpu.types = types;
  cudaMalloc((void**) &(eph_gpu.type_map_gpu), types * sizeof(int));
  cudaMemcpy(eph_gpu.type_map_gpu, type_map, types * sizeof(int), cudaMemcpyHostToDevice);
  
  eph_gpu.r_cutoff = beta_in.get_r_cutoff();
  eph_gpu.r_cutoff_sq = beta_in.get_r_cutoff_sq();
  eph_gpu.rho_cutoff = beta_in.get_rho_cutoff();
  
  eph_gpu.nlocal = 0;
  eph_gpu.nghost = 0;
  eph_gpu.n = 0;
  
  eph_gpu.type_gpu = nullptr;
  eph_gpu.mask_gpu = nullptr;
  
  eph_gpu.x_gpu = nullptr;
  eph_gpu.v_gpu = nullptr;
  eph_gpu.f_gpu = nullptr;
  
  eph_gpu.f_EPH_gpu = nullptr;
  eph_gpu.f_RNG_gpu = nullptr;
  
  eph_gpu.xi_i_gpu = nullptr;
  eph_gpu.w_i_gpu = nullptr;
  
  eph_gpu.rho_i_gpu = nullptr;
  
  return eph_gpu;
}

void deallocate_EPH_GPU(EPH_GPU& eph_gpu)
{
  cudaFree(eph_gpu.beta_gpu);
  EPH_Beta_GPU* beta = (EPH_Beta_GPU*) eph_gpu.beta;
  delete beta;
  
  cudaFree(eph_gpu.type_map_gpu);
  
  if(eph_gpu.n)
  {
    cudaFree(eph_gpu.type_gpu);
    cudaFree(eph_gpu.mask_gpu);
    
    cudaFree(eph_gpu.x_gpu);
    cudaFree(eph_gpu.v_gpu);
    cudaFree(eph_gpu.f_gpu);
    
    cudaFree(eph_gpu.f_EPH_gpu);
    cudaFree(eph_gpu.f_RNG_gpu);
    
    cudaFree(eph_gpu.xi_i_gpu);
    cudaFree(eph_gpu.w_i_gpu);
    
    cudaFree(eph_gpu.rho_i_gpu);
  }
}

void cpu_to_device_EPH_GPU(void* dst, void* src, size_t n)
{
  cudaMemcpy(dst, src, n, cudaMemcpyHostToDevice);
}

void device_to_cpu_EPH_GPU(void* dst, void* src, size_t n)
{
  cudaMemcpy(dst, src, n, cudaMemcpyDeviceToHost);
}
