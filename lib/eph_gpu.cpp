
#include <cuda.h>
#include <cuda_runtime_api.h>

#include <iostream>

#define FIX_EPH_GPU
#include "eph_gpu.h"

void calculate_environment_gpu();

// TODO: consider moving class member functions to .h instead

// EPH_Spline_GPU
void EPH_Spline_GPU::clean_memory()
{
  if(c_gpu != nullptr)
  {
    cudaFree(c_gpu);
    c_gpu = nullptr;
    n_points = 0;
  }
}

void EPH_Spline_GPU::allocate_and_copy(size_t n)
{
  clean_memory();
  
  n_points = n;
  cudaMalloc((void**) &c_gpu, n_points * sizeof(Coefficients));
  cudaMemcpy(c_gpu, c.data(), n_points * sizeof(Coefficients), cudaMemcpyHostToDevice);
}

void EPH_Spline_GPU::allocate_and_copy(size_t n, const Coefficients *c_device)
{
  clean_memory();
  
  n_points = n;
  cudaMalloc((void**) &c_gpu, n_points * sizeof(Coefficients));
  cudaMemcpy(c_gpu, c_device, n_points * sizeof(Coefficients), cudaMemcpyDeviceToDevice);
}

EPH_Spline_GPU::EPH_Spline_GPU() 
{
  //std::cout << "Empty constructor\n";
  n_points = 0;
  c_gpu = nullptr;
}

// TODO: remove this as this will not be used 
EPH_Spline_GPU::EPH_Spline_GPU(double dx, std::vector<double> &y) :
  Spline(dx, y)
{
  //std::cout << "Vector allocation constructor\n";
  // TODO: error checks  
  allocate_and_copy(c.size());
  // cleanup unused memory
  c.resize(0);
}

EPH_Spline_GPU::EPH_Spline_GPU(const Spline& spl) :
  Spline(spl)
{
  //std::cout << "Construct from spline\n";  
  allocate_and_copy(c.size());
  // cleanup unused memory
  c.resize(0);
}

EPH_Spline_GPU::EPH_Spline_GPU(const EPH_Spline_GPU& spl)
{
  //std::cout << "EPH_Spline_GPU copy construct\n";
  inv_dx = spl.inv_dx;
  allocate_and_copy(spl.n_points, spl.c_gpu);
  c.resize(0);
}

EPH_Spline_GPU::EPH_Spline_GPU(EPH_Spline_GPU&& spl)
{
  //std::cout << "EPH_Spline_GPU move construct\n";
  n_points = spl.n_points;
  inv_dx = spl.inv_dx;
  c_gpu = spl.c_gpu;
  spl.c_gpu = nullptr;
  c.resize(0);
}

EPH_Spline_GPU::~EPH_Spline_GPU()
{
  //std::cout << "Destroy EPH_Spline_GPU\n";
  clean_memory();
}

EPH_Spline_GPU& EPH_Spline_GPU::operator=(const EPH_Spline_GPU& spl)
{
  //std::cout << "Copy assignment\n";
  clean_memory();
  
  inv_dx = spl.inv_dx;
  allocate_and_copy(spl.n_points, spl.c_gpu);
  
  return *this;
}

EPH_Spline_GPU& EPH_Spline_GPU::operator=(const EPH_Spline_GPU&& spl)
{
  //std::cout << "Move assignment\n";
  clean_memory();
  
  n_points = spl.n_points;
  inv_dx = spl.inv_dx;
  c_gpu = spl.c_gpu;
  
  return *this;
}

// EPH_Beta_GPU
void EPH_Beta_GPU::clean_memory()
{
  // deallocate rho
  if(rho_gpu != nullptr)
  {
    delete[] rho_gpu;
    rho_gpu = nullptr;
  }
  
  if(rho_gpu_device != nullptr)
  {
    cudaFree(rho_gpu_device);
    rho_gpu_device = nullptr;
  }
  
  // deallocate rho_r_sq
  if(rho_r_sq_gpu != nullptr)
  {
    delete[] rho_r_sq_gpu;
    rho_r_sq_gpu = nullptr;
  }
  
  if(rho_r_sq_gpu_device != nullptr)
  {
    cudaFree(rho_r_sq_gpu_device);
    rho_r_sq_gpu_device = nullptr;
  }
  
  // deallocate alpha
  if(alpha_gpu != nullptr)
  {
    delete[] alpha_gpu;
    alpha_gpu = nullptr;
  }
  
  if(alpha_gpu_device != nullptr)
  {
    cudaFree(alpha_gpu_device);
    alpha_gpu_device = nullptr;
  }
  
   // deallocate beta
  if(beta_gpu != nullptr)
  {
    delete[] beta_gpu;
    beta_gpu = nullptr;
  }
  
  if(beta_gpu_device != nullptr)
  {
    cudaFree(beta_gpu_device);
    beta_gpu_device = nullptr;
  }
}

void EPH_Beta_GPU::allocate_and_copy(
  EPH_Spline_GPU** dst_device, EPH_Spline_GPU** dst, 
  const std::vector<Spline> &src)
{
  cudaMalloc((void**) dst_device, n_elements * sizeof(EPH_Spline_GPU));
  (*dst) = new EPH_Spline_GPU[n_elements];
  
  for(size_t i = 0; i != n_elements; ++i) 
  {
    (*dst)[i] = EPH_Spline_GPU(src[i]);
  }
  
  cudaMemcpy((*dst_device), (*dst), n_elements * sizeof(EPH_Spline_GPU), cudaMemcpyHostToDevice);
}

void EPH_Beta_GPU::allocate_and_copy()
{
  allocate_and_copy(&rho_gpu_device, &rho_gpu, rho);
  allocate_and_copy(&rho_r_sq_gpu_device, &rho_r_sq_gpu, rho_r_sq);
  allocate_and_copy(&alpha_gpu_device, &alpha_gpu, alpha);
  allocate_and_copy(&beta_gpu_device, &beta_gpu, beta);
}

void EPH_Beta_GPU::allocate_and_copy(
  EPH_Spline_GPU** dst_device, EPH_Spline_GPU** dst, 
  const EPH_Spline_GPU* src_device, const EPH_Spline_GPU* src
)
{
  cudaMalloc((void**) dst_device, n_elements * sizeof(EPH_Spline_GPU));
  (*dst) = new EPH_Spline_GPU[n_elements];
  
  for(size_t i = 0; i != n_elements; ++i) 
  {
    (*dst)[i] = src[i];
  }
  
  cudaMemcpy((*dst_device), (*dst), n_elements * sizeof(EPH_Spline_GPU), cudaMemcpyHostToDevice);
}

EPH_Beta_GPU::EPH_Beta_GPU()
{
  rho_gpu = nullptr;
  rho_gpu_device = nullptr;
  
  rho_r_sq_gpu = nullptr;
  rho_r_sq_gpu_device = nullptr;
  
  alpha_gpu = nullptr;
  alpha_gpu_device = nullptr;
  
  beta_gpu = nullptr;
  beta_gpu_device = nullptr;
}

EPH_Beta_GPU::EPH_Beta_GPU(Beta& beta) :
  Beta(beta)
{
  std::cout << "EPH_Beta_GPU constructor from EPH_Beta\n";
  allocate_and_copy();
  rho.resize(0);
}

EPH_Beta_GPU::EPH_Beta_GPU(const EPH_Beta_GPU &beta)
{
  std::cout << "EPH_Beta_GPU copy constructor\n";
  
  n_elements = beta.n_elements;
  
  allocate_and_copy(&rho_gpu_device, &rho_gpu, 
    beta.rho_gpu_device, beta.rho_gpu);
    
  allocate_and_copy(&rho_r_sq_gpu_device, &rho_r_sq_gpu, 
    beta.rho_r_sq_gpu_device, beta.rho_r_sq_gpu);
  
  allocate_and_copy(&alpha_gpu_device, &alpha_gpu, 
    beta.alpha_gpu_device, beta.alpha_gpu);
  
  allocate_and_copy(&beta_gpu_device, &beta_gpu, 
    beta.beta_gpu_device, beta.beta_gpu);
}

// TODO: this is not working correctly
EPH_Beta_GPU::~EPH_Beta_GPU() 
{
  clean_memory();
}

EPH_Beta_GPU& EPH_Beta_GPU::operator=(const EPH_Beta_GPU& beta)
{
  clean_memory();
}
