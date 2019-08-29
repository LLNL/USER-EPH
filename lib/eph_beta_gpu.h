#ifdef FIX_EPH_GPU
#ifndef EPH_BETA_GPU_H
#define EPH_BETA_GPU_H

#include <vector>

#include <iostream>

#include "eph_spline_gpu.h"
#include "../eph_beta.h"

// beta_rho class
using Beta = EPH_Beta<double, std::allocator, std::vector>;

class EPH_Beta_GPU : public Beta 
{
  public:
    EPH_Beta_GPU() : 
      Beta(),
      rho_gpu_device {nullptr},
      rho_r_sq_gpu_device {nullptr},
      alpha_gpu_device {nullptr},
      beta_gpu_device {nullptr} {}
      
    EPH_Beta_GPU(const Beta& beta_input) :
      Beta(beta_input)
    {
      //std::cout << "EPH_Beta_GPU constructor from EPH_Beta\n";
      n_elements = rho.size();
      
      for(size_t i = 0; i != n_elements; ++i)
      {
        rho_gpu.push_back(EPH_Spline_GPU(rho[i]));
        rho_r_sq_gpu.push_back(EPH_Spline_GPU(rho_r_sq[i]));
        alpha_gpu.push_back(EPH_Spline_GPU(alpha[i]));
        beta_gpu.push_back(EPH_Spline_GPU(beta[i]));
      }
      
      cudaMalloc((void**) &rho_gpu_device, n_elements * sizeof(EPH_Spline_GPU));
      cudaMemcpy(rho_gpu_device, rho_gpu.data(), n_elements * sizeof(EPH_Spline_GPU), cudaMemcpyHostToDevice);
      
      cudaMalloc((void**) &rho_r_sq_gpu_device, n_elements * sizeof(EPH_Spline_GPU));
      cudaMemcpy(rho_r_sq_gpu_device, rho_r_sq_gpu.data(), n_elements * sizeof(EPH_Spline_GPU), cudaMemcpyHostToDevice);
      
      cudaMalloc((void**) &alpha_gpu_device, n_elements * sizeof(EPH_Spline_GPU));
      cudaMemcpy(alpha_gpu_device, alpha_gpu.data(), n_elements * sizeof(EPH_Spline_GPU), cudaMemcpyHostToDevice);
      
      cudaMalloc((void**) &beta_gpu_device, n_elements * sizeof(EPH_Spline_GPU));
      cudaMemcpy(beta_gpu_device, beta_gpu.data(), n_elements * sizeof(EPH_Spline_GPU), cudaMemcpyHostToDevice);
    }
    
    EPH_Beta_GPU(const EPH_Beta_GPU &beta_input) 
      : EPH_Beta_GPU((Beta) beta_input) {}
    
    ~EPH_Beta_GPU()
    {
      // deallocate rho
      if(rho_gpu_device != nullptr) cudaFree(rho_gpu_device);
      
      // deallocate rho_r_sq
      if(rho_r_sq_gpu_device != nullptr) cudaFree(rho_r_sq_gpu_device);
  
      // deallocate alpha
      if(alpha_gpu_device != nullptr) cudaFree(alpha_gpu_device);
  
      // deallocate beta
      if(beta_gpu_device != nullptr) cudaFree(beta_gpu_device);
    }
    
    __device__ double get_rho(size_t index, double r) 
    {
      #ifdef __CUDA_ARCH__
      return rho_gpu_device[index](r);
      #else
      return rho[index](r);
      #endif
    }
    
    __device__ double get_rho_r_sq(size_t index, double r_sq)
    {
      #ifdef __CUDA_ARCH__
      return rho_r_sq_gpu_device[index](r_sq);
      #else
      return rho_r_sq[index](r_sq);
      #endif
    }
    
    __device__ double get_beta(size_t index, double rho_i)
    {
      #ifdef __CUDA_ARCH__
      return beta_gpu_device[index](rho_i);
      #else
      return beta[index](rho_i);
      #endif
    }
    
    __device__ double get_alpha(size_t index, double rho_i)
    {
      #ifdef __CUDA_ARCH__
      return alpha_gpu_device[index](rho_i);
      #else
      return alpha[index](rho_i);
      #endif
    }
  
  private:
    std::vector<EPH_Spline_GPU> rho_gpu; // gpu splines
    std::vector<EPH_Spline_GPU> rho_r_sq_gpu; // gpu splines
    std::vector<EPH_Spline_GPU> alpha_gpu; // gpu splines
    std::vector<EPH_Spline_GPU> beta_gpu; // gpu splines

    EPH_Spline_GPU* rho_gpu_device; // pointer to spline array on gpu
    EPH_Spline_GPU* rho_r_sq_gpu_device;
    EPH_Spline_GPU* alpha_gpu_device;
    EPH_Spline_GPU* beta_gpu_device;
    
    // TODO: to be removed
    void clean_memory()
    {
      // deallocate rho
      if(rho_gpu_device != nullptr) cudaFree(rho_gpu_device);
      
      // deallocate rho_r_sq
      if(rho_r_sq_gpu_device != nullptr) cudaFree(rho_r_sq_gpu_device);
  
      // deallocate alpha
      if(alpha_gpu_device != nullptr) cudaFree(alpha_gpu_device);
  
      // deallocate beta
      if(beta_gpu_device != nullptr) cudaFree(beta_gpu_device);
    }
    
    void allocate_and_copy()
    {
      n_elements = rho.size();
      
      for(size_t i = 0; i != n_elements; ++i)
      {
        rho_gpu.push_back(EPH_Spline_GPU(rho[i]));
        rho_r_sq_gpu.push_back(EPH_Spline_GPU(rho_r_sq[i]));
        alpha_gpu.push_back(EPH_Spline_GPU(alpha[i]));
        beta_gpu.push_back(EPH_Spline_GPU(beta[i]));
      }
      
      cudaMalloc((void**) &rho_gpu_device, n_elements * sizeof(EPH_Spline_GPU));
      cudaMemcpy(rho_gpu_device, rho_gpu.data(), n_elements * sizeof(EPH_Spline_GPU), cudaMemcpyHostToDevice);
      
      cudaMalloc((void**) &rho_r_sq_gpu_device, n_elements * sizeof(EPH_Spline_GPU));
      cudaMemcpy(rho_r_sq_gpu_device, rho_r_sq_gpu.data(), n_elements * sizeof(EPH_Spline_GPU), cudaMemcpyHostToDevice);
      
      cudaMalloc((void**) &alpha_gpu_device, n_elements * sizeof(EPH_Spline_GPU));
      cudaMemcpy(alpha_gpu_device, alpha_gpu.data(), n_elements * sizeof(EPH_Spline_GPU), cudaMemcpyHostToDevice);
      
      cudaMalloc((void**) &beta_gpu_device, n_elements * sizeof(EPH_Spline_GPU));
      cudaMemcpy(beta_gpu_device, beta_gpu.data(), n_elements * sizeof(EPH_Spline_GPU), cudaMemcpyHostToDevice);
    }
};

// TODO: remove later
void test_beta_rho_gpu(EPH_Beta_GPU &beta);

#endif
#endif
