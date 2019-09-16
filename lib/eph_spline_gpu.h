
#ifdef FIX_EPH_GPU
#ifndef EPH_SPLINE_GPU_H
#define EPH_SPLINE_GPU_H

// external headers
#include <cuda.h>
#include <cuda_runtime.h>

#include <iostream>

// internal headers
#include "eph_spline.h"

class EPH_Spline_GPU : public Spline
{
  public:
    EPH_Spline_GPU() : 
      Spline(),
      n_points {0},
      c_gpu {nullptr} 
    { }
      
    EPH_Spline_GPU(const Spline& spl) :
      Spline(spl)
    {
      n_points = c.size();
      cudaMalloc((void**) &c_gpu, n_points * sizeof(Coefficients));
      cudaMemcpy(c_gpu, c.data(), n_points * sizeof(Coefficients), cudaMemcpyHostToDevice);
    }
    
    EPH_Spline_GPU(const EPH_Spline_GPU& spl)
      : EPH_Spline_GPU((Spline) spl) 
    { }
    
    
    ~EPH_Spline_GPU()
    {
      if(c_gpu != nullptr) cudaFree(c_gpu);
    }
    
    __device__ double operator() (double x)
    {
      int index = x * inv_dx;
      #ifdef __CUDA_ARCH__
      return c_gpu[index].a + x * (c_gpu[index].b + x * (c_gpu[index].c + x * c_gpu[index].d));
      #else
      return c[index].a + x * (c[index].b + x * (c[index].c + x * c[index].d));;
      #endif
    }
  
  private:
    size_t n_points;
    Coefficients *c_gpu;
    
    // TODO: to be removed
    void clean_memory()
    {
      if(c_gpu != nullptr)
      {
        cudaFree(c_gpu);
        c_gpu = nullptr;
        n_points = 0;
      }
    }
    
    void allocate_and_copy(size_t n)
    {
      n_points = n;
      cudaMalloc((void**) &c_gpu, n_points * sizeof(Coefficients));
      cudaMemcpy(c_gpu, c.data(), n_points * sizeof(Coefficients), cudaMemcpyHostToDevice);
    }
};

#endif
#endif
