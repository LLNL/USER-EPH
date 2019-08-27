
#ifdef FIX_EPH_GPU

#ifndef EPH_GPU_H
#define EPH_GPU_H

// external headers
#include <memory>
#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

// internal headers
#include "../eph_spline.h"
#include "../eph_beta.h"

// spline interpolator; defined and build elsewhere
using Spline = EPH_Spline<double, std::allocator, std::vector>;

// spline interpolator
class EPH_Spline_GPU : public Spline
{
  public:
    EPH_Spline_GPU();
    EPH_Spline_GPU(double dx, std::vector<double> &y);
    EPH_Spline_GPU(const Spline& spl);
    EPH_Spline_GPU(const EPH_Spline_GPU& spl);
    EPH_Spline_GPU(EPH_Spline_GPU&& spl);
    ~EPH_Spline_GPU();
    
    EPH_Spline_GPU& operator=(const EPH_Spline_GPU& spl);
    EPH_Spline_GPU& operator=(const EPH_Spline_GPU&& spl);
    
    __device__ double operator() (double x)
    {
      #ifdef __CUDA_ARCH__
      int index = x * inv_dx;
      return c_gpu[index].a + x * (c_gpu[index].b + x * (c_gpu[index].c + x * c_gpu[index].d));
      #else
      return 0;
      #endif
    }
  
  private:
    size_t n_points;
    Coefficients *c_gpu;
    
    void clean_memory();
    void allocate_and_copy(size_t n);
    void allocate_and_copy(size_t n, const Coefficients *c_device);
};

// beta_rho class
using Beta = EPH_Beta<double, std::allocator, std::vector>;

class EPH_Beta_GPU : public Beta {
  public:
    EPH_Beta_GPU();
    EPH_Beta_GPU(Beta& beta);
    EPH_Beta_GPU(const EPH_Beta_GPU &beta);
    ~EPH_Beta_GPU();
    
    EPH_Beta_GPU& operator=(const EPH_Beta_GPU& beta);
    
    __device__ double get_rho(size_t index, double r) 
    {
      return rho_gpu_device[index](r);
    }
    
    __device__ double get_rho_r_sq(size_t index, double r_sq)
    {
      return rho_r_sq_gpu_device[index](r_sq);
    }
    
    __device__ double get_beta(size_t index, double rho_i)
    {
      return beta_gpu_device[index](rho_i);
    }
    
    __device__ double get_alpha(size_t index, double rho_i)
    {
      return alpha_gpu_device[index](rho_i);
    }
  
  private:
    EPH_Spline_GPU* rho_gpu; // gpu splines
    EPH_Spline_GPU* rho_gpu_device; // pointer to spline array on gpu
    
    EPH_Spline_GPU* rho_r_sq_gpu;
    EPH_Spline_GPU* rho_r_sq_gpu_device;
    
    EPH_Spline_GPU* alpha_gpu;
    EPH_Spline_GPU* alpha_gpu_device;
    
    EPH_Spline_GPU* beta_gpu;
    EPH_Spline_GPU* beta_gpu_device;
    
    void clean_memory();
    void allocate_and_copy();
    void allocate_and_copy(
      EPH_Spline_GPU** dst_device, EPH_Spline_GPU** dst,
      const std::vector<Spline> &src);
    void allocate_and_copy(
      EPH_Spline_GPU** dst_device, EPH_Spline_GPU** dst, 
      const EPH_Spline_GPU* src_device, const EPH_Spline_GPU* src);
};

// TODO: remove these
void test_interpolator_gpu(EPH_Spline_GPU &spl);
void test_beta_rho_gpu(EPH_Beta_GPU &beta);

#endif
#endif
