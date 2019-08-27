
#ifdef FIX_EPH_GPU

#ifndef EPH_GPU_H
#define EPH_GPU_H

// external headers
#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

// internal headers
#include "../eph_spline.h"
#include "../eph_beta.h"

using Spline = EPH_Spline<double, std::allocator, std::vector>;

// spline interpolator; defined and build elsewhere
using Spline = EPH_Spline<double, std::allocator, std::vector>;

// spline interpolator
class EPH_Spline_GPU : public Spline
{
  public:
    EPH_Spline_GPU(double dx, std::vector<double> &y);
    EPH_Spline_GPU(Spline& spl);
    ~EPH_Spline_GPU();
    
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
    Coefficients *c_gpu;
};

// beta_rho class
/*
using Beta = EPH_Beta<double, std::allocator, std::vector>;

class EPH_Beta_GPU : public Beta {
  public:
    EPH_Beta_GPU() {};
    ~EPH_Beta_GPU() {};
    size_t get_n_elements();
    Float get_r_cutoff();
    Float get_r_cutoff_sq();
    Float get_rho_cutoff();
    uint8_t get_element_number(size_t index)
    std::string get_element_name(size_t index);
    Float get_rho(size_t index, Float r);
    Float get_rho_r_sq(size_t index, Float r_sq);
    Float get_beta(size_t index, Float rho_i);
    Float get_alpha(size_t index, Float rho_i);
  
  private:
    static constexpr unsigned int max_line_length = 1024; // this is for parsing
    
    Float r_cutoff; // cutoff for locality
    Float r_cutoff_sq; // cutoff sq for locality mostly unused
    Float rho_cutoff; // cutoff for largest site density
    size_t n_elements; // number of elements
    
    Container<uint8_t, Allocator<uint8_t>> element_number;
    Container<std::string, Allocator<std::string>> element_name;
    Container<Spline, Allocator<Spline>> rho;
    Container<Spline, Allocator<Spline>> rho_r_sq;
    Container<Spline, Allocator<Spline>> alpha;
    Container<Spline, Allocator<Spline>> beta;
};
*/



// TODO: remove these
void call_dummy(int, int);
void test_interpolator_gpu(EPH_Spline_GPU spl, double *values, int n_values);


#endif
#endif
