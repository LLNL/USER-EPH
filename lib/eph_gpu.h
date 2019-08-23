
#ifdef FIX_EPH_GPU

#ifndef EPH_GPU_H
#define EPH_GPU_H

// external headers
#include <vector>

// internal headers
#include "../eph_spline.h"

void call_dummy(int, int);

using Spline = EPH_Spline<double, std::allocator, std::vector>;

// spline interpolator
class EPH_Spline_GPU : public Spline
{
  public:
    EPH_Spline_GPU(double dx, std::vector<double> &y);
    ~EPH_Spline_GPU();
    
    double operator() (double x) 
    {
      int index = x * (*inv_dx_gpu);
      return c[index].a + x * (c[index].b + x * (c[index].c + x * c[index].d));
    }
  
  private:
    
    double *inv_dx_gpu;
    Coefficients *c_gpu;
};

#endif
#endif
