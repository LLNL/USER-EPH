
#ifndef EPH_GPU_CUH
#define EPH_GPU_CUH

#define FIX_EPH_GPU

// external headers
#include <vector>

// internal headers
#include "eph_gpu.h"

#include "../eph_spline.h"
#include "../eph_beta.h"

// TODO: remove this file

__device__ double EPH_Spline_GPU::operator() (double x)
{
  int index = x * inv_dx;
  return c_gpu[index].a + x * (c_gpu[index].b + x * (c_gpu[index].c + x * c_gpu[index].d));
}



#endif
