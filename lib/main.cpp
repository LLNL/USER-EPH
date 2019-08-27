

#include "eph_gpu.h"

#include <memory>
#include <vector>
#include <iostream>
#include <cstdio>

#include <mpi.h>

using namespace std;

using Spline = EPH_Spline<double, allocator, vector>;
using SplineGPU = EPH_Spline_GPU;

int main(int args, char **argv) 
{
  MPI_Init(&args, &argv);

  int myID;
  int nrPS;
  
  MPI_Comm_size(MPI_COMM_WORLD, &nrPS);
  MPI_Comm_rank(MPI_COMM_WORLD, &myID);
  
  vector<double> y {0,1,2,3,4};
  
  Spline spl(1.0, y);
  Spline spl2(spl);
  SplineGPU splgpu(1.0, y);
  SplineGPU splgpu2(spl);
  
  printf("CPU: %.2f %.2f %.2f %.2f\n", spl(0.0), spl(1.0), spl(1.5), spl(2.0));
  printf("CPU: %.2f %.2f %.2f %.2f\n", spl2(0.0), spl2(1.0), spl2(1.5), spl2(2.0));
  printf("CPU: %.2f %.2f %.2f %.2f\n", splgpu(0.0), splgpu(1.0), splgpu(1.5), splgpu(2.0));  
  test_interpolator_gpu(splgpu, 0, 0);
  test_interpolator_gpu(splgpu2, 0, 0);
  
  MPI_Finalize();
  return 0;
}

