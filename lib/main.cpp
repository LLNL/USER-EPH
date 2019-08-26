

#include "eph_gpu.h"

#include <memory>
#include <vector>
#include <iostream>

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

  call_dummy(myID, nrPS);
  
  vector<double> y {0,1,2,3,4};
  
  Spline spl(1.0, y);
  SplineGPU splgpu(1.0, y);
  
  cout << spl(0.0) << ' ' << spl(1.0) << ' ' << spl(1.5) << ' ' << spl(2.0) << '\n';
  
  test_interpolator_gpu(splgpu, 0, 0);
  
  MPI_Finalize();
  return 0;
}

