

#include "eph_spline.h"
#include "eph_beta.h"

#include "eph_spline_gpu.h"
#include "eph_beta_gpu.h"
#include "eph_gpu.h"

#include <memory>
#include <vector>
#include <iostream>
#include <cstdio>

#include <mpi.h>

using namespace std;

using SplineGPU = EPH_Spline_GPU;
using BetaGPU = EPH_Beta_GPU;

int main(int args, char **argv) 
{
  MPI_Init(&args, &argv);

  int myID;
  int nrPS;
  
  MPI_Comm_size(MPI_COMM_WORLD, &nrPS);
  MPI_Comm_rank(MPI_COMM_WORLD, &myID);
  
  Beta beta("Beta_Rho.beta");
  
  int types = 1;
  int type_map[] = {0};
  
  EPH_GPU eph_gpu = allocate_EPH_GPU(beta, types, type_map);
  
  deallocate_EPH_GPU(eph_gpu);
  
  MPI_Finalize();
  return 0;
}

