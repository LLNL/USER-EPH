

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
  
  std::cout << "### Testing Spline ###\n";
  vector<double> y {0,1,2,3,4};
  Spline spl(1.0, y);
  SplineGPU splgpu(spl);
  
  printf("CPU: %.2f %.2f %.2f %.2f\n", spl(0.0), spl(1.0), spl(1.5), spl(2.0));
  test_interpolator_gpu(splgpu);
  test_interpolator_gpu(splgpu);
  
  std::cout << "### Testing Beta ###\n";
  
  Beta beta("Beta_Rho.beta");
  BetaGPU betagpu(beta); 

  printf("CPU: %.3f %.3f %.3f %.3f\n", 
    beta.get_rho(0, 1.0), beta.get_rho_r_sq(0, 1.0), 
    beta.get_alpha(0, 0.1), beta.get_beta(0, 0.1));
  test_beta_rho_gpu(betagpu);
  test_beta_rho_gpu(betagpu);
  
  std::cout << "### Testing EPH_GPU ###\n";
  int types = 3;
  int type_map[3] = {0, 1, 0};
  
  EPH_GPU eph_gpu = allocate_EPH_GPU(beta, types, type_map);
  calculate_environment_gpu(eph_gpu);
  
  MPI_Finalize();
  return 0;
}

