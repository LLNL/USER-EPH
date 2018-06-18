
#include <iostream>
#include <cmath>
#include <mpi.h>

#include "eph_fdm.h"

/* 
 * This example will create an electronic system and heat it from one end with
 * a laser source.
 */

// grid size
constexpr unsigned int n_x {4000};
constexpr unsigned int n_y {1};
constexpr unsigned int n_z {1};
constexpr unsigned int n_total {n_x*n_y*n_z};

// grid extent
constexpr double x_0 { 0.0};
constexpr double x_1 { 2000.0}; // 1000 Ang
constexpr double y_0 {-1.0};
constexpr double y_1 { 1.0};
constexpr double z_0 {-1.0};
constexpr double z_1 { 1.0};
constexpr double d_x {(x_1-x_0)/n_x};
constexpr double d_y {(y_1-y_0)/n_y};
constexpr double d_z {(z_1-z_0)/n_z};
constexpr double d_V {d_x*d_y*d_z};

constexpr double c_e {1.0};

constexpr double T0 {1.0};
constexpr double T1 {2.0};

constexpr double kA { 0.1};
constexpr double kB { 1.0};
constexpr double kC {10.0};

// electrons as an FDM system
EPH_FDM sysA {n_x, n_y, n_z};
EPH_FDM sysB {n_x, n_y, n_z};
EPH_FDM sysC {n_x, n_y, n_z};

int main(int args, char **argv) {
  // do the MPI_Init
  int my_id;
  int nr_ps;
  
  MPI_Init(&args, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &nr_ps);
  
  sysA.setComm(MPI_COMM_WORLD, my_id, nr_ps);
  sysA.setBox(x_0, x_1, y_0, y_1, z_0, z_1);
  sysA.setConstants(T0, c_e, 1.0, kA);
  sysA.setSteps(1);
  sysA.setDt(1.0);
  sysA.setT(0u, 0u, 0u, T1); 
  
  sysB.setComm(MPI_COMM_WORLD, my_id, nr_ps);
  sysB.setBox(x_0, x_1, y_0, y_1, z_0, z_1);
  sysB.setConstants(T0, c_e, 1.0, kB);
  sysB.setSteps(1);
  sysB.setDt(1.0);
  sysB.setT(0u, 0u, 0u, T1); 
  
  sysC.setComm(MPI_COMM_WORLD, my_id, nr_ps);
  sysC.setBox(x_0, x_1, y_0, y_1, z_0, z_1);
  sysC.setConstants(T0, c_e, 1.0, kC);
  sysC.setSteps(1);
  sysC.setDt(1.0);
  sysC.setT(0u, 0u, 0u, T1);

  for(unsigned int i = 0; i <= 1000; ++i) {
    sysA.solve();
    sysB.solve();
    sysC.solve();
  }
  
  sysA.saveTemperature("Out/sysA", 0);
  sysB.saveTemperature("Out/sysB", 0);
  sysC.saveTemperature("Out/sysC", 0);
  
  // try to find the distance traveled by energy
  for(unsigned int i = 0; i < n_x; ++i) {
    if(sysA.getT(i, 0, 0) <= T0) {
      std::cout << "sysA: " << i*d_x << std::endl; 
      break;
    }
  }
  
  for(unsigned int i = 0; i < n_x; ++i) {
    if(sysB.getT(i, 0, 0) <= T0) {
      std::cout << "sysB: " << i*d_x << std::endl; 
      break;
    }
  }
  
  for(unsigned int i = 0; i < n_x; ++i) {
    if(sysC.getT(i, 0, 0) <= T0) {
      std::cout << "sysC: " << i*d_x << std::endl; 
      break;
    }
  }
  
  // do the MPI_Finalise
  MPI_Finalize();
  
  return 0;
}

