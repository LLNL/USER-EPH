
#include <iostream>
#include <mpi.h>

#include "eph_fdm.h"

/* 
 * This example will create an electronic system and heat it from one end with
 * a laser source.
 */

// grid size
constexpr unsigned int n_x {10000};
constexpr unsigned int n_y {1};
constexpr unsigned int n_z {1};

// grid extent
constexpr double x_0 {-5000.0};
constexpr double x_1 {5000.0};
constexpr double y_0 {0.0};
constexpr double y_1 {28.16};
constexpr double z_0 {0.0};
constexpr double z_1 {28.16};
constexpr double d_x {(x_1-x_0)/n_x};
constexpr double d_y {(y_1-y_0)/n_y};
constexpr double d_z {(z_1-z_0)/n_z};

// electronic system properties
constexpr double c_e {3.5e-6}; // in eV
constexpr double rho_e {1.0}; // scaling factor
constexpr double kappa_e {0.1248};
constexpr double T_e {300.0}; // initial temperature

// solver proeprties
constexpr unsigned int min_steps {1};
constexpr double dt {0.0001};

// simulation properties
constexpr unsigned int max_steps {10000};
constexpr unsigned int save_freq {1000};

constexpr const char* save {"Out/T_out"};

// energy transfer from laser to electronic system
double dE[n_x]; // per timestep in eV

// electrons as an FDM system
EPH_FDM electrons {n_x, n_y, n_z};

// laser profile
// v1 delta function with 
void laser_hit() {
  
}

// laser profile to energy input into the system
void add_energy() {
  
}

int main(int args, char **argv) {
  // do the MPI_Init
  int my_id;
  int nr_ps;
  
  MPI_Init(&args, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &nr_ps);
  
  electrons.setComm(MPI_COMM_WORLD, my_id, nr_ps);
  electrons.setBox(x_0, x_1, y_0, y_1, z_0, z_1);
  
  electrons.setConstants(T_e, c_e, rho_e, kappa_e);
  electrons.setSteps(min_steps); // minimum nr. of steps
  electrons.setDt(dt); // in ps
  
  laser_hit();
  
  for(unsigned int i = 0; i <= max_steps; ++i) {
    add_energy();
    electrons.solve();
    
    if(i%save_freq == 0) {
      printf("Saving step %06d; t = %8.3f ps; T = %8.3f K;\n", i, dt*i, electrons.calcTtotal());
      electrons.saveTemperature(save, i/save_freq);
    }
  }
  // do the MPI_Finalise
  MPI_Finalize();
  
  return 0;
}

