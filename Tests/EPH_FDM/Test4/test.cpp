
#include <iostream>
#include <cmath>
#include <mpi.h>

#include "eph_fdm.h"

/* 
 * This example will create an electronic system and heat it from one end with
 * a laser source.
 */

// grid size
constexpr unsigned int n_x {1000};
constexpr unsigned int n_y {1};
constexpr unsigned int n_z {1};

// grid extent
constexpr double x_0 {-10.0};
constexpr double x_1 { 10.0};
constexpr double y_0 {-1.0};
constexpr double y_1 { 1.0};
constexpr double z_0 {-1.0};
constexpr double z_1 { 1.0};
constexpr double d_x {(x_1-x_0)/n_x};
constexpr double d_y {(y_1-y_0)/n_y};
constexpr double d_z {(z_1-z_0)/n_z};

// electronic system properties
constexpr double c_e {1.0}; // in eV
constexpr double rho_e {1.0}; // scaling factor
constexpr double kappa_e {1.0};
constexpr double T_e {1.0}; // initial temperature

constexpr double t_0 {0.0}; // laser turned on
constexpr double t_1 {1.0}; // laser turned off
constexpr double Q {10.0}; // laser power

// solver proeprties
constexpr unsigned int min_steps {1};
constexpr double dt {0.0001};

// simulation properties
constexpr double max_T {64.0};
constexpr double freq_t {1.0};

constexpr unsigned int max_steps {static_cast<unsigned int> (max_T/dt)};
constexpr unsigned int save_freq {static_cast<unsigned int> (freq_t/dt)};

constexpr const char* save {"Out/T_out"};

// energy transfer from laser to electronic system
double dE[n_x]; // per timestep in eV

// electrons as an FDM system
EPH_FDM electrons {n_x, n_y, n_z};

// laser profile
// v1 delta function with 
void laser_hit() {
  double d_V = d_x*d_y*d_z;
  for(int i = 0; i < n_x; ++i) {
    if((x_0 + i*d_x) == 0) dE[i] = Q / d_V;
    else dE[i] = 0.0;
  }
}

// laser profile to energy input into the system
void add_energy(double t) {
  if(t >= t_0 && t < (t_1)) {
    for(unsigned int i = 0; i < n_x; ++i) {
      electrons.setS(i, 0, 0, dE[i]);
    }
  }
  else {
    for(unsigned int i = 0; i < n_x; ++i) {
      electrons.setS(i, 0, 0, 0.0);
    }
  }
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
  
  /** define C_e as 1+sin(x) **/
  for(unsigned int i = 0; i < n_x; ++i) {
    electrons.set_C_e(i, 0, 0, 2.0 * c_e + c_e * sin(i*2.0*M_PI / n_x));
    electrons.set_kappa_e(i, 0, 0, 2.0 * kappa_e - kappa_e * sin(i*2.0*M_PI / n_x));
  }
  
  electrons.setFlag(n_x/4, 0, 0, 2);
  electrons.setFlag(3*n_x/4, 0, 0, 2);
  
  //electrons.saveState("T_before.data");
  
  for(unsigned int i = 0; i <= max_steps; ++i) {
    add_energy(i*dt);
    electrons.solve();
    
    if(i%save_freq == 0) {
      printf("Saving step %06d; t = %8.3f ps; T = %8.3f K;\n", i, dt*i, electrons.calcTtotal());
      switch (i/save_freq) {
        case 0:
        case 1:
        case 2:
        case 4:
        case 8:
        case 16:
        case 32:
        case 64:
          electrons.saveTemperature(save, i/save_freq);
          break;
        default:
          break;
      }
    }
  }
  
  //electrons.saveState("T_after.data");
  
  // do the MPI_Finalise
  MPI_Finalize();
  
  return 0;
}

