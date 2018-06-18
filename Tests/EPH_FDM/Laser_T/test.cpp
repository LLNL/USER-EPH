
#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>

#include "eph_fdm.h"
#include "eph_spline.h"

/* 
 * This example will create an electronic system and heat it from one end with
 * a laser source.
 */

// grid size
constexpr unsigned int n_x {10000};
constexpr unsigned int n_y {1};
constexpr unsigned int n_z {1};
constexpr unsigned int n_total {n_x*n_y*n_z};

// grid extent
constexpr double x_0 {-1000.0};
constexpr double x_1 { 1000000.0};
constexpr double y_0 { 0.0};
constexpr double y_1 { 28.16};
constexpr double z_0 { 0.0};
constexpr double z_1 { 28.16};
constexpr double d_x {(x_1-x_0)/n_x};
constexpr double d_y {(y_1-y_0)/n_y};
constexpr double d_z {(z_1-z_0)/n_z};
constexpr double d_V {d_x*d_y*d_z};

// electronic system properties
constexpr double c_e {3.5e-6}; // in eV
constexpr double rho_e {1.0}; // scaling factor
constexpr double kappa_e {5.68e-2};
constexpr double T_e {300.0}; // initial temperature

constexpr double c_l {2.3e-5}; // in eV
constexpr double rho_l {1.0}; // scaling factor
constexpr double kappa_l {6.2e-3};
constexpr double T_l {300.0}; // initial temperature

constexpr double g_eph {6.6e-6}; // electron-phonon coupling parameter

// this are special electronic conductivity numbers from Zhigilei
constexpr double vF2 {2.25e+8}; // (A/ps)^2
constexpr double A {1.4e-6}; // K^2/ps
constexpr double B {16.24}; // K/ps

constexpr double t_0 {0.0}; // laser turned on
constexpr double t_1 {15000.0}; // laser turned off
//constexpr double Q {62.415}; // laser power eV/ps/Ang^2
//constexpr double Q {5.0e+4}; // laser power in eV/ps
constexpr double Q {5.0e+4}; // laser power in eV/ps

// solver proeprties
constexpr unsigned int min_steps {10};
constexpr double dt {0.1};

// simulation properties
constexpr double max_T {5000000.0};
constexpr double freq_t {1000.0};

constexpr unsigned int max_steps {static_cast<unsigned int> (max_T/dt)};
constexpr unsigned int save_freq {static_cast<unsigned int> (freq_t/dt)};

constexpr const char* save_e {"Out/T_e_out"};
constexpr const char* save_l {"Out/T_l_out"};

// energy transfer from laser to electronic system
double dE[n_x]; // per timestep in eV
int on_off {0};

// temperature dependent parameters
EPH_Spline c_e_T;
EPH_Spline g_eph_T;

// electrons as an FDM system
EPH_FDM electrons {n_x, n_y, n_z};
EPH_FDM lattice {n_x, n_y, n_z};

// laser profile
// v1 delta function with at x_0 = 0
void laser_hit() {
  for(int i = 0; i < n_x; ++i) {
    if((x_0 + i*d_x) >= 0.0 && (x_0 + i*d_x) < d_x) dE[i] = Q / d_V;
    else dE[i] = 0.0;
  }
}

// laser profile to energy input into the system
void add_energy(double t) {
  if((t >= t_0 && t < (t_1)) && !on_off) {
    for(unsigned int i = 0; i < n_x; ++i) {
      electrons.setS(i, 0, 0, dE[i]);
    }
    on_off = 1;
  }
  else if((t < t_0 || t >= t_1) && on_off) {
    for(unsigned int i = 0; i < n_x; ++i) {
      electrons.setS(i, 0, 0, 0.0);
    }
    on_off = 0;
  }
}

void electron_phonon() {
  for(unsigned int k = 0; k < n_z; ++k) {
    for(unsigned int j = 0; j < n_y; ++j) {
      for(unsigned int i = 0; i < n_x; ++i) {
        double varTe = electrons.getT(i,j,k);
        double dT = varTe - lattice.getT(i, j, k);
        double varG = g_eph_T.GetValue(varTe);
        double dE = dT * varG * d_V * dt;
          
        electrons.insertEnergy(i, j, k, -dE);
        lattice.insertEnergy(i, j, k, dE);
      }
    }
  }
}

void change_parameters() {
  for(unsigned int k = 0; k < n_z; ++k) {
    for(unsigned int j = 0; j < n_y; ++j) {
      for(unsigned int i = 0; i < n_x; ++i) {
        double varTe = electrons.getT(i, j, k);
        double varTl = lattice.getT(i, j, k);
        double varC = c_e_T.GetValue(varTe);
        
        if((varTe > 0.0) && (varTl > 0.0)) {
          electrons.set_C_e(i,j,k, varC);
        
          double varK = vF2 * varC * 1.0/(A*varTe*varTe+B*varTl);
          electrons.set_kappa_e(i,j,k, varK);
          
          //std::cout << varTe << " " << varTl << " " << varC << " " << varK << std::endl;
        }
      }
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
  
  // initialise splines
  {
  std::ifstream fd("Data/C_e.data");
  
  double valT, valC;
  c_e_T.SetDiscretisation(0.0, 100.0);
  
  while(! (fd.eof())) {
    fd >> valT >> valC;
    c_e_T << valC;
  }
  
  c_e_T << true;
  }
  
  {
  std::ifstream fd("Data/C_e.data");
  
  double valT, valG;
  g_eph_T.SetDiscretisation(0.0, 100.0);
  
  while(! (fd.eof())) {
    fd >> valT >> valG;
    g_eph_T << valG;
  }
  
  g_eph_T << true;
  }
  
  electrons.setComm(MPI_COMM_WORLD, my_id, nr_ps);
  electrons.setBox(x_0, x_1, y_0, y_1, z_0, z_1);
  
  electrons.setConstants(T_e, c_e, rho_e, kappa_e);
  electrons.setSteps(min_steps); // minimum nr. of steps
  electrons.setDt(dt); // in ps
  
  lattice.setComm(MPI_COMM_WORLD, my_id, nr_ps);
  lattice.setBox(x_0, x_1, y_0, y_1, z_0, z_1);
  lattice.setConstants(T_l, c_l, rho_l, kappa_l);
  lattice.setSteps(min_steps); // minimum nr. of steps
  lattice.setDt(dt); // in ps
  
  laser_hit();
  
  /** define C_e as 1+sin(x) **/
  /*
  for(unsigned int i = 0; i < n_x; ++i) {
    electrons.set_C_e(i, 0, 0, 2.0 * c_e + c_e * sin(i*2.0*M_PI / n_x));
    electrons.set_kappa_e(i, 0, 0, 2.0 * kappa_e - kappa_e * sin(i*2.0*M_PI / n_x));
  }*/
  
  electrons.setFlag(0u, 0u, 0u, 2);
  lattice.setFlag(0u, 0u, 0u, 2);
  
  //electrons.saveState("T_before.data");
  
  for(unsigned int i = 0; i <= max_steps; ++i) {
    change_parameters();
    add_energy(i*dt);
    electron_phonon();
    electrons.solve();
    lattice.solve();
    
    //std::cout << c_e << " " << kappa_e << std::endl;
    //throw;
    
    if(i%save_freq == 0) {
      printf("Saving step %06d in %06d; t = %10.3f ps; T_e = %8.3f K; T_l = %8.3f K;\n", i, i/save_freq, dt*i, electrons.calcTtotal(), lattice.calcTtotal());
      switch (i/save_freq) {
        case 0:
        case 1:
        case 10:
        case 100:
        case 500:
        case 1000:
        case 1500:
        case 2000:
        case 2500:
        case 3000:
        case 3500:
        case 4000:
        case 4500:
        case 5000:
          electrons.saveTemperature(save_e, i/save_freq);
          lattice.saveTemperature(save_l, i/save_freq);
          break;
        default:
          electrons.saveTemperature(save_e, 999999);
          lattice.saveTemperature(save_l, 999999);
          break;
      }
    }
  }
  
  //electrons.saveState("T_after.data");
  
  // do the MPI_Finalise
  MPI_Finalize();
  
  return 0;
}

