/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: arturt@ut.ee
 */

#ifdef DEBUG_EPH
#include <iostream>
#endif

// external headers
#include <fstream>
#include <string>
#include <cmath>

// internal headers
#include "eph_fdm.h"


EPH_FDM::EPH_FDM(const char* file) {
  std::ifstream fd;
  fd.open(file);
  
  if(!fd.is_open()) {
    #ifdef DEBUG_EPH
    std::cout << "Opening file failed: " << file << std::endl;
    #endif
    return;
  }
  
  char line[lineLength];
  
  // get number of grid points
  fd >> nx;
  fd >> ny;
  fd >> nz;
  
  init();
  
  for(int i = 0; i < ntotal; i++) {
    
  }
}

EPH_FDM::EPH_FDM(unsigned int nx, unsigned int ny, unsigned int nz) {
  this->nx = nx;
  this->ny = ny;
  this->nz = nz;
  
  init();
}

void EPH_FDM::init() {
  ntotal = nx*ny*nz;
  
  x0 = y0 = z0 = 0.0;
  x1 = y1 = z1 = 1.0;
  
  dx = (x1 - x0) / nx;
  dy = (y1 - y0) / ny;
  dz = (z1 - z0) / nz;
  
  dV = dx * dy * dz;
  
  T_e = new double[ntotal];
  dT_e = new double[ntotal];
  dT_e_x = new double[ntotal];
  dT_e_y = new double[ntotal];
  dT_e_z = new double[ntotal];
  ddT_e = new double[ntotal];
  S_e = new double[ntotal];
  
  for(int i = 0; i < ntotal; i++) {
    T_e[i] = 0.0;
    dT_e[i] = 0.0;
    dT_e_x[i] = 0.0;
    dT_e_y[i] = 0.0;
    dT_e_z[i] = 0.0;
    ddT_e[i] = 0.0;
    S_e[i] = 0.0;
  }
  
  steps = 1;
  dt = 1.0;
}

EPH_FDM::~EPH_FDM() {
  delete T_e;
  delete dT_e;
  delete dT_e_x;
  delete dT_e_y;
  delete dT_e_z;
  delete ddT_e;
  
  delete S_e;
}

double EPH_FDM::calcEtotal() { // calculate the total energy of the electronic system
  double result = 0.0;
  
  for(int i = 0; i < ntotal; i++) {
    result += rho_e * C_e * T_e[i];
  }
  
  return result;
}

unsigned int EPH_FDM::solve() {
  syncBefore();
  
  // check for stability
  double inner_dt = dt / steps;
  
  double dtdxdydz = inner_dt * (1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz);
  double r = dtdxdydz / C_e / rho_e * kappa_e;
  
  unsigned int new_steps = steps;
  
  if(r > 0.5) {
    inner_dt = 0.5 * inner_dt / r;
    new_steps = ((unsigned int)(dt / inner_dt)) + 1;
    inner_dt = dt / new_steps;
  }
  
  for(int n = 0; n < new_steps; n++) {
    // calculate derivatives
    for(int k = 0; k < nz; k++) {
      for(int j = 0; j < ny; j++) {
        for(int i = 0; i < nx; i++) {
          unsigned int p, q;
          p = i + j*nx + k*nx*ny;
          q = (i+1)%nx + j*nx + k*nx*ny;
          dT_e_x[p] = (T_e[q] - T_e[p])/dx;
          
          q = i + ((j+1)%ny)*nx + k*nx*ny;
          dT_e_y[p] = (T_e[q] - T_e[p])/dy;
          
          q = i + j*nx + ((k+1)%nz)*nx*ny;
          dT_e_z[p] = (T_e[q] - T_e[p])/dz;
        }
      }
    }
    
    // calculate second derivative
    for(int k = 0; k < nz; k++) {
      for(int j = 0; j < ny; j++) {
        for(int i = 0; i < nx; i++) {
          unsigned int p, q;
          p = i + j*nx + k*nx*nz;
          q = (i+1)%nx + j*nx + k*nx*ny;
          ddT_e[q] += kappa_e * (dT_e_x[q] - dT_e_x[p])/dx;
          
          q = i + ((j+1)%ny)*nx + k*nx*ny;
          ddT_e[q] += kappa_e * (dT_e_y[q] - dT_e_y[p])/dy;
          
          q = i + j*nx + ((k+1)%nz)*nx*ny;
          ddT_e[q] += kappa_e * (dT_e_z[q] - dT_e_z[p])/dz;
        }
      }
    }
    
    // do the actual step
    double prescaler = rho_e * C_e;
    
    for(int i = 0; i < ntotal; i++) {
      //std::cout << i << " " << T_e[i] << " " << ddT_e[i] / prescaler * inner_dt << std::endl;
      if(prescaler > 0.0)
        T_e[i] += ddT_e[i] / prescaler * inner_dt + dT_e[i] * inner_dt + S_e[i] * inner_dt;
      
      // energy conservation issues
      if(T_e[i] < 0.0)
        T_e[i] = 0.0;
    }
  
    // zero arrays
    for(int i = 0; i < ntotal; i++) {
      dT_e_x[i] = 0.0;
      dT_e_y[i] = 0.0;
      dT_e_z[i] = 0.0;
      ddT_e[i] = 0.0;
    }
  }
  
  for(int i = 0; i < ntotal; i++) {
    dT_e[i] = 0.0;
  }
  
  //std::cout << "STEPS: " << new_steps << std::endl;
  
  syncAfter();
  
  return new_steps;
}

void EPH_FDM::saveState(const char* file) {
  
}

void EPH_FDM::syncBefore() {
  
}

void EPH_FDM::syncAfter() {
  
}
