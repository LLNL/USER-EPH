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

#include <mpi.h>

// change this to iostream
#include <cstdio>

// internal headers
#include "eph_fdm.h"

// TODO: change to iostream, have been doing too much C coding
EPH_FDM::EPH_FDM(const char* file) {
  std::ifstream fd;
  
  fd.open(file);
  
  if(!fd.is_open()) {
    return; // break the code
  }
  
  // 3 first lines are comments
  char line[lineLength];
  fd.getline(line, lineLength);
  fd.getline(line, lineLength);
  fd.getline(line, lineLength);
  
  // next line defines grid size
  fd >> nx;
  fd >> ny;
  fd >> nz;
  
  init();
  
  fd >> steps;
  
  // define box size
  fd >> x0; fd >> x1;
  fd >> y0; fd >> y1;
  fd >> z0; fd >> z1;
  
  setBox(x0, x1, y0, y1, z0, z1);
  
  // read grid values
  for(int i = 0; i < ntotal; ++i) {
    int lx, ly, lz;
    fd >> lx; fd >> ly; fd >> lz;
    unsigned int index = lx + ly * nx + lz * nx * ny;   
    fd >> T_e[index];
    fd >> S_e[index];
    fd >> rho_e[index];
    fd >> C_e[index];
    fd >> kappa_e[index];
    fd >> flag[index];
  }
  
  fd.close();
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
  
  rho_e = new double[ntotal];
  C_e = new double[ntotal];
  kappa_e = new double[ntotal];
  flag = new int[ntotal];
  
  for(int i = 0; i < ntotal; i++) {
    T_e[i] = 0.0;
    dT_e[i] = 0.0;
    dT_e_x[i] = 0.0;
    dT_e_y[i] = 0.0;
    dT_e_z[i] = 0.0;
    ddT_e[i] = 0.0;
    S_e[i] = 0.0;
    
    rho_e[i] = 0.0;
    C_e[i] = 0.0;
    kappa_e[i] = 0.0;
    flag[i] = -1;
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
  
  delete C_e;
  delete rho_e;
  delete kappa_e;
  delete flag;
}

void EPH_FDM::solve() {
  syncBefore();
  
  if(myID == 0) {  
    // this is strongly inspired by fix_ttm
    // check for stability
    double inner_dt = dt / steps;
    
    double dtdxdydz = inner_dt * (1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz);
    
    // this does not work corectly
    double r = dtdxdydz / C_e[0] / rho_e[0] * kappa_e[0];
    
    unsigned int new_steps = steps;
    
    if(r > 0.5) {
      inner_dt = 0.5 * inner_dt / r;
      new_steps = ((unsigned int)(dt / inner_dt)) + 1;
      inner_dt = dt / new_steps;
    }
    
    for(int n = 0; n < new_steps; n++) {
      // calculate derivatives
      // TODO: add OMP parallelism here
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
            ddT_e[q] += (kappa_e[q] * dT_e_x[q] - kappa_e[p] * dT_e_x[p])/dx;
            
            q = i + ((j+1)%ny)*nx + k*nx*ny;
            ddT_e[q] += (kappa_e[q] * dT_e_y[q] - kappa_e[p] * dT_e_y[p])/dy;
            
            q = i + j*nx + ((k+1)%nz)*nx*ny;
            ddT_e[q] += (kappa_e[q] * dT_e_z[q] - kappa_e[p] * dT_e_z[p])/dz;
          }
        }
      }
      
      /* TODO: there might be an issue with grid volume here */
      // do the actual step
      for(int i = 0; i < ntotal; i++) {
        //std::cout << i << " " << T_e[i] << " " << ddT_e[i] / prescaler * inner_dt << std::endl;
        double prescaler = rho_e[i] * C_e[i];
        
        if(prescaler > 0.0 && flag[i] > 0)
          T_e[i] += (ddT_e[i] + dT_e[i] + S_e[i]) / prescaler * inner_dt;
        
        // energy conservation issues
        /* Add a sanity check somewhere for this
        if(T_e[i] < 0.0)
          T_e[i] = 0.0;
        */
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
  }
  
  syncAfter();
}

void EPH_FDM::saveState(const char* file) {
  FILE *fd = fopen(file, "w");
  
  if(!fd) return; // give an error?
  // 3 first lines are comments
  fprintf(fd, "# A comment\n");
  fprintf(fd, "#\n");
  fprintf(fd, "#\n");
  // next line is grid size and min number of steps
  fprintf(fd, "%d %d %d %d\n", nx, ny, nz, steps);
  // next we have box size
  fprintf(fd, "%.6e %.6e\n", x0, x1);
  fprintf(fd, "%.6e %.6e\n", y0, y1);
  fprintf(fd, "%.6e %.6e\n", z0, z1);
  
  // finally we have grid values
  for(int k = 0; k < nz; ++k) {
    for(int j = 0; j < ny; ++j) {
      for(int i = 0; i < nx; ++i) {
        unsigned int index = i + j * nx + k * nx * ny;        
        fprintf(fd, "%d %d %d %.6e %.6e %.6e %.6e %.6e %d\n", 
          i, j, k, T_e[index], S_e[index], rho_e[index], C_e[index], kappa_e[index], flag[index]);
      }
    }
  }
  
  fclose(fd);
}

void EPH_FDM::saveTemperature(const char* file, int n) {
  char fn[512];
  sprintf(fn, "%s_%06d", file, n);
  FILE *fd = fopen(fn, "w");
  
  if(!fd) return; // give an error?
  
  for(int k = 0; k < nz; ++k) {
    for(int j = 0; j < ny; ++j) {
      for(int i = 0; i < nx; ++i) {
        unsigned int index = i + j * nx + k * nx * ny;
        
        double x = x0 + i * dx;
        double y = y0 + j * dy;
        double z = z0 + k * dz;
        
        fprintf(fd, "%.6e %.6e %.6e %.6e\n", x, y, z, T_e[index]);
      }
    }
  }
  
  fclose(fd);
}

void EPH_FDM::syncBefore() {
  // sync energy transfer
  if(nrPS > 0)
    MPI_Allreduce(MPI_IN_PLACE, &(dT_e[0]), ntotal, MPI_DOUBLE, MPI_SUM, world);
}

void EPH_FDM::syncAfter() {
  // zero arrays
  for(int i = 0; i < ntotal; i++) {
    dT_e_x[i] = 0.0;
    dT_e_y[i] = 0.0;
    dT_e_z[i] = 0.0;
    ddT_e[i] = 0.0;
    dT_e[i] = 0.0;
  }
  
  // synchronize electronic temperature
  if(nrPS > 0)
    MPI_Bcast(&(T_e[0]), ntotal, MPI_DOUBLE, 0, world);
}
