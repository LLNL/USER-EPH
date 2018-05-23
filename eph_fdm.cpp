/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

// external headers
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>

#include <mpi.h>

// change this to iostream
#include <cstdio>

// internal headers
#include "eph_fdm.h"

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
  
  ntotal = nx*ny*nz;
  resize_vectors();
  
  //init();
  
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

EPH_FDM::EPH_FDM(const unsigned int nx, const unsigned int ny, const unsigned int nz) {
  this->nx = nx;
  this->ny = ny;
  this->nz = nz;
  
  ntotal = nx*ny*nz;
  resize_vectors();
  
  setDt(1.0);
  setSteps(1);
  setBox(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
}

void EPH_FDM::resize_vectors() {
  T_e.resize(ntotal, 0.0);
  dT_e.resize(ntotal, 0.0);
  dT_e_x.resize(ntotal, 0.0);
  dT_e_y.resize(ntotal, 0.0);
  dT_e_z.resize(ntotal, 0.0);
  ddT_e.resize(ntotal, 0.0);
  
  C_e.resize(ntotal, 0.0);
  rho_e.resize(ntotal, 0.0);
  kappa_e.resize(ntotal, 0.0);
  
  S_e.resize(ntotal, 0.0);
  flag.resize(ntotal, 1);  
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
    
    //if(r > 0.5) {
    if(r > 0.4) {
      inner_dt = 0.4 * inner_dt / r; // be more conservative
      //inner_dt = 0.5 * inner_dt / r; // get new timestep
      //new_steps = ((unsigned int)(dt / inner_dt)) + 1;
      new_steps = std::floor(dt / inner_dt) + 1;
      inner_dt = dt / new_steps;
    }
    
    for(int n = 0; n < new_steps; n++) {
      // calculate derivatives
      #ifdef EPH_TESTING
      #pragma omp parallel 
      #endif
      {
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
      }
      
      // calculate second derivative
      #ifdef EPH_TESTING
      #pragma omp parallel 
      #endif
      {
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
      }
      
      /* TODO: there might be an issue with grid volume here */
      // do the actual step
      #ifdef EPH_TESTING
      #pragma omp parallel 
      #endif
      {
      for(int i = 0; i < ntotal; i++) {
        double prescaler = rho_e[i] * C_e[i];
        
        if(prescaler > 0.0 && flag[i] > 0)
          T_e[i] += (ddT_e[i] + dT_e[i] + S_e[i]) / prescaler * inner_dt;
        
        // energy conservation issues
        /* Add a sanity check somewhere for this */
        if(T_e[i] < 0.0)
          T_e[i] = 0.0;
      }
      }
    
      // zero arrays
      std::fill(dT_e_x.begin(), dT_e_x.end(), 0.0); 
      std::fill(dT_e_y.begin(), dT_e_y.end(), 0.0); 
      std::fill(dT_e_z.begin(), dT_e_z.end(), 0.0); 
      std::fill(ddT_e.begin(), ddT_e.end(), 0.0); 
    }
    
    // zero energy exchange array
    //std::fill(dT_e.begin(), dT_e.end(), 0.0);
  }
  
  syncAfter();
}

void EPH_FDM::saveState(const char* file) {
  FILE *fd = fopen(file, "w");
  
  if(!fd) throw std::runtime_error("eph_fdm: could not open supplied file"); 
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
  
  if(!fd) throw std::runtime_error("eph_fdm: could not open supplied file"); // give an error?
  
  // this is needed for visit Point3D
  fprintf(fd, "x y z Te\n");
  
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
    MPI_Allreduce(MPI_IN_PLACE, dT_e.data(), ntotal, MPI_DOUBLE, MPI_SUM, world);
}

void EPH_FDM::syncAfter() {
  // zero arrays
  // TODO: substitute with std::fill
  std::fill(dT_e_x.begin(), dT_e_x.end(), 0.0); 
  std::fill(dT_e_y.begin(), dT_e_y.end(), 0.0); 
  std::fill(dT_e_z.begin(), dT_e_z.end(), 0.0); 
  std::fill(ddT_e.begin(), ddT_e.end(), 0.0);
  std::fill(dT_e.begin(), dT_e.end(), 0.0);
  
  // synchronize electronic temperature
  if(nrPS > 0)
    MPI_Bcast(T_e.data(), ntotal, MPI_DOUBLE, 0, world);
}
