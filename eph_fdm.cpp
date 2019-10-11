/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */
#if 0
// external headers
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>

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
    /* find smallest C_e and rho_e and largest kappa */
    double c_min = C_e[0];
    double rho_min = rho_e[0];
    double kappa_max = kappa_e[0];
    
    for(unsigned int i = 1; i < ntotal; ++i) {
      if(flag[i] > 0) {
        if(C_e[i] < c_min) c_min = C_e[i];
        if(rho_e[i] < rho_min) rho_min = rho_e[i];
        if(kappa_e[i] > kappa_max) kappa_max = kappa_e[i];
      }
    }
    
    double r = dtdxdydz / c_min / rho_min * kappa_max;
    
    unsigned int new_steps = steps;
    
    // This will become unstable if there are any large fluctuations 
    // during the solving process; calling this at every step is expensive
    //if(r > 0.5) {
    if(r > 0.4) {
      inner_dt = 0.4 * inner_dt / r; // be more conservative
      //inner_dt = 0.5 * inner_dt / r; // get new timestep
      //new_steps = ((unsigned int)(dt / inner_dt)) + 1;
      //new_steps = std::floor(dt / inner_dt) + 1;
      new_steps = static_cast<unsigned int> ((dt / inner_dt)+1.0);
      inner_dt = dt / new_steps;
    }
    
    for(int n = 0; n < new_steps; n++) {
      std::fill(ddT_e.begin(), ddT_e.end(), 0.0);
      
      #ifdef EPH_OMP
      #pragma omp parallel for collapse(3)
      #endif
      //for(unsigned int r = 0 ; r < ntotal; ++r) {
      for(unsigned int k = 0; k < nz; ++k) {
        for(unsigned int j = 0; j < ny; ++j) {
          for(unsigned int i = 0; i < nx; ++i) {
            //unsigned int i, j, k;
            //i = r%nx;
            //j = (r - i)%ny;
            //k = (r - i - j)%nz;
            
            unsigned int q, p;
            unsigned int r;
            r = i + j*nx + k*nx*ny;
            
            ddT_e[r] = 0.0;
            if(flag[r] == 2) continue;
            
            // +- dx
            if(i > 0) p = (i-1) + j*nx + k*nx*ny;
            else p = (nx-1) + j*nx + k*nx*ny;
            
            if(i < (nx - 1)) q = (i+1) + j*nx + k*nx*ny;
            else q = j*nx + k*nx*ny;
            
            if(flag[q] == 2) q = r;
            else if(flag[p] == 2) p = r;
            
            ddT_e[r] += (kappa_e[q]-kappa_e[p]) * (T_e[q] - T_e[p]) / dx / dx / 4.0;
            ddT_e[r] += kappa_e[r] * ((T_e[q]+T_e[p]-2.0*T_e[r]) / dx / dx);
            
            // +- dy
            if(j > 0) p = i + (j-1)*nx + k*nx*ny;
            else p = i + (ny-1)*nx + k*nx*ny;
            
            if(j < (ny - 1)) q = i + (j+1)*nx + k*nx*ny;
            else q = i + k*nx*ny;
            
            if(flag[q] == 2) q = r;
            else if(flag[p] == 2) p = r;
            
            ddT_e[r] += (kappa_e[q]-kappa_e[p]) * (T_e[q] - T_e[p]) / dy / dy / 4.0;
            ddT_e[r] += kappa_e[r] * ((T_e[q]+T_e[p]-2.0*T_e[r]) / dy / dy);
            
            // +- dz
            if(k > 0) p = i + j*nx + (k-1)*nx*ny;
            else p = i + j*nx + (nz-1)*nx*ny;
            
            if(k < (nz - 1)) q = i + j*nx + (k+1)*nx*ny;
            else q = i + j*nx;
            
            if(flag[q] == 2) q = r;
            else if(flag[p] == 2) p = r;
            
            ddT_e[r] += (kappa_e[q]-kappa_e[p]) * (T_e[q] - T_e[p]) / dz / dz / 4.0;
            ddT_e[r] += kappa_e[r] * ((T_e[q]+T_e[p]-2.0*T_e[r]) / dz / dz);
          }
        }
      }
      
      /* TODO: there might be an issue with grid volume here */
      // do the actual step
      #ifdef EPH_OMP
      #pragma omp parallel for
      #endif
      for(int i = 0; i < ntotal; i++) {
        double prescaler = rho_e[i] * C_e[i];
        
        if(prescaler > 0.0) {
          switch(flag[i]) {
            case 1:
              T_e[i] += (ddT_e[i] + dT_e[i] + S_e[i]) / prescaler * inner_dt;
              break;
            case 2:
              T_e[i] = 0.0;
              break;
          }
        }
        else {
          throw std::runtime_error("eph_fdm.cpp solve(): negative prescaler");
        }
        
        // energy conservation issues
        /* Add a sanity check somewhere for this */
        if(T_e[i] < 0.0)
          T_e[i] = 0.0;
      }
    }
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
  std::fill(dT_e.begin(), dT_e.end(), 0.0);
  
  // synchronize electronic temperature
  if(nrPS > 0)
    MPI_Bcast(T_e.data(), ntotal, MPI_DOUBLE, 0, world);
}
#endif
