/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

#ifndef EPH_FDM_H
#define EPH_FDM_H

#include "eph_spline.h"

#include <iostream>
#include <cassert>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <cstring>
#include <numeric>

#include <mpi.h>

class EPH_FDM
{
  public:
    // default constructor to create a grid with one point
    EPH_FDM() = default;
    
    EPH_FDM(
      size_t in_nx, size_t in_ny, size_t in_nz,
      double in_x0, double in_x1, 
      double in_y0, double in_y1,
      double in_z0, double in_z1,
      double in_T_e, double in_C_e, double in_rho_e, double in_kappa_e) : 
      nx {in_nx},
      ny {in_ny},
      nz {in_nz}
    {
      // TODO: fix me
      resize_vectors(nx, ny, nz);
      
      set_box_dimensions(in_x0, in_x1, in_y0, in_y1, in_z0, in_z1);
      set_constants(in_T_e, in_C_e, in_rho_e, in_kappa_e);
      set_steps(1);
      set_dt(1);
      parameter_filename = "NULL";
    }
    
    EPH_FDM(const char *in_filename) 
    {
      std::ifstream fd {in_filename}; assert(fd); // break code here
      
      // 3 first lines are comments
      char line[lineLength];
      fd.getline(line, lineLength); assert(line[0] == '#');
      fd.getline(line, lineLength); assert(line[0] == '#');
      fd.getline(line, lineLength); assert(line[0] == '#');
      
      // next line defines grid size
      fd >> nx >> ny >> nz;
      resize_vectors(nx, ny, nz);
      
      fd >> steps;
      
      // define box size
      fd >> x0 >> x1
         >> y0 >> y1
         >> z0 >> z1;
      
      set_box_dimensions(x0, x1, y0, y1, z0, z1);
      
      fd >> parameter_filename;
      
      // load temperature dependent parameters 
      if(parameter_filename != "NULL")
      {
        std::ifstream fd {parameter_filename}; assert(fd);
        
        // 3 first lines are comments
        char line[lineLength];
        fd.getline(line, lineLength); assert(line[0] == '#');
        fd.getline(line, lineLength); assert(line[0] == '#');
        fd.getline(line, lineLength); assert(line[0] == '#');
        
        size_t n; fd >> n;
        double dT; fd >> dT;
        
        std::vector<double> in_C_e_T(n);
        std::vector<double> in_kappa_e_T(n);
        
        for(size_t i = 0; i < n; ++i)
        {
          fd >> in_C_e_T[i] >> in_kappa_e_T[i]; // read data from file
        }
        
        C_e_T = Spline(dT, in_C_e_T);
        kappa_e_T = Spline(dT, in_kappa_e_T);
      }
      
      // read grid values
      for(size_t i = 0; i != ntotal; ++i) {
        int lx, ly, lz;
        fd >> lx >> ly >> lz;
        size_t index = lx + ly * nx + lz * nx * ny;   
        fd >> T_e[index]
           >> S_e[index]
           >> rho_e[index]
           >> C_e[index]
           >> kappa_e[index]
           >> flag[index]
           >> T_dynamic_flag[index];
      }
    }
    
    // set box dimensions in Ang
    void set_box_dimensions(
      double in_x0, double in_x1, 
      double in_y0, double in_y1,
      double in_z0, double in_z1)
    {
      x0 = in_x0; x1 = in_x1;
      y0 = in_y0; y1 = in_y1;
      z0 = in_z0; z1 = in_z1;
      
      assert(x0 < x1);
      assert(y0 < y1);
      assert(z0 < z1);
      
      dx = (in_x1 - in_x0)/nx;
      dy = (in_y1 - in_y0)/ny;
      dz = (in_z1 - in_z0)/nz;
      
      dV = dx*dy*dz;
    }
    
    // set constant values for as grid parameters
    void set_constants(
      double in_T_e, double in_C_e, double in_rho_e, double in_kappa_e)
    {
      std::fill(T_e.begin(), T_e.end(), in_T_e);
      std::fill(rho_e.begin(), rho_e.end(), in_rho_e);
      std::fill(C_e.begin(), C_e.end(), in_C_e);
      std::fill(kappa_e.begin(), kappa_e.end(), in_kappa_e);
      std::fill(flag.begin(), flag.end(), 1);
      std::fill(T_dynamic_flag.begin(), T_dynamic_flag.end(), false);
    }
    
    void set_dt(double in_dt) 
    {
      dt = in_dt;
    }
    
    void set_steps(size_t in_steps)
    {
      steps = in_steps;
    }
    
    void set_comm(MPI_Comm in_comm, int in_myID, int in_nrPS) 
    {
      world = in_comm;
      myID = in_myID;
      nrPS = in_nrPS;
    } 
    
    // add energy into a cell
    void insert_energy(double x, double y, double z, double E) 
    {
      unsigned int index = get_index(x, y, z);
      double prescale = dV * dt;
      
      // convert energy into power per area
      dT_e[index] += E / prescale;
    }
    
    // get temperature of a cell
    double get_T(double x, double y, double z) const 
    {
      unsigned int index = get_index(x, y, z);
      
      return T_e[index];
    }
    
    double get_T_total() const 
    {
      double result {std::accumulate(T_e.begin(), T_e.end(), 0.)};
      
      result /= ntotal; // this calculates the average temperature
  
      return result;
    }
    
    void save_temperature(const char* in_filename, int in_n) const 
    {
      char fn[512];
      sprintf(fn, "%s_%06d", in_filename, in_n);
      FILE *fd = fopen(fn, "w");
      
      assert(fd > 0);
      
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
    
    void save_state(const char* in_filename) const 
    {
      FILE *fd = fopen(in_filename, "w");
      
      assert(fd > 0);
      
      // 3 first lines are comments
      fprintf(fd, "# A comment\n");
      fprintf(fd, "#\n");
      fprintf(fd, "#\n");
      
      // next line is grid size and min number of steps
      fprintf(fd, "%ld %ld %ld %ld\n", nx, ny, nz, steps);
      
      // next we have box size
      fprintf(fd, "%.6e %.6e\n", x0, x1);
      fprintf(fd, "%.6e %.6e\n", y0, y1);
      fprintf(fd, "%.6e %.6e\n", z0, z1);
      
      // filename for temperature dependent parameters
      fprintf(fd, "%s\n", parameter_filename.c_str());
      
      // finally we have grid values
      for(int k = 0; k < nz; ++k) 
      {
        for(int j = 0; j < ny; ++j) 
        {
          for(int i = 0; i < nx; ++i) 
          {
            unsigned int index = i + j * nx + k * nx * ny;
            fprintf(fd, "%d %d %d %.6e %.6e %.6e %.6e %.6e %d %d\n", 
              i, j, k, T_e[index], S_e[index], 
              rho_e[index], C_e[index], kappa_e[index], 
              flag[index], T_dynamic_flag[index]);
          }
        }
      }

      fclose(fd);
    }
    
    void solve() 
    {
      sync_before();
      
      if(myID == 0) // solving is done only on task 0 
      {  
        // this is strongly inspired by fix_ttm
        // check for stability
        double inner_dt = dt / steps;
        
        double dtdxdydz = inner_dt * (1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz);
        
        // update temperature dependent parameters
        for(size_t i = 0; i < ntotal; ++i)
        {
          if(T_dynamic_flag[i])
          {
            C_e[i] = C_e_T(T_e[i]);
            kappa_e[i] = kappa_e_T(T_e[i]);
          }
        }
        
        /* find smallest C_e and rho_e and largest kappa */
        double c_min = C_e[0];
        double rho_min = rho_e[0];
        double kappa_max = kappa_e[0];
        
        for(size_t i = 1; i < ntotal; ++i) {
          if(flag[i] != CONSTANT_VALUE) {
            if(C_e[i] < c_min) c_min = C_e[i];
            if(rho_e[i] < rho_min) rho_min = rho_e[i];
            if(kappa_e[i] > kappa_max) kappa_max = kappa_e[i];
          }
        }
        
        double r = dtdxdydz / c_min / rho_min * kappa_max;
        
        unsigned int new_steps = steps;
        
        // This will become unstable if there are any large fluctuations 
        // during the solving process; calling this at every step is expensive
        if(r > 0.4) 
        {
          inner_dt = 0.4 * inner_dt / r; // get new stable timestep
          new_steps = std::max(static_cast<unsigned int> (dt / inner_dt), 1u);
          inner_dt = dt / new_steps;
        }
        
        for(int n = 0; n < new_steps; ++n) 
        {
          std::fill(ddT_e.begin(), ddT_e.end(), 0.0);
          
          for(unsigned int k = 0; k < nz; ++k) {
            for(unsigned int j = 0; j < ny; ++j) {
              for(unsigned int i = 0; i < nx; ++i) {                
                unsigned int q, p;
                unsigned int r = i + j*nx + k*nx*ny;
                
                if(flag[r] == ZERO_DERIVATIVE) continue;
                
                // +- dx
                if(i > 0) p = (i-1) + j*nx + k*nx*ny;
                else p = (nx-1) + j*nx + k*nx*ny;
                
                if(i < (nx - 1)) q = (i+1) + j*nx + k*nx*ny;
                else q = j*nx + k*nx*ny;
                
                if(flag[q] == ZERO_DERIVATIVE) q = r;
                else if(flag[p] == ZERO_DERIVATIVE) p = r;
                
                ddT_e[r] += (kappa_e[q]-kappa_e[p]) * (T_e[q] - T_e[p]) / dx / dx / 4.0;
                ddT_e[r] += kappa_e[r] * ((T_e[q]+T_e[p]-2.0*T_e[r]) / dx / dx);
                
                // +- dy
                if(j > 0) p = i + (j-1)*nx + k*nx*ny;
                else p = i + (ny-1)*nx + k*nx*ny;
                
                if(j < (ny - 1)) q = i + (j+1)*nx + k*nx*ny;
                else q = i + k*nx*ny;
                
                if(flag[q] == ZERO_DERIVATIVE) q = r;
                else if(flag[p] == ZERO_DERIVATIVE) p = r;
                
                ddT_e[r] += (kappa_e[q]-kappa_e[p]) * (T_e[q] - T_e[p]) / dy / dy / 4.0;
                ddT_e[r] += kappa_e[r] * ((T_e[q]+T_e[p]-2.0*T_e[r]) / dy / dy);
                
                // +- dz
                if(k > 0) p = i + j*nx + (k-1)*nx*ny;
                else p = i + j*nx + (nz-1)*nx*ny;
                
                if(k < (nz - 1)) q = i + j*nx + (k+1)*nx*ny;
                else q = i + j*nx;
                
                if(flag[q] == ZERO_DERIVATIVE) q = r;
                else if(flag[p] == ZERO_DERIVATIVE) p = r;
                
                ddT_e[r] += (kappa_e[q]-kappa_e[p]) * (T_e[q] - T_e[p]) / dz / dz / 4.0;
                ddT_e[r] += kappa_e[r] * ((T_e[q]+T_e[p]-2.0*T_e[r]) / dz / dz);
              }
            }
          }
          
          /* TODO: there might be an issue with grid volume here */
          // do the actual step
          for(int i = 0; i < ntotal; i++) {
            double prescaler = rho_e[i] * C_e[i];
            assert(prescaler > 0);
            
            switch(flag[i]) {
              case DYNAMIC:
                T_e[i] += (ddT_e[i] + dT_e[i] + S_e[i]) / prescaler * inner_dt;
                break;
              default:
                break;
            }
            
            // energy conservation issues
            /* Add a sanity check somewhere for this */
            if(T_e[i] < 0.0)
            {
              T_e[i] = 0.0;
            }
          }
        }
      }
      
      sync_after();
    }
    
  private:
    static constexpr unsigned int lineLength = 1024;
    
    size_t nx, ny, nz; // number of nodes in x,y,z
    size_t ntotal; // total number of nodes
    
    double x0, x1; // box dimensions in x
    double y0, y1; // box dimensions in y
    double z0, z1; // box dimensions in z
    
    double dx, dy, dz;
    double dV; // volume of the element
    
    std::vector<double> T_e; // current electronic temperature grid 
    std::vector<double> dT_e; // source/sink term from atoms 
    std::vector<double> ddT_e; // grid to store temporary values (almost second derivative)
    
    // temperature dependence will be added later
    std::vector<double> C_e; // specific heat at each point
    std::vector<double> rho_e; // electronic density at each point
    std::vector<double> kappa_e; // electronic heat conduction
    
    std::vector<double> S_e; // external sink and source term
    
    /*
     * -1 -> uninitialised
     *  0 -> constant
     *  1 -> dynamic
     *  2 -> derivative 0
     */
    // TODO: change into enum
    enum : signed short {
      CONSTANT_VALUE = 0,
      DYNAMIC = 1,
      ZERO_DERIVATIVE = 2
    };
    
    std::vector<signed short> flag; // node property
    
    /*
     * 0 -> no temperature dependent parameters (C_e kappa_e)
     * 1 -> temperature dependent parameters
     */
    // T_dynamic_flag
    std::vector<unsigned short> T_dynamic_flag; // temperature dependence of properties
    
    // filename for the file where temperature dependent properties are saved
    std::string parameter_filename; // NULL is special value
    
    Spline C_e_T; // temperature dependent interpolation for C_e
    Spline kappa_e_T; // tempearture dependent interpolation for kappa_e
    
    size_t steps; // number of steps 
    double dt; // value of global timestep
    
    MPI_Comm world; // communicator
    int myID;
    int nrPS;
    
    void resize_vectors(size_t in_nx, size_t in_ny, size_t in_nz)
    {
      ntotal = in_nx * in_ny * in_nz;
      
      T_e.resize(ntotal, 0);
      dT_e.resize(ntotal, 0);
      ddT_e.resize(ntotal, 0);
      
      C_e.resize(ntotal, 0);
      rho_e.resize(ntotal, 0);
      kappa_e.resize(ntotal, 0);
      
      S_e.resize(ntotal, 0);
      flag.resize(ntotal, 1);
      T_dynamic_flag.resize(ntotal, false);
    }
    
    void sync_before() // this is for MPI sync before solve is called
    {
      MPI_Allreduce(MPI_IN_PLACE, dT_e.data(), ntotal, MPI_DOUBLE, MPI_SUM, world);
    }
    
    void sync_after() // this is for MPI sync after solve is called
    {
      // zero arrays
      std::fill(dT_e.begin(), dT_e.end(), 0.0);
  
      // synchronize electronic temperature
      MPI_Bcast(T_e.data(), ntotal, MPI_DOUBLE, 0, world);
    }
    
    // possible source of error if nx*ny*nz does not fit into int
    size_t get_index(double x, double y, double z) const 
    {
      int lx = std::floor((x-x0) / dx);
      int px = std::floor( ((double) lx) / nx);
      lx -= px * nx;

      int ly = std::floor((y-y0) / dy);
      int py = std::floor( ((double) ly) / ny);
      ly -= py * ny;

      int lz = std::floor((z-z0) / dz);
      int pz = std::floor( ((double) lz) / nz);
      lz -= pz * nz;
      
      return lx + ly*nx + lz*nx*ny;
    }
    
};

#endif

