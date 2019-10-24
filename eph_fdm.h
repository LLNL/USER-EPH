/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

#ifndef EPH_FDM_H
#define EPH_FDM_H

#include <iostream>
#include <cassert>
#include <vector>
#include <stdexcept>
#include <cmath>

#include <mpi.h>

class EPH_FDM
{
  public:
    // default constructor to create a grid with one point
    EPH_FDM() : EPH_FDM(1, 1, 1) {}
    
    EPH_FDM(size_t in_nx, size_t in_ny, size_t in_nz) :
      nx {in_nx},
      ny {in_ny},
      nz {in_nz}
    {
      resize_vectors();
      set_box_dimensions(0, 1, 0, 1, 0, 1);
      set_constants(0, 1, 1, 1);
      set_dt(1);
    }
    
    EPH_FDM(const char *fn) 
    {
      
      
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
      std::fill(T_flag.begin(), T_flag.end(), false);
    }
    
    void set_dt(double in_dt) 
    {
      dt = in_dt;
    }
    
    void set_steps(size_t in_steps)
    {
      steps = in_steps;
    }
    
    void set_comm(MPI_Comm in_comm, int in_myID, int in_nrPS) {
      world = in_comm;
      myID = in_myID;
      nrPS = in_nrPS;
    } 
    
    void save_temperature(const char* file, int n) {
      char fn[512];
      sprintf(fn, "%s_%06d", file, n);
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
    
  private:
    size_t nx, ny, nz; // number of nodes in x,y,z
    size_t ntotal; // total number of nodes
    
    double x0, x1; // box dimensions in x
    double y0, y1; // box dimensions in y
    double z0, z1; // box dimensions in z
    
    double dx, dy, dz;
    double dV; // volume of the element
    
    std::vector<double> T_e;
    std::vector<double> dT_e;
    std::vector<double> ddT_e;
    
    // temperature dependence will be added later
    std::vector<double> C_e; // specific heat at each point
    std::vector<double> rho_e; // electronic density at each point
    std::vector<double> kappa_e; // electronic heat conduction
    
    std::vector<double> S_e; // sink and source term
    
    /*
     * -1 -> uninitialised
     *  0 -> constant
     *  1 -> dynamic
     *  2 -> derivative 0
     **/
    std::vector<signed short> flag; // node property
    std::vector<bool> T_flag; // temperature dependence of properties
    
    size_t steps; // number of steps 
    double dt; // value of global timestep
    
    MPI_Comm world; // communicator
    int myID;
    int nrPS;
    
    void resize_vectors()
    {
      ntotal = nx * ny * nz;
      
    }
    
    void sync_before() // this is for MPI sync before solve is called
    {
      if(nrPS > 0)
        MPI_Allreduce(MPI_IN_PLACE, dT_e.data(), ntotal, MPI_DOUBLE, MPI_SUM, world);
    }
    
    void sync_after() // this is for MPI sync after solve is called
    {
      // zero arrays
      std::fill(dT_e.begin(), dT_e.end(), 0.0);
  
      // synchronize electronic temperature
      if(nrPS > 0)
        MPI_Bcast(T_e.data(), ntotal, MPI_DOUBLE, 0, world);
    }
};


#if 0
class EPH_FDM {
  private:
    unsigned int nx; // number of nodes in x
    unsigned int ny; // number of nodes in y
    unsigned int nz; // number of nodes in z
    unsigned int ntotal; // total number of nodes
    
    double x0, x1; // box dimensions in x
    double y0, y1; // box dimensions in y
    double z0, z1; // box dimensions in z
    
    double dx, dy, dz; // element size
    double dV; // volume of the element
    
    std::vector<double> T_e; // electronic temperature at a grid point
    std::vector<double> dT_e; // energy transfer between atoms and electrons
    std::vector<double> ddT_e; // second derivative
    
    // temperature dependence will be added later
    std::vector<double> C_e; // specific heat at each point
    std::vector<double> rho_e; // electronic density at each point
    std::vector<double> kappa_e; // electronic heat conduction
    
    std::vector<double> S_e; // sink and source term
    
    /*
     * -1 -> uninitialised
     *  0 -> constant
     *  1 -> dynamic
     *  2 -> derivative 0
     **/
    std::vector<signed short> flag; // point property
    std::vector<bool> T_flag; // temperature dependence

    unsigned int steps; // number of steps 
    double dt; // value of global timestep
    
    MPI_Comm world; // communicator
    int myID;
    int nrPS;
    
  public:
    EPH_FDM() : nx {0}, ny {0}, nz {0}, ntotal {0} {}; // no default constructor
    EPH_FDM(const char* file); // initialise class from an input file
    EPH_FDM(unsigned int nx, unsigned int ny, unsigned int nz); // initialise class manually
    
    // set box size
    void set_box(double x0, double x1, double y0, double y1, double z0, double z1) {
      this->x0 = x0;
      this->x1 = x1;
      this->y0 = y0;
      this->y1 = y1;
      this->z0 = z0;
      this->z1 = z1;
      
      dx = (x1-x0)/nx;
      dy = (y1-y0)/ny;
      dz = (z1-z0)/nz;
      
      dV = dx*dy*dz;

      assert(dx > 0);
      assert(dy > 0);
      assert(dz > 0);
    }
    
    // set system properties
    void set_constants(double T_e, double C_e, double rho_e, double kappa_e) {
      for(int i = 0; i < ntotal; i++) {
        this->T_e[i] = T_e;
        
        this->rho_e[i] = rho_e;
        this->C_e[i] = C_e;
        this->kappa_e[i] = kappa_e;
        this->flag[i] = 1;
      }
    }
    
    // individual setters
    void set_C_e(unsigned int x, unsigned int y, unsigned int z, double C_e) {
      unsigned int index = (x%nx) + (y%ny) * nx + (z%nz) * nx * ny;
      this->C_e[index] = C_e;
    }
    
    void set_rho_e(unsigned int x, unsigned int y, unsigned int z, double rho_e) {
      unsigned int index = (x%nx) + (y%ny) * nx + (z%nz) * nx * ny;
      this->rho_e[index] = rho_e;
    }
    
    void set_kappa_e(unsigned int x, unsigned int y, unsigned int z, double kappa_e) {
      unsigned int index = (x%nx) + (y%ny) * nx + (z%nz) * nx * ny;
      this->kappa_e[index] = kappa_e;
    }
    
    // define number of minimum steps for integration
    void set_steps(unsigned int steps) {
      this->steps = steps;
    }
    
    // define timestep
    void set_dt(double in_dt) {
      dt = in_dt;
    }
    
    // get temperature of a cell
    double get_T(double x, double y, double z) const {
      unsigned int index = get_index(x, y, z);
      
      #ifndef DNDEBUG
      assert(index >= 0);
      assert(index < (nx*ny*nz));
      #endif
      return T_e[index];
    }
    
    // this has not been checked thoroughly
    double get_T(unsigned int x, unsigned int y, unsigned int z) const {
      unsigned int index = (x%nx) + (y%ny) * nx + (z%nz) * nx * ny;
      return T_e[index];
    }
    
    bool insertEnergy(double x, double y, double z, double E) {
      unsigned int index = get_index(x, y, z);
      double prescale = dV * dt;
      
      // convert energy into power per area
      #ifndef DNDEBUG
      assert(prescale >= 0.0);
      assert(index >= 0);
      assert(index < (nx*ny*nz));
      #endif
      dT_e[index] += E / prescale;
      
      return true;
    }
    
    bool insert_Energy(unsigned int x, unsigned y, unsigned z, double E) {
      unsigned int index = x + y*nx + z*nx*ny;
      double prescale = dV * dt;
      
      // convert energy into power per area
      #ifndef DNDEBUG
      assert(prescale >= 0.0);
      assert(index >= 0);
      assert(index < (nx*ny*nz));
      #endif
      dT_e[index] += E / prescale;
      
      return true;
    }
    
    void set_T(double x, double y, double z, double T) {
      unsigned int index = get_index(x, y, z);
      T_e[index] = T;
    }
    
    void set_T(const unsigned int x, const unsigned int y, const unsigned int z, const double T) {
      unsigned int index = (x%nx) + (y%ny) * nx + (z%nz) * nx * ny;
      T_e[index] = T;
    }
    
    void set_S(double x, double y, double z, double S) {
      unsigned int index = get_index(x, y, z);
      S_e[index] = S;
    }
    
    void set_S(unsigned int x, unsigned int y, unsigned int z, double S) {
      unsigned int index = (x%nx) + (y%ny) * nx + (z%nz) * nx * ny;
      S_e[index] = S;
    }
    
    void set_Flag(double x, double y, double z, signed char f) {
      unsigned int index = get_index(x, y, z);
      flag[index] = f;
    }
    
    void set_Flag(unsigned int x, unsigned int y, unsigned int z, signed char f) {
      unsigned int index = x + y * nx + z * nx * ny;
      flag[index] = f;
    }
    
    void set_Comm(MPI_Comm comm, int myID, int nrPS) {
      world = comm;
      this->myID = myID;
      this->nrPS = nrPS;
    } 
    
    void solve(); // evolve electronic system
    
    // save final state of electronic structure for continuation
    void save_State(const char* file);
    
    // save current temperature map
    void save_Temperature(const char* file, int n);
    
    double calc_Ttotal() const {
      double result {0.0};
  
      for(int i = 0; i < ntotal; i++) {
        result += T_e[i];
      }
      result /= ntotal; // this calculates the average temperature
  
      return result;
    }
  
  private:
    static constexpr unsigned int lineLength = 1024;
    
    // possible source of error if nx*ny*nz does not fit into int
    unsigned int get_index(double x, double y, double z) const {
      int lx = std::floor((x-x0) / dx);
      int px = std::floor( ((double) lx) / nx);
      lx -= px*nx;

      int ly = std::floor((y-y0) / dy);
      int py = std::floor( ((double) ly) / ny);
      ly -= py * ny;

      int lz = std::floor((z-z0) / dz);
      int pz = std::floor( ((double) lz) / nz);
      lz -= pz * nz;
      
      return lx + ly*nx + lz*nx*ny;
    }

    void resize_vectors();
    void sync_before(); // this is for MPI sync before solve is called
    void sync_after(); // this is for MPI sync after solve is called
};
#endif
#endif
