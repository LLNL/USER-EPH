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

class EPH_FDM {
  private:
    unsigned int nx; // number of nodes in x
    unsigned int ny; // number of nodes in y
    unsigned int nz; // number of nodes in z
    unsigned int ntotal; // total number of nodes
    
    double x0, x1; // box dimensions in x
    double y0, y1; // box dimensions in y
    double z0, z1; // box dimensions in z
    //double lx, ly, lz; // box size ; possible source of errors
    
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
    void setBox(double x0, double x1, double y0, double y1, double z0, double z1) {
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
      assert(dz < 0);
    }
    
    // set system properties
    void setConstants(double T_e, double C_e, double rho_e, double kappa_e) {
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
    void setSteps(unsigned int steps) {
      this->steps = steps;
    }
    
    // define timestep
    void setDt(double dt) {
      this->dt = dt;
    }
    
    // get temperature of a cell
    double getT(double x, double y, double z) const {
      unsigned int index = get_index(x, y, z);
      
      #ifndef DNDEBUG
      assert(index >= 0);
      assert(index < (nx*ny*nz));
      #endif
      return T_e[index];
    }
    
    // this has not been checked thoroughly
    double getT(unsigned int x, unsigned int y, unsigned int z) const {
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
    
    bool insertEnergy(unsigned int x, unsigned y, unsigned z, double E) {
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
    
    void setT(double x, double y, double z, double T) {
      unsigned int index = get_index(x, y, z);
      T_e[index] = T;
    }
    
    void setT(const unsigned int x, const unsigned int y, const unsigned int z, const double T) {
      unsigned int index = (x%nx) + (y%ny) * nx + (z%nz) * nx * ny;
      T_e[index] = T;
    }
    
    void setS(double x, double y, double z, double S) {
      unsigned int index = get_index(x, y, z);
      S_e[index] = S;
    }
    
    void setS(unsigned int x, unsigned int y, unsigned int z, double S) {
      unsigned int index = (x%nx) + (y%ny) * nx + (z%nz) * nx * ny;
      S_e[index] = S;
    }
    
    void setFlag(double x, double y, double z, signed char f) {
      unsigned int index = get_index(x, y, z);
      flag[index] = f;
    }
    
    void setFlag(unsigned int x, unsigned int y, unsigned int z, signed char f) {
      unsigned int index = x + y * nx + z * nx * ny;
      flag[index] = f;
    }
    
    void setComm(MPI_Comm comm, int myID, int nrPS) {
      world = comm;
      this->myID = myID;
      this->nrPS = nrPS;
    } 
    
    void solve(); // evolve electronic system
    
    // save final state of electronic structure for continuation
    void saveState(const char* file);
    
    // save current temperature map
    void saveTemperature(const char* file, int n);
    
    double calcTtotal() const {
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
    void syncBefore(); // this is for MPI sync before solve is called
    void syncAfter(); // this is for MPI sync after solve is called
};

#endif
