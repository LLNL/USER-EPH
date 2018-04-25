/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

#ifndef EPH_FDM_H
#define EPH_FDM_H

#include <iostream>
#include <cassert>
#include <vector>

#include <mpi.h>

/**
 * This is first stupid solution
 **/

class EPH_FDM {
  private:
    unsigned int nx; // number of nodes in x
    unsigned int ny; // number of nodes in y
    unsigned int nz; // number of nodes in z
    unsigned int ntotal; // total number of nodes
    
    double x0, x1; // box dimensions in x
    double y0, y1; // box dimensions in y
    double z0, z1; // box dimensions in z
    double lx, ly, lz; // box size
    
    double dx, dy, dz; // element size
    double dV; // volume of the element
    
    std::vector<double> T_e; // electronic temperature at a grid point
    std::vector<double> dT_e; // energy transfer between atoms and electrons
    std::vector<double> dT_e_x; // derivative in x
    std::vector<double> dT_e_y; // derivative in y
    std::vector<double> dT_e_z; // derivative in z
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
     **/
    std::vector<int> flag; // point prorperty
    
    unsigned int steps; // number of steps 
    double dt; // value of global timestep
    
    MPI_Comm world; // communicator
    int myID;
    int nrPS;
    
  public:
    EPH_FDM() = delete;
    EPH_FDM(const char* file);
    EPH_FDM(const unsigned int nx, const unsigned int ny, const unsigned int nz); // construct a system with 
    
    void setBox(const double x0, const double x1, const double y0, const double y1, const double z0, const double z1) {
      this->x0 = x0;
      this->x1 = x1;
      this->y0 = y0;
      this->y1 = y1;
      this->z0 = z0;
      this->z1 = z1;
      
      lx = x1-x0;
      ly = y1-y0;
      lz = z1-z0;
      
      dx = lx/nx;
      dy = ly/ny;
      dz = lz/nz;
      
      dV = dx*dy*dz;
    }
    
    void setConstants(const double T_e, const double C_e, const double rho_e, const double kappa_e) {
      for(int i = 0; i < ntotal; i++) {
        this->T_e[i] = T_e;
        
        this->rho_e[i] = rho_e;
        this->C_e[i] = C_e;
        this->kappa_e[i] = kappa_e;
        this->flag[i] = 1;
      }
    }
    
    void setSteps(const unsigned int steps) {
      this->steps = steps;
    }
    
    void setDt(const double dt) {
      this->dt = dt;
    }
    
    double getT(const double x, const double y, const double z) const {
      unsigned int lx, ly, lz;
      double px, py, pz; // periodicity corrected
      
      /** TODO: this is probably slow and should be redobne **/
      if(x < x0) px = x + this->lx;
      else if( x > x1) px = x - this->lx;
      else px = x;
      
      if(y < y0) py = y + this->ly;
      else if( y > y1) py = y - this->ly;
      else py = y;
      
      if(z < z0) pz = z + this->lz;
      else if( z > z1) pz = z - this->lz;
      else pz = z;
      
      lx = (unsigned int)((px-x0)/dx+1e-12);
      ly = (unsigned int)((py-y0)/dy+1e-12);
      lz = (unsigned int)((pz-z0)/dz+1e-12);
      
      lx = lx%nx;
      ly = ly%ny;
      lz = lz%nz;
      
      unsigned int index = lx + ly * nx + lz * nx * ny;
      //printf("%f %f %f %d %d %d %d\n", x, y, z, lx, ly, lz, index);
      assert(index >= 0);
      assert(index < (nx*ny*nz));
      return T_e[index];
    }
    
    double getT(const unsigned int x, const unsigned int y, const unsigned int z) const {
      unsigned int index = (x%nx) + (y%ny) * nx + (z%nz) * nx * ny;
      return T_e[index];
    }
    
    bool insertEnergy(const double x, const double y, const double z, const double E) {
      unsigned int lx, ly, lz;
      double px, py, pz; // periodicity corrected
      
      /* wrap around */
      /** TODO: this is probably slow and should be redobne **/
      if(x < x0) px = x + this->lx;
      else if( x > x1) px = x - this->lx;
      else px = x;
      
      if(y < y0) py = y + this->ly;
      else if( y > y1) py = y - this->ly;
      else py = y;
      
      if(z < z0) pz = z + this->lz;
      else if( z > z1) pz = z - this->lz;
      else pz = z;
      
      lx = (unsigned int)((px-x0)/dx+1e-12);
      ly = (unsigned int)((py-y0)/dy+1e-12);
      lz = (unsigned int)((pz-z0)/dz+1e-12);
      
      lx = lx%nx;
      ly = ly%ny;
      lz = lz%nz;
      
      unsigned int index = lx + ly * nx + lz * nx * ny;
      //printf("%f %f %f %d %d %d %d\n", px, py, pz, lx, ly, lz, index);
      double prescale = dV * dt;
      
      // convert energy into power per area
      if(prescale > 0.0) {
        assert(index >= 0);
        assert(index < (nx*ny*nz));
        dT_e[index] += E / prescale;
      }
      else
        return false;
    }
    
    void setT(const double x, const double y, const double z, const double T) {
      unsigned int lx, ly, lz;
      double px, py, pz; // periodicity corrected
      
      if(x < x0) px = x + this->lx;
      else px = x;
      
      if(y < y0) py = y + this->ly;
      else py = y;
      
      if(z < z0) pz = z + this->lz;
      else pz = z;
      
      lx = (unsigned int)((px-x0)/dx+1e-12);
      ly = (unsigned int)((py-y0)/dy+1e-12);
      lz = (unsigned int)((pz-z0)/dz+1e-12);
      
      lx = lx%nx;
      ly = ly%ny;
      lz = lz%nz;
      
      unsigned int index = lx + ly * nx + lz * nx * ny;
      
      T_e[index] = T;
    }
    
    void setT(const unsigned int x, const unsigned int y, const unsigned int z, const double T) {
      unsigned int index = (x%nx) + (y%ny) * nx + (z%nz) * nx * ny;
      
      T_e[index] = T;
    }
    
    void setS(double x, double y, double z, double S) {
      unsigned int lx, ly, lz;
      double px, py, pz; // periodicity corrected
      
      if(x < x0) px = x + this->lx;
      else px = x;
      
      if(y < y0) py = y + this->ly;
      else py = y;
      
      if(z < z0) pz = z + this->lz;
      else pz = z;
      
      lx = (unsigned int)((px-x0)/dx+1e-12);
      ly = (unsigned int)((py-y0)/dy+1e-12);
      lz = (unsigned int)((pz-z0)/dz+1e-12);
      
      lx = lx%nx;
      ly = ly%ny;
      lz = lz%nz;
      
      unsigned int index = lx + ly * nx + lz * nx * ny;
      
      S_e[index] = S;
    }
    
    void setS(const unsigned int x, const unsigned int y, const unsigned int z, double S) {
      unsigned int index = (x%nx) + (y%ny) * nx + (z%nz) * nx * ny;
      
      S_e[index] = S;
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
      double result = 0.0;
  
      for(int i = 0; i < ntotal; i++) {
        result += T_e[i];
      }
      result /= ntotal; // this calculates the average temperature
  
      return result;
    }
  
  private:
    static constexpr unsigned int lineLength = 1024;
    
    void resize_vectors();
    void syncBefore(); // this is for MPI sync before solve is called
    void syncAfter(); // this is for MPI sync after solve is called
};

#endif
