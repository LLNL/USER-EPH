/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

#ifndef EPH_FDM_H
#define EPH_FDM_H

#include <iostream>
#include <mpi.h>

/**
 * This is first stupid solution
 **/

class EPH_FDM {
  private:
    unsigned int nx; // number of nodes in x
    unsigned int ny; // number of nodes in y
    unsigned int nz; // number of nodes in z
    
    double x0, x1; // box dimensions in x
    double y0, y1; // box dimensions in y
    double z0, z1; // box dimensions in z
    double lx, ly, lz; // box size
    
    double dx, dy, dz; // element size
    double dV; // volume of the element
    
    unsigned int ntotal; // total number of nodes
    
    double *T_e; // electronic temperature grid
    double *dT_e; // energy transfer
    double *dT_e_x; // derivative in x direction
    double *dT_e_y; // derivative in y direction
    double *dT_e_z; // derivative in x direction
    double *ddT_e; // second derivative
    
    // temperature dependence will be added later
    double *C_e; // specific heat at each point; 
    double *rho_e; // electronic density at each point
    double *kappa_e; // electronic heat conduction
    
    double *S_e; // source or sink term for the electronic system
    
    /*
     * -1 -> uninitialised
     *  0 -> constant
     *  1 -> dynamic
     **/
    int *flag; // term to create walls inside the contiinum
    
    unsigned int steps; // number of steps 
    double dt; // value of global timestep
    
    MPI_Comm world; // communicator
    int myID;
    int nrPS;
    
  public:
    EPH_FDM(const char* file);
    EPH_FDM(unsigned int nx, unsigned int ny, unsigned int nz); // construct a system with 
    ~EPH_FDM();
    
    void setBox(double x0, double x1, double y0, double y1, double z0, double z1) {
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
    
    void setConstants(double T_e, double C_e, double rho_e, double kappa_e) {
      for(int i = 0; i < ntotal; i++) {
        this->T_e[i] = T_e;
        
        this->rho_e[i] = rho_e;
        this->C_e[i] = C_e;
        this->kappa_e[i] = kappa_e;
        this->flag[i] = 1;
      }
    }
    
    void setSteps(unsigned int steps) {
      this->steps = steps;
    }
    
    void setDt(double dt) {
      this->dt = dt;
    }
    
    double getT(double x, double y, double z) {
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
      
      return T_e[lx + ly * nx + lz * nx * ny];
    }
    
    double getT(unsigned int lx, unsigned int ly, unsigned int lz) {
      lx = lx%nx;
      ly = ly%ny;
      lz = lz%nz;
      
      return T_e[lx + ly * nx + lz * nx * ny];
    }
    
    bool insertEnergy(double x, double y, double z, double E) {
      unsigned int lx, ly, lz;
      double px, py, pz; // periodicity corrected
      
      /* wrap around */
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
      //printf("%f %f %f %d %d %d %d\n", px, py, pz, lx, ly, lz, index);
      double prescale = dV * dt;
      
      // convert energy into power per area
      if(prescale > 0.0)
        dT_e[index] += E / prescale;
      else
        return false;
    }
    
    void setT(double x, double y, double z, double T) {
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
    
    void setT(unsigned int lx, unsigned int ly, unsigned int lz, double T) {
      lx = lx%nx;
      ly = ly%ny;
      lz = lz%nz;
      
      unsigned int index = lx + ly * nx + lz * nx * ny;
      
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
    
    void setS(unsigned int lx, unsigned int ly, unsigned int lz, double S) {
      lx = lx%nx;
      ly = ly%ny;
      lz = lz%nz;
      
      unsigned int index = lx + ly * nx + lz * nx * ny;
      
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
    
    double calcTtotal() {
      double result = 0.0;
  
      for(int i = 0; i < ntotal; i++) {
        result += T_e[i];
      }
      result /= ntotal; // this calculates the average temperature
  
      return result;
    }
  
  private:
    static constexpr unsigned int lineLength = 1024;
    
    void init();
    void syncBefore(); // this is for MPI sync before solve is called
    void syncAfter(); // this is for MPI sync after solve is called
};

#endif
