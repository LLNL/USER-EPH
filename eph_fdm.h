/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: arturt@ut.ee
 */

#ifndef EPH_FDM_H
#define EPH_FDM_H

#include <iostream>

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
    double C_e; // specific heat at each point; 
    double rho_e; // electronic density at each point
    double kappa_e; // electronic heat conduction
    
    double *S_e; // source or sink term for the electronic system
    unsigned int *flag; // term to create walls inside the contiinum
    
    unsigned int steps; // number of steps 
    double dt; // value of global timestep
    
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
      
      dx = (x1-x0)/nx;
      dy = (y1-y0)/ny;
      dz = (z1-z0)/nz;
      
      dV = dx*dy*dz;
    }
    
    void setConstants(double T_e, double C_e, double rho_e, double kappa_e) {
      this->C_e = C_e;
      this->rho_e = rho_e;
      this->kappa_e = kappa_e;
        
      for(int i = 0; i < ntotal; i++) {
        this->T_e[i] = T_e;
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
      
      lx = (unsigned int)((x-x0)/dx+1e-12);
      ly = (unsigned int)((y-y0)/dy+1e-12);
      lz = (unsigned int)((z-z0)/dz+1e-12);
      
      lx = lx%nx;
      ly = ly%ny;
      lz = lz%nz;
      
      return T_e[lx + ly * nx + lz * nx * ny];
    }
    
    double getT(int lx, int ly, int lz) {
      lx = lx%nx;
      ly = ly%ny;
      lz = lz%nz;
      
      return T_e[lx + ly * nx + lz * nx * ny];
    }
    
    bool insertEnergy(double x, double y, double z, double E) {
      unsigned int lx, ly, lz;
      
      lx = (unsigned int)((x-x0)/dx+1e-12);
      ly = (unsigned int)((y-y0)/dy+1e-12);
      lz = (unsigned int)((z-z0)/dz+1e-12);
      
      unsigned int index = lx + ly * nx + lz * nx * ny;
      double prescale = rho_e * C_e * dV;
      
      if(prescale > 0.0)
        dT_e[index] += E / prescale;
      else
        return false;
    }
    
    void setT(double x, double y, double z, double T) {
      unsigned int lx, ly, lz;
      
      lx = (unsigned int)((x-x0)/dx+1e-12);
      ly = (unsigned int)((y-y0)/dy+1e-12);
      lz = (unsigned int)((z-z0)/dz+1e-12);
      
      lx = lx%nx;
      ly = ly%ny;
      lz = lz%nz;
      
      unsigned int index = lx + ly * nx + lz * nx * ny;
      
      T_e[index] = T;
    }
    
    void setT(int lx, int ly, int lz, double T) {
      lx = lx%nx;
      ly = ly%ny;
      lz = lz%nz;
      
      unsigned int index = lx + ly * nx + lz * nx * ny;
      
      T_e[index] = T;
    }
    
    void setS(double x, double y, double z, double S) {
      unsigned int lx, ly, lz;
      
      lx = (unsigned int)((x-x0)/dx+1e-12);
      ly = (unsigned int)((y-y0)/dy+1e-12);
      lz = (unsigned int)((z-z0)/dz+1e-12);
      
      lx = lx%nx;
      ly = ly%ny;
      lz = lz%nz;
      
      unsigned int index = lx + ly * nx + lz * nx * ny;
      
      S_e[index] = S;
    }
    
    void setS(int lx, int ly, int lz, double S) {
      lx = lx%nx;
      ly = ly%ny;
      lz = lz%nz;
      
      unsigned int index = lx + ly * nx + lz * nx * ny;
      
      S_e[index] = S;
    }

    unsigned int solve(); // evolve electronic system
    
    void saveState(const char* file);
    
    double calcEtotal();
  
  private:
    static constexpr unsigned int lineLength = 1024;
    
    void init();
    void syncBefore(); // this is for MPI sync before solve is called
    void syncAfter(); // this is for MPI sync after solve is called
};

#endif
