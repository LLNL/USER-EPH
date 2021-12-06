/*
 * Authors of the extension Artur Tamm, Alfredo Correa
 * e-mail: artur.tamm.work@gmail.com
 */

#ifdef FIX_CLASS
FixStyle(eph/coloured,FixEPHColoured)
#else

#ifndef LMP_FIX_EPH_COLOURED_H
#define LMP_FIX_EPH_COLOURED_H

// external headers
#include <memory>
#include <vector>
#include <cstddef>

// lammps headers
#include "fix.h"

// internal headers
#include "eph_beta.h"
#include "eph_fdm.h"

namespace LAMMPS_NS {
class FixEPHColoured : public Fix {
 public:
    // enumeration for tracking fix state, this is used in comm forward
    enum class FixState : unsigned int {
      NONE,
      XIX,
      XIY,
      XIZ,
      RHO,
      WX,
      WY,
      WZ,
      XI,
      WI
    };
    
    // enumeration for selecting fix functionality
    enum Flag : int {
      FRICTION = 0x01,
      RANDOM = 0x02,
      FDM = 0x04,
      NOINT = 0x08, // disable integration
      NOFRICTION = 0x10, // disable effect of friction force
      NORANDOM = 0x20 // disable effect of random force
    };
    
    FixEPHColoured(class LAMMPS *, int, char **); // constructor
    ~FixEPHColoured(); // destructor
    
    void init() override; // called by lammps after constructor
    void init_list(int id, NeighList *ptr) override; // called by lammps after constructor
    int setmask() override; // called by lammps
    void post_force(int) override; // called by lammps after pair_potential
    void end_of_step() override; // called by lammps before next step
    void reset_dt() override; // called by lammps if dt changes
    void grow_arrays(int) override; // called by lammps if number of atoms changes for some task
    double compute_vector(int) override; // called by lammps if a value is requested
    double memory_usage() override; // prints the memory usage // TODO
    void post_run() override; // called by lammps after run ends
    
    /* integrator functionality */
    void initial_integrate(int) override; // called in the beginning of a step
    void final_integrate() override; // called in the end of the step
    
    // forward communication copies information of owned local atoms to ghost
    // atoms, reverse communication does the opposite
    int pack_forward_comm(int, int *, double *, int, int *) override;
    void unpack_forward_comm(int, int, double *) override;
  
  protected:
    static constexpr size_t max_file_length = 256; // max filename length
  
    int myID; // mpi rank for current instance
    int nrPS; // number of processes
    
    FixState state; // tracks the state of the fix
    
    int eph_flag; // term flags
    int eph_model; // model selection
    
    int types; // number of different types
    int* type_map; // TODO: type map // change this to vector
    //Container<uint8_t, Allocator<uint8_t> type_map; // type map // change this to vector
    
    Beta beta; // instance for beta(rho) parametrisation
    EPH_FDM fdm; // electronic FDM grid
    
    /** integrator functionality **/
    double dtv;
    double dtf;
    
    double r_cutoff; // cutoff for rho(r)
    double r_cutoff_sq;
    double rho_cutoff; // cutoff for beta(rho)
    
    int T_freq; // frequency for printing electronic temperatures to files 
    char T_out[max_file_length]; // this will print temperature heatmap
    char T_state[max_file_length]; // this will store the final state into file
    
    double eta_factor; // this is for the conversion from energy/ps -> force
    
    int seed; // seed for random number generator
    class RanMars *random; // rng
    
    // Neighbor list
    class NeighList *list;
    
    // energy of the electronic system
    double Ee;
    
    size_t n; // size of peratom arrays
    
    // friction force
    double **f_EPH; // size = [nlocal][3] // TODO: try switching to vector
    
    // random force
    double **f_RNG; // size = [nlocal][3] // TODO: try switching to vector

    // Electronic density at each atom
    double* rho_i; // size = [nlocal] // TODO: try switching to vector
    
    // Inverse of electronic density in order to avoid 1./rho_i
    double* inv_rho_i; // size = [nlocal] // TODO: try switching to vector
    
    // dissipation vector W_ij v_j
    double** w_i; // size = [nlocal][3] // TODO: try switching to vector
    
    // random numbers
    double **xi_i; // size = [nlocal][3] // TODO: try switching to vector
    
    // electronic temperature per atom
    double* T_e_i; // size = [nlocal + nghost]
    
    // per atom array
    double **array; // size = [nlocal][8] // TODO: try switching to vector
    
    // private member functions
    void calculate_environment(); // calculate the site density and coupling for every atom
    void force_prl(); // PRL model with full functionality
    
    // TODO: remove
    static Float get_scalar(const Float* x, const Float* y) 
    {
      return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
    }
    
    static Float get_norm(const Float* x) 
    {
      return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    }
    
    static Float get_distance_sq(const Float* x, const Float* y) 
    {
      Float dxy[3];
      dxy[0] = x[0] - y[0];
      dxy[1] = x[1] - y[1];
      dxy[2] = x[2] - y[2]; 
      
      return get_norm(dxy);
    }
    
    // TODO: add restrict
    static Float get_difference_sq(const Float* x, const Float* y, Float* __restrict z) 
    {
      z[0] = x[0] - y[0];
      z[1] = x[1] - y[1];
      z[2] = x[2] - y[2];
      
      return get_norm(z);
    }
};
}
#endif
#endif

