/**
 * This fix is a rewrite of our previous USER-EPH (a mod of fix_ttm) code and is based on our 
 * PRB 94, 024305 (2016) paper
 **/

/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

#ifdef FIX_CLASS
FixStyle(eph,FixEPH)
#else

#ifndef LMP_FIX_EPH_H
#define LMP_FIX_EPH_H

// external headers
#include "fix.h"

// internal headers
#include "eph_beta.h"
#include "eph_fdm.h"

namespace LAMMPS_NS {

class FixEPH : public Fix {
 public:
    enum class FixState : unsigned int {
      NONE,
      XIX,
      XIY,
      XIZ,
      RHO,
      BETA,
      WX,
      WY,
      WZ
    };
    
    enum Flag : unsigned int {
      FRICTION = 0x01,
      RANDOM = 0x02,
      FDM = 0x04
    };
    
    enum Model : unsigned int {
      NONE = 0, // no friction at all (just calculates densities, gradients)
      TTM = 1, // two-temperature like model
      PRB = 2, // model in PRB 94, 024305 (2016)
      PRBMOD = 3, // random force idea
      ETA = 4, // random force idea with angular momentum
      GAP = 5, 
      GAPB = 6,
      TESTING = 9
    };
    
    FixEPH(class LAMMPS *, int, char **); // constructor
    ~FixEPH(); // destructor
    
    void init();
    void init_list(int id, NeighList *ptr);
    int setmask();
    void post_force(int);
    void end_of_step();
    void reset_dt();
    void grow_arrays(int);
    double compute_vector(int);
    double memory_usage();
    void post_run();
    
    /** integrator functionality **/
    void initial_integrate(int);
    void final_integrate();
    
    // forward communication copies information of owned local atoms to ghost
    // atoms, reverse communication does the opposite
    int pack_forward_comm(int, int *, double *, int, int *);
    void unpack_forward_comm(int, int, double *);
    //int pack_reverse_comm(int, int, double *);
    //void unpack_reverse_comm(int, int *, double *);
  
  private:
    int myID; // mpi rank for current instance
    int nrPS; // number of processes
    
    FixState state;
    
    char eph_flag; // term flags
    char eph_model; // model selection
    
    unsigned int types; // number of different types
    unsigned int* typeMap; // type map
    EPH_Beta* beta; // instance for
    EPH_FDM* fdm; // electronic FDM grid
    
    /** integrator functionality **/
    double dtv;
    double dtf;
    
    double rcutoff; // cutoff for beta
    
    // TODO: fix this; confusing names and 256
    int Tfreq;
    char Tout[256]; // this will print temperature heatmap
    char Tstate[256]; // this will store the final state into file
    
    double beta_factor; // this is for the conversion from amu/ps -> force
    double eta_factor; // this is for the conversion from energy/ps -> force
    
    int seed; // seed for random number generator
    class RanMars *random; // rng
    
    // Neighbor list
    class NeighList *list;
    
    // energy of the electronic system
    double Ee;
    
    // friction force
    double **f_EPH; // size = [nlocal][3]
    
    // random force
    double **f_RNG; // size = [nlocal][3]
    
    // stopping power for each atom
    double* beta_i; // size = [nlocal]
    
    // Electronic density at each atom
    double* rho_i; // size = [nlocal]
    
    // dissipation vector W_ij v_j
    double** w_i; // size = [nlocal][3]
    
    // locally stored individual densities to lower the number of calls to interpolation
    int rho_neigh; // 512 by default
    double** rho_ij; // size = [nlocal][rho_neigh] contribution by atom i to site j
    double** rho_ji; // size = [nlocal][rho_neigh] contribution by atom j to site i
    
    // gradient of the density
    double** grad_rho_i; // size [nlocal][3]
    
    // random numbers
    double **xi_i; // size = [nlocal][3]
    
    // per atom array
    double **array; // size = [nlocal][5]
    
    // these are temporary
    double v_alpha;
    double v_struc;
    
    // private member functions
    void calculate_environment();
    void force_ttm();
    void force_prb();
    void force_prbmod();
    void force_eta();
    void force_gap();
    void force_gapb();
};

}
#endif
#endif

