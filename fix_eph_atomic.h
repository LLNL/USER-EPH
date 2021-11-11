/**
 * This fix is a rewrite of our previous USER-EPH (a mod of fix_ttm) code and is based on our
 * PRB 94, 024305 (2016) paper
 **/

/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

#ifdef FIX_CLASS
FixStyle(eph/atomic, FixEPHAtomic)
#else

#ifndef LMP_FIX_EPH_ATOMIC_H
#define LMP_FIX_EPH_ATOMIC_H

// external headers
#include <memory>
#include <vector>
#include <cstddef>

// lammps headers
#include "fix.h"

// internal headers
#include "eph_beta.h"
#include "eph_kappa.h"

namespace LAMMPS_NS {

class FixEPHAtomic : public Fix {
 public:
    // enumeration for tracking fix state, this is used in comm forward
    enum class FixState : unsigned int {
      NONE,
      RHO,
      XI,
      WI,
      EI // update temperatures // update energies
    };

    // enumeration for selecting fix functionality
    enum Flag : int {
      FRICTION = 0x01,
      RANDOM = 0x02,
      HEAT = 0x04,
      NOINT = 0x08, // disable integration
      NOFRICTION = 0x10, // disable effect of friction force
      NORANDOM = 0x20 // disable effect of random force
    };

    FixEPHAtomic(class LAMMPS *, int, char **); // constructor
    ~FixEPHAtomic(); // destructor

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
    
    // needed to distribute electronic energy per atom
    int pack_exchange(int, double *) override;
    int unpack_exchange(int, double*) override;
    
  protected:
    static constexpr size_t max_file_length = 256; // max filename length

    int my_id; // mpi rank for current instance
    int nr_ps; // number of processes

    FixState state; // tracks the state of the fix
    int eph_flag; // term flags remove

    int types; // number of different types
    Container<int, Allocator<int>> type_map_beta; // type map for beta file
    Container<int, Allocator<int>> type_map_kappa; // type map for kappa file

    Beta beta; // instance for beta(rho) parametrisation
    Kappa kappa; // heat diffusion parametrisation

    /** integrator functionality **/
    double dtv;
    double dtf;

    double r_cutoff; // cutoff for rho(r)
    double r_cutoff_sq;
    double rho_cutoff;

    double eta_factor; // this is for the conversion from energy/ps -> force
    double kB;
    
    int seed; // seed for random number generator
    class RanMars *random; // rng

    // Neighbor list
    class NeighList *list;

    double Ee; // energy of the electronic system
    double Te; // average electronic temperature

    size_t n; // size of peratom arrays

    // friction force
    double **f_EPH; // size = [nlocal][3]

    // random force
    double **f_RNG; // size = [nlocal][3]

    // electronic density at each atom
    double* rho_i; // size = [nlocal]

    // inverse of electronic density in order to avoid 1./rho_i
    double* inv_rho_i; // size = [nlocal]

    // dissipation vector W_ij v_j
    double** w_i; // size = [nlocal][3]

    // random numbers
    double **xi_i; // size = [nlocal][3]

    // electronic temperature per atom
    double* rho_a_i; // this specifies correlations size = [nlocal + nghost]
    double* E_a_i; // electronic energy per atom placeholder for future
    double* dE_a_i; // energy deposition by stochastic forces
    double* T_a_i; // this is for convenience

    int inner_loops; // automatic loop selection override

    // per atom array
    double **array; // size = [nlocal][8]

    // private member functions
    void calculate_environment(); // calculate the site density and coupling for every atom
    void force_prl(); // PRL model with full functionality
    void heat_solve(); // atomic heat diffusion solving
    void populate_array(); // populate per atom array with values

    // TODO: remove
    static Float get_scalar(Float const* x, Float const* y) {
      return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
    }

    static Float get_norm(Float const* x) {
      return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    }

    static Float get_distance_sq(Float const* x, Float const* y) {
      Float dxy[3];
      dxy[0] = x[0] - y[0];
      dxy[1] = x[1] - y[1];
      dxy[2] = x[2] - y[2];

      return get_norm(dxy);
    }

    static Float get_difference_sq(Float const* x, Float const* y, Float* __restrict z) {
      z[0] = x[0] - y[0];
      z[1] = x[1] - y[1];
      z[2] = x[2] - y[2];

      return get_norm(z);
    }
};

}
#endif
#endif

