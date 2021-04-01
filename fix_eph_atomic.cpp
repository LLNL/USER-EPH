/*
 * Authors of the extension Artur Tamm, Alfredo Correa
 * e-mail: artur.tamm.work@gmail.com
 */

// external headers
#include <iostream>
#include <cstring> // TODO: remove
#include <string>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <cassert>

// lammps headers
#include "error.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "atom.h"
#include "memory.h"
#include "random_mars.h"
#include "force.h"
#include "update.h"
#include "comm.h"

// internal headers
#include "fix_eph_atomic.h"
#include "eph_beta.h"
#include "eph_kappa.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/*
 * TODO: implement write restart
 */

/**
   * FixEPH arguments
   * arg[ 0] <- fix ID
   * arg[ 1] <- group
   * arg[ 2] <- name
   * arg[ 3] <- rng seed
   * arg[ 4] <- eph parameter; 0 disable all terms; 1 enable friction term; 2 enable random force; 4 enable fde;
   * arg[ 5] <- initial electronic temperature
   * arg[ 6] <- input file for initial temperatures; how do we load temperatures?
   * arg[ 7] <- number of inner loops n < 1 -> automatic selection
   * arg[ 8] <- output file for temperatures
   * arg[ 9] <- input file for eph model functions
   * arg[10] <- input file for kappa model functions
   * arg[11] <- element name for type 0
   * arg[12] <- element name for type 1
   * ...
   **/

// constructor
FixEPHAtomic::FixEPHAtomic(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if (narg < 12) error->all(FLERR, "fix_eph_atomic: too few arguments");
  if (atom->natoms < 1) error->all(FLERR, "fix_eph_atomic: error no atoms in simulation");
  MPI_Comm_rank(world, &my_id);
  MPI_Comm_size(world, &nr_ps);

  state = FixState::NONE;
  { // setup fix properties
    vector_flag = 1; // fix is able to output a vector compute
    size_vector = 2; // 2 elements in the vector
    global_freq = 1; // frequency for vector data
    extvector = 1; // external vector allocated by this fix???
    nevery = 1; // call end_of_step every step
    peratom_flag = 1; // fix provides per atom values
    size_peratom_cols = 11; // per atom has 8 dimensions
    peratom_freq = 1; // per atom values are provided every step
    //ghostneigh = 1; // neighbours of neighbours

    comm_forward = 1; // forward communication is needed
    comm->ghost_velocity = 1; // special: fix requires velocities for ghost atoms
  }

  { // setup arrays
    f_EPH = nullptr;
    f_RNG = nullptr;

    w_i = nullptr;

    rho_i = nullptr;
    array = nullptr;

    xi_i = nullptr;

    rho_a_i = nullptr;
    E_a_i = nullptr;
    dE_a_i = nullptr;

    list = nullptr;

    // NO ARRAYS BEFORE THIS
    grow_arrays(atom->nmax);
    atom->add_callback(0);

    // zero arrays, so they would not contain garbage
    size_t nlocal = atom->nlocal;
    size_t ntotal = atom->nghost + nlocal;

    std::fill_n(&(rho_i[0]), ntotal, 0);
    std::fill_n(&(xi_i[0][0]), 3 * ntotal, 0);
    std::fill_n(&(w_i[0][0]), 3 * ntotal, 0);

    std::fill_n(&(f_EPH[0][0]), 3 * ntotal, 0);
    std::fill_n(&(f_RNG[0][0]), 3 * ntotal, 0);

    std::fill_n(&(array[0][0]), size_peratom_cols * ntotal, 0);

    std::fill_n(&(rho_a_i[0]), ntotal, 0);
    std::fill_n(&(E_a_i[0]), ntotal, 0);
    std::fill_n(&(dE_a_i[0]), ntotal, 0);
  }

  {
    types = atom->ntypes;
    eta_factor = sqrt(2.0 * force->boltz / update->dt);
  }

  { /** integrator functionality **/
    dtv = update->dt;
    dtf = 0.5 * update->dt * force->ftm2v;
  }

  // initialise rng
  seed = atoi(arg[3]);
  random = new RanMars(lmp, seed + my_id);

  // read model behaviour parameters
  eph_flag = strtol(arg[4], NULL, 0);

  // print enabled fix functionality
  if(my_id == 0) {
    std::cout << '\n';
    std::cout << "Flag read: " << arg[4] << " -> " << eph_flag << '\n';
    if(eph_flag & Flag::FRICTION) { std::cout << "Friction evaluation: ON\n"; }
    if(eph_flag & Flag::RANDOM) { std::cout << "Random evaluation: ON\n"; }
    if(eph_flag & Flag::HEAT) { std::cout << "Heat diffusion solving: ON\n"; }
    if(eph_flag & Flag::NOINT) { std::cout << "No integration: ON\n"; }
    if(eph_flag & Flag::NOFRICTION) { std::cout << "No friction application: ON\n"; }
    if(eph_flag & Flag::NORANDOM) { std::cout << "No random application: ON\n"; }
    std::cout << '\n';
  }

  // argument 5  and 6 are handled below

  { // setup automagic inner loops or not
    inner_loops = atoi(arg[7]);

    if(inner_loops < 1) { inner_loops = 0; }
  }

  { // setup output
    if(strcmp("NULL" , arg[8]) == 0) { } // do nothing for now
  }

  int n_elem = 11; // location where element names start

  if(types > (narg - n_elem)) {
    error->all(FLERR, "fix_eph_atomic: number of types larger than provided in fix");
  }

  type_map_beta.resize(types);
  type_map_kappa.resize(types);

  beta = Beta(arg[9]);
  kappa = Kappa(arg[10]);

  if(beta.get_n_elements() < 1) {
    error->all(FLERR, "fix_eph_atomic: no elements found in beta file");
  }

  if(kappa.n_elements < 1) {
    error->all(FLERR, "fix_eph_atomic: no elements found in kappa file");
  }

  r_cutoff = beta.get_r_cutoff();
  r_cutoff_sq = beta.get_r_cutoff_sq();
  rho_cutoff = beta.get_rho_cutoff();

  // do element mapping for beta
  for(size_t i = 0; i < types; ++i) {
    type_map_beta[i] = std::numeric_limits<int>::max();
    type_map_kappa[i] = std::numeric_limits<int>::max();

    for(size_t j = 0; j < beta.get_n_elements(); ++j) {
      if((beta.get_element_name(j)).compare(arg[n_elem + i]) == 0) {
        type_map_beta[i] = j;
        break;
      }
    }

    for(size_t j = 0; j < kappa.n_elements; ++j) {
      if((kappa.element_name[j]).compare(arg[n_elem + i]) == 0) {
        type_map_kappa[i] = j;
        break;
      }
    }

    if(type_map_beta[i] > types || type_map_kappa[i] > types) {
      error->all(FLERR, "fix_eph_atomic: elements not found in input file");
    }
  }

  { // setup temperatures per atom
    double v_Te = atof(arg[5]);

    if(strcmp("NULL" , arg[6]) == 0) {
      for(size_t i = 0; i < atom->nlocal; ++i) {
        E_a_i[i] = kappa.T_E_atomic[type_map_kappa[atom->type[i] - 1]].reverse(v_Te);
      }
    }
    //~ else {
      //~ std::fill_n(&(T_a_i[0]), atom->nlocal, v_Te); // TODO: placeholder
    //~ }
  }

  // I think we will switch to keeping track of energy instead of tracking temperature
  Ee = 0.0; // electronic energy is zero in the beginning
}

// destructor
FixEPHAtomic::~FixEPHAtomic() {
  delete random;

  atom->delete_callback(id, 0);
  memory->destroy(rho_i);

  memory->destroy(array);

  memory->destroy(f_EPH);
  memory->destroy(f_RNG);
  memory->destroy(xi_i);
  memory->destroy(w_i);

  memory->destroy(rho_a_i);
  memory->destroy(E_a_i);
  memory->destroy(dE_a_i);
}

void FixEPHAtomic::init() {
  if (domain->dimension == 2) {
    error->all(FLERR,"Cannot use fix eph with 2d simulation");
  }
  if (domain->nonperiodic != 0) {
    error->all(FLERR,"Cannot use nonperiodic boundares with fix eph");
  }
  if (domain->triclinic) {
    error->all(FLERR,"Cannot use fix eph with triclinic box");
  }

  /* copy paste from vcsgc */
  /** we are a fix and we need full neighbour list **/
  int irequest = neighbor->request((void*)this, this->instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->ghost = 1;

  neighbor->requests[irequest]->cutoff = r_cutoff;

  reset_dt();
}

void FixEPHAtomic::init_list(int id, NeighList *ptr) {
  this->list = ptr;
}

int FixEPHAtomic::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  /* integrator functionality */
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;

  return mask;
}

/* integrator functionality */
void FixEPHAtomic::initial_integrate(int) {
  if(eph_flag & Flag::NOINT) return;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (size_t i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      double dtfm = dtf / mass[type[i]];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
  }
}

void FixEPHAtomic::final_integrate() {
  if(eph_flag & Flag::NOINT) return;

  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (size_t i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      double dtfm = dtf / mass[type[i]];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
    }
  }
}

void FixEPHAtomic::end_of_step() {
  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double E_local = 0.0;

  // calculate the energy transferred to electronic system
  // this is a potential source of errors due to using velocity verlet for integration
  // friction force depends on the velocities and therefore acceleration at
  //   next timestep depends on the velocity at next time step
  // this leads to errors of the order of dt^2
  if(eph_flag & Flag::FRICTION) {
    for(size_t i = 0; i < nlocal; ++i) {
      if(mask[i] & groupbit) {
        double dE_i = 0.0;
        dE_i -= f_EPH[i][0] * v[i][0] * update->dt;
        dE_i -= f_EPH[i][1] * v[i][1] * update->dt;
        dE_i -= f_EPH[i][2] * v[i][2] * update->dt;

        //~ fdm.insert_energy(x[i][0], x[i][1], x[i][2], dE_i);
        E_local += dE_i;
      }
    }
  }

  if(eph_flag & Flag::RANDOM) {
    for(size_t i = 0; i < nlocal; ++i) {
      if(mask[i] & groupbit) {
        double dE_i = 0.0;
        dE_i -= f_RNG[i][0] * v[i][0] * update->dt;
        dE_i -= f_RNG[i][1] * v[i][1] * update->dt;
        dE_i -= f_RNG[i][2] * v[i][2] * update->dt;

        //~ fdm.insert_energy(x[i][0], x[i][1], x[i][2], dE_i);
        E_local += dE_i;
      }
    }
  }

  if(eph_flag & Flag::HEAT) {
    heat_solve();
  }

  // this is for checking energy conservation
  MPI_Allreduce(MPI_IN_PLACE, &E_local, 1, MPI_DOUBLE, MPI_SUM, world);

  Ee += E_local;

  for(size_t i = 0; i < nlocal; ++i) {
    if(mask[i] & groupbit) {
      int itype = type[i];
      array[i][ 0] = rho_i[i];
      array[i][ 1] = beta.get_beta(type_map_beta[itype - 1], rho_i[i]);
      array[i][ 2] = f_EPH[i][0];
      array[i][ 3] = f_EPH[i][1];
      array[i][ 4] = f_EPH[i][2];
      array[i][ 5] = f_RNG[i][0];
      array[i][ 6] = f_RNG[i][1];
      array[i][ 7] = f_RNG[i][2];
      array[i][ 8] = rho_a_i[i];
      array[i][ 9] = E_a_i[i];
      array[i][10] = dE_a_i[i];
    }
    else {
      array[i][ 0] = 0.0;
      array[i][ 1] = 0.0;
      array[i][ 2] = 0.0;
      array[i][ 3] = 0.0;
      array[i][ 4] = 0.0;
      array[i][ 5] = 0.0;
      array[i][ 6] = 0.0;
      array[i][ 7] = 0.0;
      array[i][ 8] = 0.0;
      array[i][ 9] = 0.0;
      array[i][10] = 0.0;
    }
  }
}

void FixEPHAtomic::calculate_environment() {
  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // loop over atoms and their neighbours and calculate rho and beta(rho)
  for(size_t i = 0; i != nlocal; ++i)
  {
    rho_i[i] = 0;
    rho_a_i[i] = 0;

    // check if current atom belongs to fix group and if an atom is local
    if(mask[i] & groupbit)
    {
      int itype = type[i];
      int *jlist = firstneigh[i];
      int jnum = numneigh[i];

      for(size_t j = 0; j != jnum; ++j) {
        int jj = jlist[j];
        jj &= NEIGHMASK;

        int jtype = type[jj];
        double r_sq = get_distance_sq(x[jj], x[i]);

        if(r_sq < r_cutoff_sq) {
          rho_i[i] += beta.get_rho_r_sq(type_map_beta[jtype - 1], r_sq);
        }

        if(r_sq < kappa.r_cutoff_sq) {
          rho_a_i[i] += kappa.rho_r_sq[type_map_kappa[jtype - 1]](r_sq);
        }
      }
    }
  }
}

void FixEPHAtomic::force_prl() {
  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // create friction forces
  if(eph_flag & Flag::FRICTION)
  {
    // w_i = W_ij^T v_j
    for(size_t i = 0; i != nlocal; ++i)
    {
      if(mask[i] & groupbit) {
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];

        if(!(rho_i[i] > 0)) continue;

        double alpha_i = beta.get_alpha(type_map_beta[itype - 1], rho_i[i]);

        for(size_t j = 0; j != jnum; ++j)
        {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];

          // calculate the e_ij vector TODO: change these
          double e_ij[3];
          double e_r_sq = get_difference_sq(x[jj], x[i], e_ij);

          // first sum
          if(e_r_sq >= r_cutoff_sq) continue;

          double v_rho_ji = beta.get_rho_r_sq(type_map_beta[jtype - 1], e_r_sq);
          double prescaler = alpha_i * v_rho_ji / (rho_i[i] * e_r_sq);

          double e_v_v1 = get_scalar(e_ij, v[i]);
          double var1 = prescaler * e_v_v1;

          double e_v_v2 = get_scalar(e_ij, v[jj]);
          double var2 = prescaler * e_v_v2;

          double dvar = var1 - var2;
          w_i[i][0] += dvar * e_ij[0];
          w_i[i][1] += dvar * e_ij[1];
          w_i[i][2] += dvar * e_ij[2];
        }
      }
    }

    state = FixState::WX;
    comm->forward_comm_fix(this);
    state = FixState::WY;
    comm->forward_comm_fix(this);
    state = FixState::WZ;
    comm->forward_comm_fix(this);

    // now calculate the forces
    // f_i = W_ij w_j
    for(size_t i = 0; i != nlocal; ++i) {
      if(mask[i] & groupbit) {
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];

        if( not(rho_i[i] > 0) ) continue;

        double alpha_i = beta.get_alpha(type_map_beta[itype - 1], rho_i[i]);

        for(size_t j = 0; j != jnum; ++j)
        {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];

          // calculate the e_ij vector
          double e_ij[3];
          double e_r_sq = get_difference_sq(x[jj], x[i], e_ij);

          if(e_r_sq >= r_cutoff_sq or not(rho_i[jj] > 0)) continue;

          double alpha_j = beta.get_alpha(type_map_beta[jtype - 1], rho_i[jj]);

          double v_rho_ji = beta.get_rho_r_sq(type_map_beta[jtype - 1], e_r_sq);
          double e_v_v1 = get_scalar(e_ij, w_i[i]);
          double var1 = alpha_i * v_rho_ji * e_v_v1 / (rho_i[i] * e_r_sq);

          double v_rho_ij = beta.get_rho_r_sq(type_map_beta[itype - 1], e_r_sq);
          double e_v_v2 = get_scalar(e_ij, w_i[jj]);
          double var2 = alpha_j * v_rho_ij * e_v_v2 / (rho_i[jj] * e_r_sq);

          double dvar = var1 - var2;
          // friction is negative!
          f_EPH[i][0] -= dvar * e_ij[0];
          f_EPH[i][1] -= dvar * e_ij[1];
          f_EPH[i][2] -= dvar * e_ij[2];
        }
      }
    }
  }

  // create random forces
  if(eph_flag & Flag::RANDOM) {
    for(size_t i = 0; i != nlocal; i++) {
      if(mask[i] & groupbit) {
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];

        if(!(rho_i[i] > 0)) continue;

        double alpha_i = beta.get_alpha(type_map_beta[itype - 1], rho_i[i]);

        for(size_t j = 0; j != jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];

          // calculate the e_ij vector
          double e_ij[3];
          double e_r_sq = get_difference_sq(x[jj], x[i], e_ij);

          if((e_r_sq >= r_cutoff_sq) || !(rho_i[jj] > 0)) continue;

          double alpha_j = beta.get_alpha(type_map_beta[jtype - 1], rho_i[jj]);

          double v_rho_ji = beta.get_rho_r_sq(type_map_beta[jtype - 1], e_r_sq);
          double e_v_xi1 = get_scalar(e_ij, xi_i[i]);
          double var1 = alpha_i * v_rho_ji * e_v_xi1 / (rho_i[i] * e_r_sq);

          double v_rho_ij = beta.get_rho_r_sq(type_map_beta[itype - 1], e_r_sq);
          double e_v_xi2 = get_scalar(e_ij, xi_i[jj]);
          double var2 = alpha_j * v_rho_ij * e_v_xi2 / (rho_i[jj] * e_r_sq);

          double dvar = var1 - var2;
          f_RNG[i][0] += dvar * e_ij[0];
          f_RNG[i][1] += dvar * e_ij[1];
          f_RNG[i][2] += dvar * e_ij[2];
        }

        double v_Te = kappa.T_E_atomic[itype - 1](E_a_i[i]);
        double var = eta_factor * sqrt(v_Te);
        f_RNG[i][0] *= var;
        f_RNG[i][1] *= var;
        f_RNG[i][2] *= var;
      }
    }
  }
}

void FixEPHAtomic::heat_solve() {
  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // test stability
  double stability = 1.0;

  double loops = 2;
  double scaling = 1.0 / loops;
  double dt = scaling * dtv; // we will do extra steps if necessary

  for(size_t i = 0; i < loops; ++i) {
    { // add small portion of energy and redistribute temperatures
      for(size_t j = 0; j < nlocal; ++j) {
        if(mask[j] & groupbit) {
          int jtype = type[j];
          E_a_i[j] += dE_a_i[j] * scaling;
        }
      }

      state = FixState::EI;
      comm->forward_comm_fix(this);
    }

    { // solve diffusion a bit
      for(size_t j = 0; j < nlocal; ++j) {
        if(mask[j] & groupbit) {
          int jtype = type[j];
          int *klist = firstneigh[j];
          int knum = numneigh[j];

          if(!(rho_a_i[j] > 0)) continue;

          for(size_t k = 0; k != knum; ++k) {
            int kk = klist[k];
            kk &= NEIGHMASK;
            int ktype = type[kk];
          }
        }
      }
    }
  }
}

void FixEPHAtomic::post_force(int vflag) {
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *numneigh = list->numneigh;

  //zero all arrays
  std::fill_n(&(w_i[0][0]), 3 * nlocal, 0);
  std::fill_n(&(xi_i[0][0]), 3 * nlocal, 0);
  std::fill_n(&(f_EPH[0][0]), 3 * nlocal, 0);
  std::fill_n(&(f_RNG[0][0]), 3 * nlocal, 0);

  // generate random forces and distribute them
  if(eph_flag & Flag::RANDOM) {
    for(size_t i = 0; i < nlocal; ++i) {
      if(mask[i] & groupbit) {
        xi_i[i][0] = random->gaussian();
        xi_i[i][1] = random->gaussian();
        xi_i[i][2] = random->gaussian();
      }
    }

    state = FixState::XIX;
    comm->forward_comm_fix(this);
    state = FixState::XIY;
    comm->forward_comm_fix(this);
    state = FixState::XIZ;
    comm->forward_comm_fix(this);
  }

  // calculate the site densities, gradients (future) and beta(rho)
  calculate_environment();

  state = FixState::RHO;
  comm->forward_comm_fix(this);

  force_prl();

  // second loop over atoms if needed
  if((eph_flag & Flag::FRICTION) && !(eph_flag & Flag::NOFRICTION)) {
    for(int i = 0; i < nlocal; i++) {
      f[i][0] += f_EPH[i][0];
      f[i][1] += f_EPH[i][1];
      f[i][2] += f_EPH[i][2];
    }
  }

  if((eph_flag & Flag::RANDOM) && !(eph_flag & Flag::NORANDOM)) {
    for(int i = 0; i < nlocal; i++) {
      f[i][0] += f_RNG[i][0];
      f[i][1] += f_RNG[i][1];
      f[i][2] += f_RNG[i][2];
    }
  }
}

void FixEPHAtomic::reset_dt() {
  eta_factor = sqrt(2.0 * force->boltz / update->dt);

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

void FixEPHAtomic::grow_arrays(int ngrow) {
  //std::cout << "NGROW NLOCAL NGHOST NMAX\n";
  //std::cout << ngrow << ' ' <<
  //  atom->nlocal << ' ' << atom->nghost << ' ' << atom->nmax << '\n';
  n = ngrow;

  memory->grow(f_EPH, ngrow, 3,"eph_atomic:fEPH");
  memory->grow(f_RNG, ngrow, 3,"eph_atomic:fRNG");

  memory->grow(rho_i, ngrow, "eph:rho_i");

  memory->grow(w_i, ngrow, 3, "eph:w_i");
  memory->grow(xi_i, ngrow, 3, "eph:xi_i");

  memory->grow(rho_a_i, ngrow, "eph:rho_a_i");
  memory->grow(E_a_i, ngrow, "eph:E_a_i");
  memory->grow(dE_a_i, ngrow, "eph:dE_a_i");

  // per atom values
  // we need only nlocal elements here
  memory->grow(array, ngrow, size_peratom_cols, "eph:array");
  array_atom = array;
}

double FixEPHAtomic::compute_vector(int i) {
  if(i == 0)
    return Ee;
  else if(i == 1) {
    return Te;
  }

  return Ee;
}

/** TODO: There might be synchronisation issues here; maybe should add barrier for sync **/
int FixEPHAtomic::pack_forward_comm(int n, int *list, double *data, int pbc_flag, int *pbc) {
  int m;
  m = 0;
  switch(state) {
    case FixState::RHO:
      for(size_t i = 0; i < n; ++i) {
        data[m++] = rho_i[list[i]];
        data[m++] = rho_a_i[list[i]];
      }
      break;
    case FixState::XI:
      for(size_t i = 0; i < n; ++i) {
        data[m++] = xi_i[list[i]][0];
        data[m++] = xi_i[list[i]][1];
        data[m++] = xi_i[list[i]][2];
      }
      break;
    case FixState::XIX:
      for(size_t i = 0; i < n; ++i) {
        data[m++] = xi_i[list[i]][0];
      }
      break;
    case FixState::XIY:
      for(size_t i = 0; i < n; ++i) {
        data[m++] = xi_i[list[i]][1];
      }
      break;
    case FixState::XIZ:
      for(size_t i = 0; i < n; ++i) {
        data[m++] = xi_i[list[i]][2];
      }
      break;
    case FixState::WI:
      for(size_t i = 0; i < n; ++i) {
        data[m++] = w_i[list[i]][0];
        data[m++] = w_i[list[i]][1];
        data[m++] = w_i[list[i]][2];
      }
      break;
    case FixState::WX:
      for(size_t i = 0; i < n; ++i) {
        data[m++] = w_i[list[i]][0];
      }
      break;
    case FixState::WY:
      for(size_t i = 0; i < n; ++i) {
        data[m++] = w_i[list[i]][1];
      }
      break;
    case FixState::WZ:
      for(size_t i = 0; i < n; ++i) {
        data[m++] = w_i[list[i]][2];
      }
      break;
    case FixState::EI:
      for(size_t i = 0; i < n; ++i) {
        data[m++] = E_a_i[list[i]];
      }
      break;
    default:
      break;
  }

  return m;
}

void FixEPHAtomic::unpack_forward_comm(int n, int first, double *data) {
  int m, last;
  m = 0;
  last = first + n;

  switch(state) {
    case FixState::RHO:
      for(size_t i = first; i < last; ++i) {
        rho_i[i] = data[m++];
        rho_a_i[i] = data[m++];
      }
      break;
    case FixState::XI:
      for(size_t i = first; i < last; ++i) {
        xi_i[i][0] = data[m++];
        xi_i[i][1] = data[m++];
        xi_i[i][2] = data[m++];
      }
      break;
    case FixState::XIX:
      for(size_t i = first; i < last; ++i) {
        xi_i[i][0] = data[m++];
      }
      break;
    case FixState::XIY:
      for(size_t i = first; i < last; ++i) {
        xi_i[i][1] = data[m++];
      }
      break;
    case FixState::XIZ:
      for(size_t i = first; i < last; ++i) {
        xi_i[i][2] = data[m++];
      }
      break;
    case FixState::WI:
      for(size_t i = first; i < last; ++i) {
        w_i[i][0] = data[m++];
        w_i[i][1] = data[m++];
        w_i[i][2] = data[m++];
      }
      break;
    case FixState::WX:
      for(size_t i = first; i < last; ++i) {
        w_i[i][0] = data[m++];
      }
      break;
    case FixState::WY:
      for(size_t i = first; i < last; ++i) {
        w_i[i][1] = data[m++];
      }
      break;
    case FixState::WZ:
      for(size_t i = first; i < last; ++i) {
        w_i[i][2] = data[m++];
      }
      break;
    case FixState::EI:
      for(size_t i = first; i < last; ++i) {
        E_a_i[i] = data[m++];
      }
      break;
    default:
      break;
  }
}

/** TODO **/
double FixEPHAtomic::memory_usage() {
    double bytes = 0;

    return bytes;
}

/* save temperature state after run */
void FixEPHAtomic::post_run() {
  // save temperatures somehow (MPI maybe)
}

