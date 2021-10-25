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
#include "fix_eph_coloured_exp.h"
#include "eph_beta.h"
#include "eph_fdm.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/**
   * FixEPHMem arguments
   * arg[ 0] <- fix ID
   * arg[ 1] <- group
   * arg[ 2] <- name
   * arg[ 3] <- rng seed
   * arg[ 4] <- eph parameter; 0 disable all terms; 1 enable friction term; 2 enable random force; 4 enable fde;
   * arg[ 5] <- tau0; coloured noise parameter for exp // in time units
   * arg[ 6] <- electronic density; this might be changed with new fdm model // this might be used in the future
   * arg[ 7] <- electronic heat capacity
   * arg[ 8] <- electronic heat conduction
   * arg[ 9] <- initial electronic temperature TODO
   * arg[10] <- FDE grid x
   * arg[11] <- FDE grid y
   * arg[12] <- FDE grid z
   * arg[13] <- input file for initial temperatures
   * arg[14] <- frequency of output file writing
   * arg[15] <- output file for temperatures
   * arg[16] <- input file for eph model functions
   * arg[17] <- element name for type 0
   * arg[18] <- element name for type 1
   * ...
   **/

// constructor
FixEPHColouredExp::FixEPHColouredExp(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if (narg < 18) error->all(FLERR, "Illegal fix eph command: too few arguments");
  if (atom->natoms < 1) error->all(FLERR, "fix_eph: error no atoms in simulation");
  MPI_Comm_rank(world, &myID);
  MPI_Comm_size(world, &nrPS);

  state = FixState::NONE;

  vector_flag = 1; // fix is able to output a vector compute
  size_vector = 2; // 2 elements in the vector
  global_freq = 1; // frequency for vector data
  extvector = 1; // external vector allocated by this fix???
  nevery = 1; // call end_of_step every step
  peratom_flag = 1; // fix provides per atom values
  size_peratom_cols = 8; // per atom has 8 dimensions
  peratom_freq = 1; // per atom values are provided every step
  
  comm_forward = 3; // forward communication is needed
  comm->ghost_velocity = 1; // special: fix requires velocities for ghost atoms
  
  //~ maxexchange = 3;
  
  // initialise rng
  seed = atoi(arg[3]);
  random = new RanMars(lmp, seed + myID);

  // read model behaviour parameters
  eph_flag = strtol(arg[4], NULL, 0);

  // print enabled fix functionality
  if(myID == 0) {
    std::cout << '\n';
    std::cout << "Flag read: " << arg[4] << " -> " << eph_flag << '\n';
    if(eph_flag & Flag::FRICTION) std::cout << "Friction evaluation: ON\n";
    if(eph_flag & Flag::RANDOM) std::cout << "Random evaluation: ON\n";
    if(eph_flag & Flag::FDM) std::cout << "FDM grid solving: ON\n";
    if(eph_flag & Flag::NOINT) std::cout << "No integration: ON\n";
    if(eph_flag & Flag::NOFRICTION) std::cout << "No friction application: ON\n";
    if(eph_flag & Flag::NORANDOM) std::cout << "No random application: ON\n";
    std::cout << '\n';
  }

  // read model selection
  tau0 = atof(arg[5]);
  
  // electronic structure parameters
  double v_rho = atof(arg[6]);
  double v_Ce = atof(arg[7]);
  double v_kappa = atof(arg[8]);
  double v_Te = atof(arg[9]);
  int nx = atoi(arg[10]);
  int ny = atoi(arg[11]);
  int nz = atoi(arg[12]);

  /** initialise FDM grid **/
  // if filename is provided use that to initialise grid everything else is ignored
  if(strcmp("NULL" , arg[13]) == 0) {
    if(nx < 1 || ny < 1 || nz < 1) {
      error->all(FLERR, "FixEPH: non-positive grid values");
    }

    double x0 = domain->boxlo[0];
    double x1 = domain->boxhi[0];
    double y0 = domain->boxlo[1];
    double y1 = domain->boxhi[1];
    double z0 = domain->boxlo[2];
    double z1 = domain->boxhi[2];

    fdm = EPH_FDM(nx, ny, nz,
      x0, x1, y0, y1, z0, z1,
      v_Te, v_Ce, v_rho, v_kappa);

    // now the FDM should be defined
    strcpy(T_state, "T.restart");
  }
  else {
    fdm = EPH_FDM(arg[13]);

    sprintf(T_state, "%s.restart", arg[13]);
  }

  T_freq = atoi(arg[14]);
  if(T_freq > 0) { sprintf(T_out, "%s", arg[15]); }

  // set the communicator
  fdm.set_comm(world, myID, nrPS);
  fdm.set_dt(update->dt);

  // initialise beta(rho)
  types = atom->ntypes;

  if(types > (narg - 17)) {
    error->all(FLERR, "Fix eph: number of types larger than provided in fix");
  }

  type_map = new int[types]; // TODO: switch to vector

  beta = Beta(arg[16]);

  if(beta.get_n_elements() < 1) {
    error->all(FLERR, "Fix eph: no elements found in input file");
  }

  r_cutoff = beta.get_r_cutoff();
  r_cutoff_sq = beta.get_r_cutoff_sq();
  rho_cutoff = beta.get_rho_cutoff();

  // do element mapping
  for(size_t i = 0; i < types; ++i) {
    type_map[i] = std::numeric_limits<int>::max();

    for(size_t j = 0; j < beta.get_n_elements(); ++j)
      if((beta.get_element_name(j)).compare(arg[17+i]) == 0) type_map[i] = j;

    if(type_map[i] > types)
      error->all(FLERR, "Fix eph: elements not found in input file");
  }

  // set force prefactors
  eta_factor = sqrt(2.0 * force->boltz / update->dt);
  zeta_factor = 1.0 - exp(- update->dt / tau0);
  /** DEBUG **/
  //~ std::cout << "DEBUG: dt: " << update->dt << " zeta_factor: " << zeta_factor << '\n';
  
  /** integrator functionality **/
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  // allocate arrays for fix_eph
  f_EPH = nullptr;
  f_RNG = nullptr;

  w_i = nullptr;

  rho_i = nullptr;
  array = nullptr;

  xi_i = nullptr;
  zi_i = nullptr;
  
  zv_i = nullptr;

  T_e_i = nullptr;

  list = nullptr;

  // NO ARRAYS BEFORE THIS
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  // zero arrays, so they would not contain garbage
  size_t nlocal = atom->nlocal;
  size_t ntotal = atom->nghost + nlocal;

  std::fill_n(&(rho_i[0]), ntotal, 0);
  std::fill_n(&(xi_i[0][0]), 3 * ntotal, 0);
  std::fill_n(&(zi_i[0][0]), 3 * ntotal, 0);
  
  std::fill_n(&(zv_i[0][0]), 3 * ntotal, 0);
  std::fill_n(&(w_i[0][0]), 3 * ntotal, 0);

  std::fill_n(&(T_e_i[0]), ntotal, 0);

  std::fill_n(&(f_EPH[0][0]), 3 * ntotal, 0);
  std::fill_n(&(f_RNG[0][0]), 3 * ntotal, 0);

  std::fill_n(&(array[0][0]), size_peratom_cols * ntotal, 0);

  Ee = 0.0; // electronic energy is zero in the beginning
}

// destructor
FixEPHColouredExp::~FixEPHColouredExp() {
  delete random;
  delete[] type_map;

  atom->delete_callback(id, 0);

  memory->destroy(rho_i);

  memory->destroy(array);

  memory->destroy(f_EPH);
  memory->destroy(f_RNG);

  memory->destroy(xi_i);
  memory->destroy(zi_i);
  
  memory->destroy(zv_i);
  memory->destroy(w_i);

  memory->destroy(T_e_i);
}

void FixEPHColouredExp::init() {
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix eph with 2d simulation");
  if (domain->nonperiodic != 0)
    error->all(FLERR,"Cannot use nonperiodic boundares with fix eph");
  if (domain->triclinic)
    error->all(FLERR,"Cannot use fix eph with triclinic box");

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

void FixEPHColouredExp::init_list(int id, NeighList *ptr) {
  this->list = ptr;
}

int FixEPHColouredExp::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  /* integrator functionality */
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;

  return mask;
}

/* integrator functionality */
void FixEPHColouredExp::initial_integrate(int) {
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

void FixEPHColouredExp::final_integrate() {
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

void FixEPHColouredExp::end_of_step() {
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

        fdm.insert_energy(x[i][0], x[i][1], x[i][2], dE_i);
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

        fdm.insert_energy(x[i][0], x[i][1], x[i][2], dE_i);
        E_local += dE_i;
      }
    }
  }

  if(eph_flag & Flag::FDM) {
    fdm.solve();
  }

  // save heatmap
  if(myID == 0 && T_freq > 0 && (update->ntimestep % T_freq) == 0) { // TODO: implement a counter instead
    fdm.save_temperature(T_out, update->ntimestep / T_freq);
  }

  // this is for checking energy conservation
  MPI_Allreduce(MPI_IN_PLACE, &E_local, 1, MPI_DOUBLE, MPI_SUM, world);

  Ee += E_local;

  for(size_t i = 0; i < nlocal; ++i) {
    if(mask[i] & groupbit) {
      int itype = type[i];
      array[i][ 0] = rho_i[i];
      array[i][ 1] = beta.get_beta(type_map[itype - 1], rho_i[i]);
      array[i][ 2] = f_EPH[i][0];
      array[i][ 3] = f_EPH[i][1];
      array[i][ 4] = f_EPH[i][2];
      array[i][ 5] = f_RNG[i][0];
      array[i][ 6] = f_RNG[i][1];
      array[i][ 7] = f_RNG[i][2];
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
    }
  }
}

void FixEPHColouredExp::calculate_environment() {
  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // loop over atoms and their neighbours and calculate rho and beta(rho)
  for(size_t i = 0; i != nlocal; ++i) {
    rho_i[i] = 0;

    // check if current atom belongs to fix group and if an atom is local
    if(mask[i] & groupbit) {
      int itype = type[i];
      int *jlist = firstneigh[i];
      int jnum = numneigh[i];

      for(size_t j = 0; j != jnum; ++j) {
        int jj = jlist[j];
        jj &= NEIGHMASK;

        int jtype = type[jj];
        double r_sq = get_distance_sq(x[jj], x[i]);

        if(r_sq < r_cutoff_sq) {
          rho_i[i] += beta.get_rho_r_sq(type_map[jtype-1], r_sq);
        }
      }
    }
  }
}

void FixEPHColouredExp::force_prl() {
  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // create friction forces
  if(eph_flag & Flag::FRICTION) {
    // coloured dissipation
    for(size_t i = 0; i < nlocal; ++i) {
      zv_i[i][0] = zv_i[i][0] * (1. - zeta_factor) + zeta_factor * v[i][0];
      zv_i[i][1] = zv_i[i][1] * (1. - zeta_factor) + zeta_factor * v[i][1];
      zv_i[i][2] = zv_i[i][2] * (1. - zeta_factor) + zeta_factor * v[i][2];
    }
    
    state = FixState::ZV;
    comm->forward_comm_fix(this);
    
    // w_i = W_ij^T v_j
    for(size_t i = 0; i < nlocal; ++i) {
      if(mask[i] & groupbit) {
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];

        if(!(rho_i[i] > 0)) continue;

        double alpha_i = beta.get_alpha(type_map[itype - 1], rho_i[i]);

        for(size_t j = 0; j != jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];

          // calculate the e_ij vector TODO: change these
          double e_ij[3];
          double e_r_sq = get_difference_sq(x[jj], x[i], e_ij);

          // first sum
          if(e_r_sq >= r_cutoff_sq) continue;

          double v_rho_ji = beta.get_rho_r_sq(type_map[jtype - 1], e_r_sq);
          double prescaler = alpha_i * v_rho_ji / (rho_i[i] * e_r_sq);

          double e_v_v1 = get_scalar(e_ij, zv_i[i]);
          double var1 = prescaler * e_v_v1;

          double e_v_v2 = get_scalar(e_ij, zv_i[jj]);
          double var2 = prescaler * e_v_v2;

          double dvar = var1 - var2;
          w_i[i][0] += dvar * e_ij[0];
          w_i[i][1] += dvar * e_ij[1];
          w_i[i][2] += dvar * e_ij[2];
        }
      }
    }

    state = FixState::WI;
    comm->forward_comm_fix(this);
    
    // now calculate the forces
    // f_i = W_ij w_j
    for(size_t i = 0; i != nlocal; ++i) {
      if(mask[i] & groupbit) {
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];

        if(not(rho_i[i] > 0.)) { continue; }

        double alpha_i = beta.get_alpha(type_map[itype - 1], rho_i[i]);

        for(size_t j = 0; j != jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];

          // calculate the e_ij vector
          double e_ij[3];
          double e_r_sq = get_difference_sq(x[jj], x[i], e_ij);

          if(e_r_sq >= r_cutoff_sq or not(rho_i[jj] > 0)) { continue; }

          double alpha_j = beta.get_alpha(type_map[jtype - 1], rho_i[jj]);

          double v_rho_ji = beta.get_rho_r_sq(type_map[jtype - 1], e_r_sq);
          double e_v_v1 = get_scalar(e_ij, w_i[i]);
          double var1 = alpha_i * v_rho_ji * e_v_v1 / (rho_i[i] * e_r_sq);

          double v_rho_ij = beta.get_rho_r_sq(type_map[itype - 1], e_r_sq);
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

        double alpha_i = beta.get_alpha(type_map[itype - 1], rho_i[i]);

        for(size_t j = 0; j != jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];

          // calculate the e_ij vector
          double e_ij[3];
          double e_r_sq = get_difference_sq(x[jj], x[i], e_ij);

          if((e_r_sq >= r_cutoff_sq) || !(rho_i[jj] > 0)) continue;

          double alpha_j = beta.get_alpha(type_map[jtype - 1], rho_i[jj]);

          double v_rho_ji = beta.get_rho_r_sq(type_map[jtype - 1], e_r_sq);
          double e_v_zi1 = get_scalar(e_ij, zi_i[i]);
          double var1 = alpha_i * v_rho_ji * e_v_zi1 / (rho_i[i] * e_r_sq);

          double v_rho_ij = beta.get_rho_r_sq(type_map[itype - 1], e_r_sq);
          double e_v_zi2 = get_scalar(e_ij, zi_i[jj]);
          double var2 = alpha_j * v_rho_ij * e_v_zi2 / (rho_i[jj] * e_r_sq);

          double dvar = var1 - var2;
          f_RNG[i][0] += dvar * e_ij[0];
          f_RNG[i][1] += dvar * e_ij[1];
          f_RNG[i][2] += dvar * e_ij[2];
        }

        double v_Te = fdm.get_T(x[i][0], x[i][1], x[i][2]);
        double var = eta_factor * sqrt(v_Te);
        f_RNG[i][0] *= var;
        f_RNG[i][1] *= var;
        f_RNG[i][2] *= var;
      }
    }
  }
}

void FixEPHColouredExp::post_force(int vflag) {
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
    
    // colourize noise
    for(size_t i = 0; i < nlocal; ++i) {
      zi_i[i][0] = zi_i[i][0] * (1. - zeta_factor) + zeta_factor * xi_i[i][0];
      zi_i[i][1] = zi_i[i][1] * (1. - zeta_factor) + zeta_factor * xi_i[i][1];
      zi_i[i][2] = zi_i[i][2] * (1. - zeta_factor) + zeta_factor * xi_i[i][2];
    }
    
    state = FixState::ZI;
    comm->forward_comm_fix(this);
  }
  
  // calculate the site densities, gradients (future) and beta(rho)
  calculate_environment();

  state = FixState::RHO;
  comm->forward_comm_fix(this);
  
  force_prl(); // calculate the dissipation and friction forces

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

void FixEPHColouredExp::reset_dt() {
  eta_factor = sqrt(2.0 * force->boltz / update->dt);
  zeta_factor = 1.0 - exp(- update->dt / tau0);
  
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  fdm.set_dt(update->dt);
}

void FixEPHColouredExp::grow_arrays(int ngrow) {
  n = ngrow;

  memory->grow(f_EPH, ngrow, 3,"EPH:fEPH");
  memory->grow(f_RNG, ngrow, 3,"EPH:fRNG");

  memory->grow(rho_i, ngrow, "eph:rho_i");

  memory->grow(w_i, ngrow, 3, "eph:w_i");
  memory->grow(xi_i, ngrow, 3, "eph:xi_i");
  memory->grow(zi_i, ngrow, 3, "eph:zi_i");
  
  memory->grow(zv_i, ngrow, 3, "eph:zv_i");

  memory->grow(T_e_i, ngrow, "eph:T_e_i");

  // per atom values
  // we need only nlocal elements here
  memory->grow(array, ngrow, size_peratom_cols, "eph:array");
  array_atom = array;
}

double FixEPHColouredExp::compute_vector(int i) {
  if(i == 0)
    return Ee;
  else if(i == 1) {
    return fdm.get_T_total();
  }

  return Ee;
}

/** TODO: There might be synchronisation issues here; maybe should add barrier for sync **/
int FixEPHColouredExp::pack_forward_comm(int n, int *list, double *data, int pbc_flag, int *pbc) {
  int m;
  m = 0;
  switch(state) {
    case FixState::RHO:
      for(size_t i = 0; i < n; ++i) {
        data[m++] = rho_i[list[i]];
      }
      break;
    case FixState::ZI:
      for(size_t i = 0; i < n; ++i) {
        data[m++] = zi_i[list[i]][0];
        data[m++] = zi_i[list[i]][1];
        data[m++] = zi_i[list[i]][2];
      }
      break;
    case FixState::ZV:
      for(size_t i = 0; i < n; ++i) {
        data[m++] = zv_i[list[i]][0];
        data[m++] = zv_i[list[i]][1];
        data[m++] = zv_i[list[i]][2];
      }
      break;
    case FixState::WI:
      for(size_t i = 0; i < n; ++i) {
        data[m++] = w_i[list[i]][0];
        data[m++] = w_i[list[i]][1];
        data[m++] = w_i[list[i]][2];
      }
      break;
    default:
      break;
  }

  return m;
}

void FixEPHColouredExp::unpack_forward_comm(int n, int first, double *data) {
  int m, last;
  m = 0;
  last = first + n;

  switch(state) {
    case FixState::RHO:
      for(size_t i = first; i < last; ++i) {
        rho_i[i] = data[m++];
      }
      break;
    case FixState::ZI:
      for(size_t i = first; i < last; ++i) {
        zi_i[i][0] = data[m++];
        zi_i[i][1] = data[m++];
        zi_i[i][2] = data[m++];
      }
      break;
    case FixState::ZV:
      for(size_t i = first; i < last; ++i) {
        zv_i[i][0] = data[m++];
        zv_i[i][1] = data[m++];
        zv_i[i][2] = data[m++];
      }
      break;
    case FixState::WI:
      for(size_t i = first; i < last; ++i) {
        w_i[i][0] = data[m++];
        w_i[i][1] = data[m++];
        w_i[i][2] = data[m++];
      }
      break;
    default:
      break;
  }
}

/** TODO **/
double FixEPHColouredExp::memory_usage() { return 0; }

/* save temperature state after run */
void FixEPHColouredExp::post_run() {
  if(myID == 0) fdm.save_state(T_state);
}

int FixEPHColouredExp::pack_exchange(int i, double *buf) {
  size_t m = 0;
  buf[m++] = zi_i[i][0];
  buf[m++] = zi_i[i][1];
  buf[m++] = zi_i[i][2];
  
  buf[m++] = zv_i[i][0];
  buf[m++] = zv_i[i][1];
  buf[m++] = zv_i[i][2];
  return m;
}

int FixEPHColouredExp::unpack_exchange(int nlocal, double *buf) {
  size_t m = 0;
  zi_i[nlocal][0] = buf[m++];
  zi_i[nlocal][1] = buf[m++];
  zi_i[nlocal][2] = buf[m++];
  
  zv_i[nlocal][0] = buf[m++];
  zv_i[nlocal][1] = buf[m++];
  zv_i[nlocal][2] = buf[m++];
  return m;
}

