/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

//#define DEBUG_EPH

// external headers
#include <iostream>

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

#include <cstring> // TODO: remove
#include <string>
#include <limits>
#include <algorithm>

// internal headers
#include "fix_eph.h"
#include "eph_beta.h"
#include "eph_fdm.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/**
   * FixEPH arguments
   * arg[ 0] <- fix ID
   * arg[ 1] <- group
   * arg[ 2] <- name
   * arg[ 3] <- rng seed
   * arg[ 4] <- eph parameter; 0 disable all terms; 1 enable friction term; 2 enable random force; 4 enable fde;
   * arg[ 5] <- eph model selection; 1 standard langevin with beta(rho); 2 PRB version; 3 new model with CM only; 4 full new model ; 9 Testing
   * arg[ 6] <- electronic density; this might be changed with new fdm model // this is in principle unnessesary
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
FixEPH::FixEPH(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if (narg < 18) error->all(FLERR, "Illegal fix eph command: too few arguments");
  if (atom->natoms < 1) error->all(FLERR, "fix_eph: error no atoms in simulation");
  MPI_Comm_rank(world, &myID);
  MPI_Comm_size(world, &nrPS);
  state = FixState::NONE;
  
  vector_flag = 1; // fix is able to output a vector compute
  size_vector = 2; // 2 elements in the vector
  global_freq = 1; // frequency for vector data
  extvector = 1; // 
  nevery = 1; // call end_of_step every step
  peratom_flag = 1;
  size_peratom_cols = 8;
  peratom_freq = 1;
  
  /** comm test **/
  comm_forward = 1;
  //comm_reverse = 1;
  
  // special needed for velocity things to work
  comm->ghost_velocity = 1;
  
  // initialise rng
  seed = atoi(arg[3]);
  random = new RanMars(lmp, seed + myID);
  
  // read model behaviour parameters
  eph_flag = atoi(arg[4]);
  eph_model = atoi(arg[5]);
  
  // TODO: magic parameters for passing values TEMPORARY
  v_alpha = 1.0;
  v_struc = 1.0;
  
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
  if(strcmp("NULL" , arg[13]) == 0) { // this might not be the best test
    if(nx < 1 || ny < 1 || nz < 1) { // TODO: negative values could be used for sth
      error->all(FLERR, "FixEPH: non-positive grid values");
    }
    fdm = new EPH_FDM(nx, ny, nz);
    double x0 = domain->boxlo[0];
    double x1 = domain->boxhi[0];
    double y0 = domain->boxlo[1];
    double y1 = domain->boxhi[1];
    double z0 = domain->boxlo[2];
    double z1 = domain->boxhi[2];
    fdm->setBox(x0, x1, y0, y1, z0, z1);
    fdm->setConstants(v_Te, v_Ce, v_rho, v_kappa);
    // now the FDM should be defined
    strcpy(Tstate, "T.restart");
  }
  else {
    fdm = new EPH_FDM(arg[13]);
    /* TODO: add error check here
    if(fdm == NULL) { // this does not do the correct thing
      error->all(FLERR, "FixEPH: loading FDM parameters from file failed.");
    }*/
    sprintf(Tstate, "%s.restart", arg[13]);
  }
  
  Tfreq = atoi(arg[14]);
  if(Tfreq > 0) {
    sprintf(Tout, "%s", arg[15]);
  }
  
  // set the communicator
  fdm->setComm(world, myID, nrPS);
  // oops
  fdm->setDt(update->dt);
  
  // initialise beta(rho)
  types = atom->ntypes;
  
  if(types > (narg - 17)) {
    error->all(FLERR, "Fix eph: number of types larger than provided in fix");
  }
  
  typeMap = new unsigned int[types];
  // TODO: add error check here
  beta = new EPH_Beta(arg[16]);
  
  if(beta->getElementsNumber() < 1) {
    error->all(FLERR, "Fix eph: no elements found in input file");
  }
  
  rcutoff = beta->getCutoff();
  
  // do element mapping
  for(unsigned int i = 0; i < types; ++i) {
    typeMap[i] = std::numeric_limits<unsigned int>::max();
    for(unsigned int j = 0; j < beta->getElementsNumber(); ++j) {
      if((beta->getName(j)).compare(arg[17+i]) == 0) typeMap[i] = j;
    }
    if(typeMap[i] > types) {
      error->all(FLERR, "Fix eph: elements not found in input file");
    }
  }
  
  // set force prefactors
  beta_factor = 1.0;
  eta_factor = sqrt(2.0 * force->boltz / update->dt);
  
  /** integrator functionality **/
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  
  // beta_factor = 1.0 / force->ftm2v;
  
  // uniform distribution
  // eta_factor = sqrt(24.0 * force->boltz / update->dt / force->mvv2e) / force->ftm2v;
  
  // gaussian distribution
  // eta_factor = sqrt(2.0 * force->boltz / update->dt / force->mvv2e) / force->ftm2v;
  
  // allocate arrays for fix_eph
  f_EPH = NULL;
  f_RNG = NULL;
  
  w_i = NULL;

  beta_i = NULL;
  rho_i = NULL;
  array = NULL;
  
  xi_i = NULL;
  
  rho_neigh = 512;
  rho_ij = NULL;
  rho_ji = NULL;
  
  grad_rho_i = NULL;
  
  list = NULL;
  
  // NO ARRAYS BEFORE THIS
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  
  // zero arrays, so they would not contain garbage
  int nlocal = atom->nlocal;
  //t_nlocal = nlocal;
  std::fill_n(&(rho_i[0]), nlocal, 0.0D);
  std::fill_n(&(beta_i[0]), nlocal, 0.0D);
  std::fill_n(&(xi_i[0][0]), 3 * nlocal, 0.0D);
  std::fill_n(&(w_i[0][0]), 3 * nlocal, 0.0D);
  std::fill_n(&(f_EPH[0][0]), 3 * nlocal, 0.0D);
  std::fill_n(&(f_RNG[0][0]), 3 * nlocal, 0.0D);
  std::fill_n(&(grad_rho_i[0][0]), 3 * nlocal, 0.0D);
  
  for(int i = 0; i < nlocal; ++i) {
    for(int j = 0; j < size_peratom_cols; ++j) {
      array[i][j] = 0.0;
    }
  }
  
  Ee = 0.0; // electronic energy is zero in the beginning
  
  // we ignore all other input parameters at the moment
}

// destructor
FixEPH::~FixEPH() {
  delete random;
  delete[] typeMap;
  delete beta;
  delete fdm;
  
  atom->delete_callback(id, 0);
  
  memory->destroy(beta_i);
  memory->destroy(rho_i);
  
  memory->destroy(array);
  
  memory->destroy(f_EPH);
  memory->destroy(f_RNG);
  
  memory->destroy(rho_ij);
  memory->destroy(rho_ji);
  
  memory->destroy(grad_rho_i);
  
  memory->destroy(xi_i);
  memory->destroy(w_i);
}

void FixEPH::init() {
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix eph with 2d simulation");
  if (domain->nonperiodic != 0)
    error->all(FLERR,"Cannot use nonperiodic boundares with fix eph");
  if (domain->triclinic)
    error->all(FLERR,"Cannot use fix eph with triclinic box");
    
  /* copy paste from vcsgc */
  /** we are a fix and we need full neighbour list **/
  int irequest = neighbor->request((void*)this, this->instance_me);
  neighbor->requests[irequest]->pair=0;
  neighbor->requests[irequest]->fix=1;
  neighbor->requests[irequest]->half=0;
  neighbor->requests[irequest]->full=1;
  
  reset_dt();
}

void FixEPH::init_list(int id, NeighList *ptr) {
  this->list = ptr;
}

int FixEPH::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  /** integrator functionality **/
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  
  return mask;
}

/** integrator functionality **/
void FixEPH::initial_integrate(int) {
  //if(myID == 0) std::cout << "DEBUG initial_integrate()" << std::endl;
  if(eph_flag & Flag::NOINT) return;

  double dtfm;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
  }
}

void FixEPH::final_integrate() {
  //sif(myID == 0) std::cout << "DEBUG final_integrate()" << std::endl;
  if(eph_flag & Flag::NOINT) return;
  
  double dtfm;

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
    }
  }
}

void FixEPH::end_of_step() {
  //if(myID == 0) std::cout << "DEBUG end_of_step()" << std::endl;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *tag = atom->tag;
  
  double E_local = 0.0;
  
  // calculate the energy transferred to electronic system
  // this is a potential source of errors due to using velocity verlet for integration
  // friction force depends on the velocities and therefore acceleration at 
  //   next timestep depends on the velocity at next time step
  // this leads to errors of the order of dt^2
  if(eph_flag & Flag::FRICTION) {
    for(unsigned int i = 0; i < nlocal; ++i) {
      if(mask[i] & groupbit) {
        double dE_i = 0.0;
        dE_i -= f_EPH[i][0]*v[i][0]*update->dt;
        dE_i -= f_EPH[i][1]*v[i][1]*update->dt;
        dE_i -= f_EPH[i][2]*v[i][2]*update->dt;
        
        fdm->insertEnergy(x[i][0], x[i][1], x[i][2], dE_i);
        E_local += dE_i;
      }
    }
  }
  
  if(eph_flag & Flag::RANDOM) {
    for(unsigned int i = 0; i < nlocal; ++i) {
      if(mask[i] & groupbit) {
        double dE_i = 0.0;
        dE_i -= f_RNG[i][0]*v[i][0]*update->dt;
        dE_i -= f_RNG[i][1]*v[i][1]*update->dt;
        dE_i -= f_RNG[i][2]*v[i][2]*update->dt;
        
        fdm->insertEnergy(x[i][0], x[i][1], x[i][2], dE_i);
        E_local += dE_i;
      }
    }
  }
  
  if(eph_flag & Flag::FDM) {
    fdm->solve();
  }
  
  // save heatmap
  if(myID == 0 && Tfreq > 0 && (update->ntimestep % Tfreq) == 0) {
    fdm->saveTemperature(Tout, update->ntimestep / Tfreq);
  }
  
  // this is for checking energy conservation
  MPI_Allreduce(MPI_IN_PLACE, &E_local, 1, MPI_DOUBLE, MPI_SUM, world);
  //if(v_alpha > 0.0)
  //  v_Te += E_local / v_Ce / v_rho; // v_Ce * v_rho [K / eV]
  
  Ee += E_local;
  
  // store per atom values into array
  /*if(nlocal != t_nlocal) {
    std::cout << "NLOCAL CHANGED " << t_nlocal << " -> " << nlocal << std::endl;
    t_nlocal = nlocal;
  }*/
  for(unsigned int i = 0; i < nlocal; ++i) {
    if(mask[i] & groupbit) {
      array[i][ 0] = rho_i[i];
      array[i][ 1] = beta_i[i];
      array[i][ 2] = w_i[i][0];
      array[i][ 3] = w_i[i][1];
      array[i][ 4] = w_i[i][2];
      array[i][ 5] = xi_i[i][0];
      array[i][ 6] = xi_i[i][1];
      array[i][ 7] = xi_i[i][2];
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

void FixEPH::calculate_environment() {
  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  
  // loop over atoms and their neighbours and calculate rho and beta(rho)
  /* TODO: inum vs nlocal */
  for(int i = 0; i < nlocal; ++i) {
    
    // check if current atom belongs to fix group and if an atom is local
    if(mask[i] & groupbit) {
      int itype = type[i];
      int *jlist = firstneigh[i];
      int jnum = numneigh[i];
      
      for(int j = 0; j < jnum; ++j) {
        int jj = jlist[j];
        jj &= NEIGHMASK;
        
        int jtype = type[jj];
        
        double delx = x[jj][0] - x[i][0];
        double dely = x[jj][1] - x[i][1];
        double delz = x[jj][2] - x[i][2];
        
        double r = sqrt(delx*delx + dely*dely + delz*delz);
        
        if(r < rcutoff) {          
          double v_rho_ji = beta->getRho(typeMap[jtype-1], r);
          double v_rho_ij = beta->getRho(typeMap[itype-1], r);
          
          rho_ij[i][j] = v_rho_ij;
          rho_ji[i][j] = v_rho_ji;
          
          rho_i[i] += v_rho_ji;
          
          double v_grad_rho = -beta->getDRho(typeMap[jtype-1], r);
          
          grad_rho_i[i][0] += v_grad_rho * delx / r;
          grad_rho_i[i][1] += v_grad_rho * dely / r;
          grad_rho_i[i][2] += v_grad_rho * delz / r;
        }
      }
      
      beta_i[i] = beta->getBeta(typeMap[itype-1], rho_i[i]);
    }
  }
  
  state = FixState::RHO;
  comm->forward_comm_fix(this);
  
  state = FixState::BETA;
  comm->forward_comm_fix(this);
}

void FixEPH::force_ttm() {
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  
  // create friction forces
  if(eph_flag & Flag::FRICTION) {
    for(int i = 0; i < nlocal; ++i) {
      if(mask[i] & groupbit) {
        double var = -beta_factor * beta_i[i];
        f_EPH[i][0] = var * v[i][0];
        f_EPH[i][1] = var * v[i][1];
        f_EPH[i][2] = var * v[i][2];
      }
    }
  }
  
  // create random forces
  if(eph_flag & Flag::RANDOM) {
    for(int i = 0; i < nlocal; i++) {
      if(mask[i] & groupbit) {
        double v_Te = fdm->getT(x[i][0], x[i][1], x[i][2]);
        double var = eta_factor * sqrt(v_Te * beta_i[i]);
        f_RNG[i][0] = var * xi_i[i][0];
        f_RNG[i][1] = var * xi_i[i][1];
        f_RNG[i][2] = var * xi_i[i][2];
      }
    }
  }  
}

void FixEPH::force_prb() {
  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  
  // create friction forces
  if(eph_flag & Flag::FRICTION) {
    for(int i = 0; i < nlocal; ++i) {
      if(mask[i] & groupbit) {
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];
        
        double var = beta_i[i];
        f_EPH[i][0] = var * v[i][0];
        f_EPH[i][1] = var * v[i][1];
        f_EPH[i][2] = var * v[i][2];
        
        for(int j = 0; j < jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];
          
          if(rho_ji[i][j] > 0.0D && rho_i[i] > 0.0D) {
            double v_rho_ji = rho_ji[i][j];
            double var = beta_i[i] * v_rho_ji / rho_i[i];
            f_EPH[i][0] -= var * v[jj][0];
            f_EPH[i][1] -= var * v[jj][1];
            f_EPH[i][2] -= var * v[jj][2];
          }
        }
        
        // unit prefactor
        f_EPH[i][0] *= -beta_factor;
        f_EPH[i][1] *= -beta_factor;
        f_EPH[i][2] *= -beta_factor;
      }
    }
  }
  
  // create random forces
  if(eph_flag & Flag::RANDOM) {
    for(int i = 0; i < nlocal; i++) {
      if(mask[i] & groupbit) {
        double v_Te = fdm->getT(x[i][0], x[i][1], x[i][2]);
        double var = eta_factor * sqrt(v_Te * beta_i[i]);
        f_RNG[i][0] = var * xi_i[i][0];
        f_RNG[i][1] = var * xi_i[i][1];
        f_RNG[i][2] = var * xi_i[i][2];
      }
    }
  }
}

void FixEPH::force_prbmod() {
  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  
  // create friction forces
  if(eph_flag & Flag::FRICTION) {
    for(int i = 0; i < nlocal; ++i) {
      if(mask[i] & groupbit) {
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];
        
        double alpha_i = sqrt(beta_i[i]);
        w_i[i][0] = alpha_i * v[i][0];
        w_i[i][1] = alpha_i * v[i][1];
        w_i[i][2] = alpha_i * v[i][2];
        
        for(int j = 0; j < jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];
          
          if (rho_ji[i][j] > 0.0D && rho_i[i] > 0.0D) {
            double v_rho_ji = rho_ji[i][j];
            double var = alpha_i * v_rho_ji / rho_i[i];
            w_i[i][0] -= var * v[jj][0];
            w_i[i][1] -= var * v[jj][1];
            w_i[i][2] -= var * v[jj][2];
          }
        }
      }
    }
    
    state = FixState::WX;
    comm->forward_comm_fix(this);
    state = FixState::WY;
    comm->forward_comm_fix(this);
    state = FixState::WZ;
    comm->forward_comm_fix(this);
    //MPI_Allreduce(MPI_IN_PLACE, &(w_i[0][0]), 3 * atom->natoms, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    // now calculate the forces
    for(int i = 0; i < nlocal; ++i) {
      if(mask[i] & groupbit) {
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];
        
        double alpha_i = sqrt(beta_i[i]);
        f_EPH[i][0] = alpha_i * w_i[i][0];
        f_EPH[i][1] = alpha_i * w_i[i][1];
        f_EPH[i][2] = alpha_i * w_i[i][2];
        
        for(int j = 0; j < jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];
          
          if(rho_ij[i][j] > 0.0D &&rho_i[jj] > 0.0D) {
            double alpha_j = sqrt(beta_i[jj]);
            double v_rho_ij = rho_ij[i][j];
            double var = alpha_j * v_rho_ij / rho_i[jj];
            f_EPH[i][0] -= var * w_i[jj][0];
            f_EPH[i][1] -= var * w_i[jj][1];
            f_EPH[i][2] -= var * w_i[jj][2];
          }
        }
        
        // unit prefactor
        f_EPH[i][0] *= -beta_factor;
        f_EPH[i][1] *= -beta_factor;
        f_EPH[i][2] *= -beta_factor;
      }
    }
  }
  
  // create random forces
  if(eph_flag & Flag::RANDOM) {
    for(int i = 0; i < nlocal; i++) {
      if(mask[i] & groupbit) {
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];
        
        double var = sqrt(beta_i[i]);
        
        f_RNG[i][0] = var * xi_i[i][0];
        f_RNG[i][1] = var * xi_i[i][1];
        f_RNG[i][2] = var * xi_i[i][2];
        
        for(int j = 0; j < jnum; j++) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          
          if(rho_ij[i][j] > 0.0D && rho_i[jj] > 0.0) {          
            double v_rho_ij = rho_ij[i][j];
            var = sqrt(beta_i[jj]);
            var *= v_rho_ij / rho_i[jj];
            
            f_RNG[i][0] -= var * xi_i[jj][0];
            f_RNG[i][1] -= var * xi_i[jj][1];
            f_RNG[i][2] -= var * xi_i[jj][2];
          }
        }
        
        double v_Te = fdm->getT(x[i][0], x[i][1], x[i][2]);
        var = eta_factor * sqrt(v_Te);
        f_RNG[i][0] *= var;
        f_RNG[i][1] *= var;
        f_RNG[i][2] *= var;
      }
    }
  }
}
  
void FixEPH::force_eta() {
  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  
  // create friction forces
  if(eph_flag & Flag::FRICTION) {
    // w_i = W_ij^T v_j
    for(int i = 0; i < nlocal; ++i) {
      if(mask[i] & groupbit) {
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];
        
        double alpha_i = sqrt(beta_i[i]);
        
        for(int j = 0; j < jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];
          
          // calculate the e_ij vector
          double e_ij_x = x[jj][0] - x[i][0];
          double e_ij_y = x[jj][1] - x[i][1];
          double e_ij_z = x[jj][2] - x[i][2];
          
          double e_r_2 = e_ij_x * e_ij_x + e_ij_y * e_ij_y + e_ij_z * e_ij_z;
          
          double alpha_j = sqrt(beta_i[jj]);
          
          // first sum
          if (rho_ji[i][j] > 0.0D && rho_i[i] > 0.0D && e_r_2 > 0.0D) {
            double v_rho_ji = rho_ji[i][j];
            double e_v_v = e_ij_x * v[i][0] + 
                          e_ij_y * v[i][1] + 
                          e_ij_z * v[i][2];
            
            double var = alpha_i * v_rho_ji / rho_i[i] * e_v_v / e_r_2;
            
            w_i[i][0] += var * e_ij_x;
            w_i[i][1] += var * e_ij_y;
            w_i[i][2] += var * e_ij_z;
          }
          
          // second sum
          if (rho_ji[i][j] > 0.0D && rho_i[i] > 0.0D && e_r_2 > 0.0D) {
            double v_rho_ji = rho_ji[i][j];
            double e_v_v = e_ij_x * v[jj][0] + 
                          e_ij_y * v[jj][1] + 
                          e_ij_z * v[jj][2];
            
            double var = alpha_i * v_rho_ji / rho_i[i] * e_v_v / e_r_2;
            
            w_i[i][0] -= var * e_ij_x;
            w_i[i][1] -= var * e_ij_y;
            w_i[i][2] -= var * e_ij_z;
          }
        }
      }
    }
    
    state = FixState::WX;
    comm->forward_comm_fix(this);
    state = FixState::WY;
    comm->forward_comm_fix(this);
    state = FixState::WZ;
    comm->forward_comm_fix(this);
    
    //MPI_Allreduce(MPI_IN_PLACE, &(w_i[0][0]), 3 * atom->natoms, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    // now calculate the forces
    // f_i = W_ij w_j
    for(int i = 0; i < nlocal; ++i) {
      if(mask[i] & groupbit) {
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];
        
        double alpha_i = sqrt(beta_i[i]);
        
        for(int j = 0; j < jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];
          
          // calculate the e_ij vector
          double e_ij_x = x[jj][0] - x[i][0];
          double e_ij_y = x[jj][1] - x[i][1];
          double e_ij_z = x[jj][2] - x[i][2];
          
          double e_r_2 = e_ij_x * e_ij_x + e_ij_y * e_ij_y + e_ij_z * e_ij_z;
          
          double alpha_j = sqrt(beta_i[jj]);
          
          // first sum
          if (rho_ji[i][j] > 0.0D && rho_i[i] > 0.0D && e_r_2 > 0.0D) {
            double v_rho_ji = rho_ji[i][j];
            double e_v_v = e_ij_x * w_i[i][0] + 
                          e_ij_y * w_i[i][1] + 
                          e_ij_z * w_i[i][2];
            
            double var = alpha_i * v_rho_ji / rho_i[i] * e_v_v / e_r_2;
            
            f_EPH[i][0] += var * e_ij_x;
            f_EPH[i][1] += var * e_ij_y;
            f_EPH[i][2] += var * e_ij_z;
          }
          
          // second sum
          if (rho_ij[i][j] > 0.0D && rho_i[jj] > 0.0D && e_r_2 > 0.0D) {
            double v_rho_ij = rho_ij[i][j];
            double e_v_v = e_ij_x * w_i[jj][0] + 
                          e_ij_y * w_i[jj][1] + 
                          e_ij_z * w_i[jj][2];
            
            double var = alpha_j * v_rho_ij / rho_i[jj] * e_v_v / e_r_2;
            
            f_EPH[i][0] -= var * e_ij_x;
            f_EPH[i][1] -= var * e_ij_y;
            f_EPH[i][2] -= var * e_ij_z;
          }
        }
        
        // unit prefactor
        f_EPH[i][0] *= -beta_factor * v_struc * v_struc;
        f_EPH[i][1] *= -beta_factor * v_struc * v_struc;
        f_EPH[i][2] *= -beta_factor * v_struc * v_struc;
      }
    }
  }
  
  // create random forces
  if(eph_flag & Flag::RANDOM) {
    for(int i = 0; i < nlocal; i++) {
      if(mask[i] & groupbit) {
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];
        
        double alpha_i = sqrt(beta_i[i]);
        
        for(int j = 0; j < jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          
          // calculate the e_ij vector
          double e_ij_x = x[jj][0] - x[i][0];
          double e_ij_y = x[jj][1] - x[i][1];
          double e_ij_z = x[jj][2] - x[i][2];
          
          double e_r_2 = e_ij_x * e_ij_x + e_ij_y * e_ij_y + e_ij_z * e_ij_z;
          
          double alpha_j = sqrt(beta_i[jj]);
          
          // first sum
          if(rho_ji[i][j] > 0.0D && rho_i[i] > 0.0D && e_r_2 > 0.0D) {
            double v_rho_ji = rho_ji[i][j];
            double e_v_xi = e_ij_x * xi_i[i][0] + 
                            e_ij_y * xi_i[i][1] + 
                            e_ij_z * xi_i[i][2];
            
            double var = alpha_i * v_rho_ji / rho_i[i] * e_v_xi / e_r_2;
            
            f_RNG[i][0] += var * e_ij_x;
            f_RNG[i][1] += var * e_ij_y;
            f_RNG[i][2] += var * e_ij_z;
          }
          
          // second sum
          if(rho_ij[i][j] > 0.0D && rho_i[jj] > 0.0D && e_r_2 > 0.0D) {
            double v_rho_ij = rho_ij[i][j];
            double e_v_xi = e_ij_x * xi_i[jj][0] + 
                            e_ij_y * xi_i[jj][1] + 
                            e_ij_z * xi_i[jj][2];
            
            double var = alpha_j * v_rho_ij / rho_i[jj] * e_v_xi / e_r_2;
            
            f_RNG[i][0] -= var * e_ij_x;
            f_RNG[i][1] -= var * e_ij_y;
            f_RNG[i][2] -= var * e_ij_z;
          }
        }
        
        double v_Te = fdm->getT(x[i][0], x[i][1], x[i][2]);
        double var = eta_factor * sqrt(v_Te);
        f_RNG[i][0] *= v_struc * var;
        f_RNG[i][1] *= v_struc * var;
        f_RNG[i][2] *= v_struc * var;
      }
    }
  }
}

void FixEPH::force_gap() {
  double **v = atom->v;
  int *type = atom->type;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  
  if(eph_flag & Flag::FRICTION) {
    for(int i = 0; i < nlocal; ++i) {
      double v_star_x = 0.0;
      double v_star_y = 0.0;
      double v_star_z = 0.0;
      
      if(mask[i] & groupbit) {
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];
        
        for(int j = 0; j < jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];
          
          v_star_x -= v[jj][0] * rho_ji[i][j];
          v_star_y -= v[jj][1] * rho_ji[i][j];
          v_star_z -= v[jj][2] * rho_ji[i][j];
        }
        
        if(rho_i[i] > 0.0) {
          v_star_x /= rho_i[i];
          v_star_y /= rho_i[i];
          v_star_z /= rho_i[i];
          
          v_star_x += v[i][0];
          v_star_y += v[i][1];
          v_star_z += v[i][2];
          
          double v_rho = pow(rho_i[i], 4.0/3.0);
          double v_v = sqrt(v_star_x * v_star_x + 
                            v_star_y * v_star_y + 
                            v_star_z * v_star_z);
          
          double v_e_x = v_star_x / v_v;
          double v_e_y = v_star_y / v_v;
          double v_e_z = v_star_z / v_v;
          
          double var = grad_rho_i[i][0] * v_e_x + 
                       grad_rho_i[i][1] * v_e_y +
                       grad_rho_i[i][2] * v_e_z;
          
          //printf("%.6f\n", var);
          var *= v_alpha / v_rho;
          //printf("%.6f %.6f %.6f ", grad_rho_i[i][0], grad_rho_i[i][1], grad_rho_i[i][2]);
          //printf("%.6f %.6f %.6f ", v_star_x, v_star_y, v_star_z);
          //printf("%.6f %.6f %.6f %.6f %.6f\n", v_alpha, grad_rho_i[i][0], v_rho, rho_i[i], var);
          var = exp(var) * beta_i[i];
          
          double v_th = grad_rho_i[i][0] * grad_rho_i[i][0] +
                        grad_rho_i[i][1] * grad_rho_i[i][1] +
                        grad_rho_i[i][2] * grad_rho_i[i][2];
          v_th *= v_alpha / rho_i[i] / v_rho; // v_beta
                    
          if(v_v > v_th) {
            v_v -= v_th;
            
            var *= v_v;
            
            f_EPH[i][0] = -beta_factor * var * v_e_x;
            f_EPH[i][1] = -beta_factor * var * v_e_y;
            f_EPH[i][2] = -beta_factor * var * v_e_z;
          }
          else {
            f_EPH[i][0] = 0.0;
            f_EPH[i][1] = 0.0;
            f_EPH[i][2] = 0.0;
          }
        }
        else {
          f_EPH[i][0] = 0.0;
          f_EPH[i][1] = 0.0;
          f_EPH[i][2] = 0.0;
        } 
      }
    }
  }
}

void FixEPH::force_gapb() {
  double **v = atom->v;
  int *type = atom->type;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  
  if(eph_flag & Flag::FRICTION) {
    for(int i = 0; i < nlocal; ++i) {
      double v_star_x = 0.0;
      double v_star_y = 0.0;
      double v_star_z = 0.0;
      
      if(mask[i] & groupbit) {
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];
        
        for(int j = 0; j < jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];
          
          v_star_x -= v[jj][0] * rho_ji[i][j];
          v_star_y -= v[jj][1] * rho_ji[i][j];
          v_star_z -= v[jj][2] * rho_ji[i][j];
        }
        
        if(rho_i[i] > 0.0) {
          v_star_x /= rho_i[i];
          v_star_y /= rho_i[i];
          v_star_z /= rho_i[i];
          
          v_star_x += v[i][0];
          v_star_y += v[i][1];
          v_star_z += v[i][2];
          
          double v_rho = pow(rho_i[i], 4.0/3.0);
          double v_v = sqrt(v_star_x * v_star_x + 
                            v_star_y * v_star_y + 
                            v_star_z * v_star_z);
          
          double v_e_x = v_star_x / v_v;
          double v_e_y = v_star_y / v_v;
          double v_e_z = v_star_z / v_v;
          
          double var = grad_rho_i[i][0] * v_e_x + 
                       grad_rho_i[i][1] * v_e_y +
                       grad_rho_i[i][2] * v_e_z;
          var *= v_alpha / rho_i[i] / rho_i[i];
          //printf("%.6f %.6f %.6f ", grad_rho_i[i][0], grad_rho_i[i][1], grad_rho_i[i][2]);
          //printf("%.6f %.6f %.6f ", v_star_x, v_star_y, v_star_z);
          //printf("%.6f\n", var);
          var = exp(var) * beta_i[i];
          
          double v_th = grad_rho_i[i][0] * grad_rho_i[i][0] +
                        grad_rho_i[i][1] * grad_rho_i[i][1] +
                        grad_rho_i[i][2] * grad_rho_i[i][2];
          v_th *= v_alpha / rho_i[i] / v_rho; // v_beta
          
          //printf("%.6f %.6f\n", v_v, v_th);
          if(v_v > v_th) {
            v_v -= v_th;
            
            var *= v_v;
            
            f_EPH[i][0] = -beta_factor * var * v_e_x;
            f_EPH[i][1] = -beta_factor * var * v_e_y;
            f_EPH[i][2] = -beta_factor * var * v_e_z;
          }
          else {
            f_EPH[i][0] = 0.0;
            f_EPH[i][1] = 0.0;
            f_EPH[i][2] = 0.0;
          }
        }
        else {
          f_EPH[i][0] = 0.0;
          f_EPH[i][1] = 0.0;
          f_EPH[i][2] = 0.0;
        } 
      }
    }
  }
}

void FixEPH::post_force(int vflag) {
  //if(myID == 0) std::cout << "DEBUG post_force()" << std::endl;
  double **f = atom->f;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  int *numneigh = list->numneigh;
  
  //zero all arrays
  std::fill_n(&(rho_i[0]), nlocal, 0.0D);
  std::fill_n(&(beta_i[0]), nlocal, 0.0D);
  std::fill_n(&(xi_i[0][0]), 3 * nlocal, 0.0D);
  std::fill_n(&(w_i[0][0]), 3 * nlocal, 0.0D);
  std::fill_n(&(f_EPH[0][0]), 3 * nlocal, 0.0D);
  std::fill_n(&(f_RNG[0][0]), 3 * nlocal, 0.0D);
  std::fill_n(&(grad_rho_i[0][0]), 3 * nlocal, 0.0D);
  
  for(int i = 0; i < nlocal; ++i) {
    rho_i[i] = 0.0;
  }
  
  // generate random forces and distribute them
  if(eph_flag & Flag::RANDOM) {
    for(int i = 0; i < nlocal; ++i) {
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
    
    //MPI_Allreduce(MPI_IN_PLACE, &(xi_i[0][0]), 3*atom->natoms, 
    //  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  
  /* do a quick check for number of neighbours */
  if(update->ntimestep % 1000 == 0) {
    bool increase_neighs = false;
    for(int i = 0; i < nlocal; ++i) {
      if(numneigh[i] > rho_neigh) {
        increase_neighs = true;
        
        while(numneigh[i] > rho_neigh) {
          rho_neigh += 128;
        }
      }
    }
    
    if(increase_neighs) {
      std::cout << "WARNING GROWING ARRAY " << rho_neigh << std::endl;
      memory->grow(rho_ij, nlocal, rho_neigh, "eph:rho_ij");
      memory->grow(rho_ji, nlocal, rho_neigh, "eph:rho_ji");
    }
  }
  
  // zero density contributions
  std::fill_n(&(rho_ij[0][0]), rho_neigh * nlocal, 0.0D);
  std::fill_n(&(rho_ji[0][0]), rho_neigh * nlocal, 0.0D);
  
  // calculate the site densities, gradients (future) and beta(rho)
  calculate_environment();
  
  /* 
   * we have separated the model specific codes to make it more readable 
   * at the expense of code duplication 
   */
  
  if(eph_model == Model::TTM) {
    force_ttm();
  }
  else if(eph_model == Model::PRB) {
    force_prb();
  }
  else if(eph_model == Model::PRBMOD) {
    force_prbmod();
  }
  else if(eph_model == Model::ETA) {
    force_eta();
  }
  else if(eph_model == Model::GAP) {
    force_gap();
  }
  else if(eph_model == Model::GAPB) {
    force_gapb();
  }
  
  // second loop over atoms if needed
  // debug
  /*
  printf("%f %f %f\n", f[0][0], f[0][1], f[0][2]);
  printf("%f %f %f\n", f_EPH[0][0], f_EPH[0][1], f_EPH[0][2]);
  printf("%f %f %f\n", f_RNG[0][0], f_RNG[0][1], f_RNG[0][2]);
  */
 
  for(int i = 0; i < nlocal; i++) {
    if(eph_flag & Flag::FRICTION) {
      f[i][0] += f_EPH[i][0];
      f[i][1] += f_EPH[i][1];
      f[i][2] += f_EPH[i][2];
    }
    
    if(eph_flag & Flag::RANDOM) {
      f[i][0] += f_RNG[i][0];
      f[i][1] += f_RNG[i][1];
      f[i][2] += f_RNG[i][2];
    }
  }
  
  /** debug **/
  /*
  double fsumx = 0.0;
  double fsumy = 0.0;
  double fsumz = 0.0;
  
  for(int i = 0; i < nlocal; i++) {
    fsumx += f_RNG[i][0];
    fsumy += f_RNG[i][1];
    fsumz += f_RNG[i][1];
  }
  
  printf("%f %f %f\n", fsumx, fsumy, fsumz);
  */
}

void FixEPH::reset_dt() {
  // DEBUG
  //std::cout << "###### Timestep changed" << std::endl;
  
  // this should be correct if beta is in eV ps / Ang^2
  beta_factor = 1.0;
  eta_factor = sqrt(2.0 * force->boltz / update->dt);
  
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  
  fdm->setDt(update->dt);
  
  // beta_factor = 1.0 / force->ftm2v;
  // this is true for uniform distribution
  // eta_factor = sqrt(24.0 * force->boltz / update->dt / force->mvv2e) / force->ftm2v;
  
  // this is true for a gaussian distribution
  // eta_factor = sqrt(2.0 * force->boltz / update->dt / force->mvv2e) / force->ftm2v;
}

void FixEPH::grow_arrays(int ngrow) {
  //if(myID == 0) std::cout << "DEBUG grow_arrays()" << std::endl;
  memory->grow(f_EPH, ngrow, 3,"EPH:fEPH");
  memory->grow(f_RNG, ngrow, 3,"EPH:fRNG");
  
  memory->grow(beta_i, ngrow, "eph:beta_i");
  memory->grow(rho_i, ngrow, "eph:rho_i");
  
  memory->grow(w_i, ngrow, 3, "eph:w_i");
  memory->grow(xi_i, ngrow, 3, "eph:xi_i");
  
  memory->grow(rho_ij, ngrow, rho_neigh, "eph:rho_ij");
  memory->grow(rho_ji, ngrow, rho_neigh, "eph:rho_ji");
  
  memory->grow(grad_rho_i, ngrow, 3, "eph:grad_rho_i");
  
  // per atom values
  // we need only nlocal elements here
  memory->grow(array, ngrow, size_peratom_cols, "eph:array");
  array_atom = array;
}

double FixEPH::compute_vector(int i) {
  //if(myID == 0) std::cout << "DEBUG compute_vector()" << std::endl;
  if(i == 0)
    return Ee;
  else if(i == 1) {
    return fdm->calcTtotal();
  }
  
  return Ee;
}

/** TODO: There might be synchronisation issues here; maybe should add barrier for sync **/
int FixEPH::pack_forward_comm(int n, int *list, double *data, int pbc_flag, int *pbc) {
  int j, m;
  m = 0;
  if(state == FixState::RHO) {
    for(int i = 0; i < n; ++i) {
      j = list[i];
      data[m++] = rho_i[j];
    }
  }
  else if(state == FixState::BETA) {
    for(int i = 0; i < n; ++i) {
      j = list[i];
      data[m++] = beta_i[j];
    }
  }
  else if(state == FixState::XIX) {
    for(int i = 0; i < n; ++i) {
      j = list[i];
      data[m++] = xi_i[j][0];
    }
  }
  else if(state == FixState::XIY) {
    for(int i = 0; i < n; ++i) {
      j = list[i];
      data[m++] = xi_i[j][1];
    }
  }
  else if(state == FixState::XIZ) {
    for(int i = 0; i < n; ++i) {
      j = list[i];
      data[m++] = xi_i[j][2];
    }
  }
  else if(state == FixState::WX) {
    for(int i = 0; i < n; ++i) {
      j = list[i];
      data[m++] = w_i[j][0];
    }
  }
  else if(state == FixState::WY) {
    for(int i = 0; i < n; ++i) {
      j = list[i];
      data[m++] = w_i[j][1];
    }
  }
  else if(state == FixState::WZ) {
    for(int i = 0; i < n; ++i) {
      j = list[i];
      data[m++] = w_i[j][2];
    }
  }
  
  return m;
}

void FixEPH::unpack_forward_comm(int n, int first, double *data) {
  int m, last;
  m = 0;
  last = first + n;
  
  if(state == FixState::RHO) {
    for(int i = first; i < last; ++i) {
      rho_i[i] = data[m++];
    }
  }
  else if(state == FixState::BETA) {
    for(int i = first; i < last; ++i) {
      beta_i[i] = data[m++];
    }
  }
  else if(state == FixState::XIX) {
    for(int i = first; i < last; ++i) {
      xi_i[i][0] = data[m++];
    }
  }
  else if(state == FixState::XIY) {
    for(int i = first; i < last; ++i) {
      xi_i[i][1] = data[m++];
    }
  }
  else if(state == FixState::XIZ) {
    for(int i = first; i < last; ++i) {
      xi_i[i][2] = data[m++];
    }
  }
  else if(state == FixState::WX) {
    for(int i = first; i < last; ++i) {
      w_i[i][0] = data[m++];
    }
  }
  else if(state == FixState::WY) {
    for(int i = first; i < last; ++i) {
      w_i[i][1] = data[m++];
    }
  }
  else if(state == FixState::WZ) {
    for(int i = first; i < last; ++i) {
      w_i[i][2] = data[m++];
    }
  }
}

/** TODO **/
double FixEPH::memory_usage() {
    double bytes = 0;
    
    return bytes;
}

/* save temperature state after run */
void FixEPH::post_run() {
  if(myID == 0) fdm->saveState(Tstate);
}

/*
int FixEPH::pack_reverse_comm(int, int, double *) {
  
}

void FixEPH::unpack_reverse_comm(int, int *, double *) {
  
}
*/

