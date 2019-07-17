/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

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
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <cassert>

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
  
  // special: needed for velocity things to work
  comm->ghost_velocity = 1;
  
  // initialise rng
  seed = atoi(arg[3]);
  random = new RanMars(lmp, seed + myID);
  
  // read model behaviour parameters
  //eph_flag = atoi(arg[4]);
  eph_flag = strtol(arg[4], NULL, 0);
  /* debug */
  if(myID == 0) {
    std::cout << std::endl;
    std::cout << "Flag read: " << arg[4] << " -> " << eph_flag << '\n';
    if(eph_flag & Flag::FRICTION) std::cout << "Friction evaluation: ON\n";
    if(eph_flag & Flag::RANDOM) std::cout << "Random evaluation: ON\n";
    if(eph_flag & Flag::FDM) std::cout << "FDM grid solving: ON\n";
    if(eph_flag & Flag::NOINT) std::cout << "No integration: ON\n";
    if(eph_flag & Flag::NOFRICTION) std::cout << "No friction application: ON\n";
    if(eph_flag & Flag::NORANDOM) std::cout << "No random application: ON\n";
    std::cout << std::endl;
  }
  
  eph_model = atoi(arg[5]);
  
  if(myID == 0) {
    std::cout << std::endl;
    std::cout << "Model read: " << arg[5] << " -> " << eph_model << '\n';
    if(eph_model & Model::NONE) std::cout << "No friction model\n";
    else if(eph_model & Model::TTM) std::cout << "Standard TTM with beta(rho)\n";
    else if(eph_model & Model::PRB) std::cout << "PRB 94, 024305 (2016)\n";
    else if(eph_model & Model::PRLCM) std::cout << "PRL 120, 185501 (2018) CM model\n";
    else if(eph_model & Model::PRL) std::cout << "PRL 120, 185501 (2018)\n";
    else if(eph_model & Model::TESTING) std::cout << "Testing model\n";
    
    std::cout << std::endl;
  }
  
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
    
    sprintf(Tstate, "%s.restart", arg[13]);
  }
  
  Tfreq = atoi(arg[14]);
  if(Tfreq > 0) {
    sprintf(Tout, "%s", arg[15]);
  }
  
  // set the communicator
  fdm->setComm(world, myID, nrPS);
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
  rcutoffsq = rcutoff * rcutoff;
  rhocutoff = beta->getRhoCutoff();
  
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
  
  // allocate arrays for fix_eph
  f_EPH = nullptr;
  f_RNG = nullptr;
  
  w_i = nullptr;

  beta_i = nullptr;
  rho_i = nullptr;
  array = nullptr;
  
  xi_i = nullptr;
  
  rho_neigh = 512;
  rho_ij = nullptr;
  rho_ji = nullptr;
  
  list = nullptr;
  
  // NO ARRAYS BEFORE THIS
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  
  // zero arrays, so they would not contain garbage
  int nlocal = atom->nlocal;
  std::fill_n(&(rho_i[0]), nlocal, 0.0);
  std::fill_n(&(beta_i[0]), nlocal, 0.0);
  std::fill_n(&(xi_i[0][0]), 3 * nlocal, 0.0);
  std::fill_n(&(w_i[0][0]), 3 * nlocal, 0.0);
  std::fill_n(&(f_EPH[0][0]), 3 * nlocal, 0.0);
  std::fill_n(&(f_RNG[0][0]), 3 * nlocal, 0.0);
  
  for(int i = 0; i < nlocal; ++i) {
    for(int j = 0; j < size_peratom_cols; ++j) {
      array[i][j] = 0.0;
    }
  }
  
  Ee = 0.0; // electronic energy is zero in the beginning
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
  if(eph_flag & Flag::NOINT) return;

  double dtfm;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
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
  if(eph_flag & Flag::NOINT) return;
  
  double dtfm;

  double **v = atom->v;
  double **f = atom->f;
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
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
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
  
  Ee += E_local;

  for(unsigned int i = 0; i < nlocal; ++i) {
    if(mask[i] & groupbit) {
      array[i][ 0] = rho_i[i];
      array[i][ 1] = beta_i[i];
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

void FixEPH::calculate_environment() {
  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  
  // loop over atoms and their neighbours and calculate rho and beta(rho)
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
        
        double r_sq = delx*delx + dely*dely + delz*delz;
        
        if(r_sq < rcutoffsq) {
          double r = sqrt(r_sq);          
          double v_rho_ji = beta->getRho(typeMap[jtype-1], r);
          double v_rho_ij = beta->getRho(typeMap[itype-1], r);
          
          rho_ij[i][j] = v_rho_ij;
          rho_ji[i][j] = v_rho_ji;
          
          rho_i[i] += v_rho_ji;
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
          
          if(rho_ji[i][j] > 0.0 && rho_i[i] > 0.0) {
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

void FixEPH::force_prlcm() {
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
          
          if (rho_ji[i][j] > 0.0 && rho_i[i] > 0.0) {
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
          
          if(rho_ij[i][j] > 0.0 &&rho_i[jj] > 0.0) {
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
          
          if(rho_ij[i][j] > 0.0 && rho_i[jj] > 0.0) {          
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
  
void FixEPH::force_prl() {
  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int *mask = atom->mask;
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
          
          //double alpha_j = sqrt(beta_i[jj]);
          
          // first sum
          //if (rho_ji[i][j] > 0.0 && rho_i[i] > 0.0 && e_r_2 < rcutoffsq && e_r_2 > 0) {
          if (e_r_2 < rcutoffsq) {
            //double e_r = sqrt(e_r_2);
            
            //double v_rho_ji = beta->getRho(typeMap[jtype - 1], e_r);
            double v_rho_ji = rho_ji[i][j];
            
            double e_v_v = e_ij_x * v[i][0] + 
                          e_ij_y * v[i][1] + 
                          e_ij_z * v[i][2];
            
            double var = alpha_i * v_rho_ji / rho_i[i] * e_v_v / e_r_2;
            
            w_i[i][0] += var * e_ij_x;
            w_i[i][1] += var * e_ij_y;
            w_i[i][2] += var * e_ij_z;
          
          //}
          
          // second sum
          //if (rho_ji[i][j] > 0.0 && rho_i[i] > 0.0 && e_r_2 < rcutoffsq && e_r_2 > 0.0) {
            //double v_rho_ji = rho_ji[i][j];
            
            e_v_v = e_ij_x * v[jj][0] + 
                          e_ij_y * v[jj][1] + 
                          e_ij_z * v[jj][2];
            
            var = alpha_i * v_rho_ji / rho_i[i] * e_v_v / e_r_2;
            
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
          
          // first sum
          //if (rho_ji[i][j] > 0.0 && rho_i[i] > 0.0 && e_r_2 > 0.0) {
          if(e_r_2 < rcutoffsq) {
            //double e_r = sqrt(e_r_2);
            double alpha_j = sqrt(beta_i[jj]);
            //double v_rho_ji = beta->getRho(typeMap[jtype - 1], e_r);
            double v_rho_ji = rho_ji[i][j];
            double e_v_v = e_ij_x * w_i[i][0] + 
                          e_ij_y * w_i[i][1] + 
                          e_ij_z * w_i[i][2];
            
            double var = alpha_i * v_rho_ji / rho_i[i] * e_v_v / e_r_2;
            
            f_EPH[i][0] += var * e_ij_x;
            f_EPH[i][1] += var * e_ij_y;
            f_EPH[i][2] += var * e_ij_z;
          //}
          
          // second sum
          //if (rho_ij[i][j] > 0.0 && rho_i[jj] > 0.0 && e_r_2 > 0.0) {
            double v_rho_ij = rho_ij[i][j];
            //double v_rho_ij = beta->getRho(typeMap[itype - 1], e_r);
            e_v_v = e_ij_x * w_i[jj][0] + 
                          e_ij_y * w_i[jj][1] + 
                          e_ij_z * w_i[jj][2];
            
            var = alpha_j * v_rho_ij / rho_i[jj] * e_v_v / e_r_2;
            
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
          if(rho_ji[i][j] > 0.0 && rho_i[i] > 0.0 && e_r_2 > 0.0) {
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
          if(rho_ij[i][j] > 0.0 && rho_i[jj] > 0.0 && e_r_2 > 0.0) {
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

void FixEPH::force_testing() {};

void FixEPH::post_force(int vflag) {
  double **f = atom->f;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  int *numneigh = list->numneigh;
  
  //zero all arrays
  std::fill_n(&(rho_i[0]), nlocal, 0.0);
  std::fill_n(&(beta_i[0]), nlocal, 0.0);
  std::fill_n(&(xi_i[0][0]), 3 * nlocal, 0.0);
  std::fill_n(&(w_i[0][0]), 3 * nlocal, 0.0);
  std::fill_n(&(f_EPH[0][0]), 3 * nlocal, 0.0);
  std::fill_n(&(f_RNG[0][0]), 3 * nlocal, 0.0);
  
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
  }
  
  /* do a quick check for number of neighbours */
  // TODO: remove this an use function calls instead
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
  std::fill_n(&(rho_ij[0][0]), rho_neigh * nlocal, 0.0);
  std::fill_n(&(rho_ji[0][0]), rho_neigh * nlocal, 0.0);
  
  // calculate the site densities, gradients (future) and beta(rho)
  calculate_environment();
  
  /* 
   * we have separated the model specific codes to make it more readable 
   * at the expense of code duplication 
   */
  switch(eph_model) {
    case Model::TTM: force_ttm();
      break;
    case Model::PRB: force_prb();
      break;
    case Model::PRLCM: force_prlcm();
      break;
    case Model::PRL: force_prl();
      break;
    case Model::TESTING: force_testing();
      break;
    default: throw;
  }
  
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

void FixEPH::reset_dt() {
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
  memory->grow(f_EPH, ngrow, 3,"EPH:fEPH");
  memory->grow(f_RNG, ngrow, 3,"EPH:fRNG");
  
  memory->grow(beta_i, ngrow, "eph:beta_i");
  memory->grow(rho_i, ngrow, "eph:rho_i");
  
  memory->grow(w_i, ngrow, 3, "eph:w_i");
  memory->grow(xi_i, ngrow, 3, "eph:xi_i");
  
  memory->grow(rho_ij, ngrow, rho_neigh, "eph:rho_ij");
  memory->grow(rho_ji, ngrow, rho_neigh, "eph:rho_ji");
  
  // per atom values
  // we need only nlocal elements here
  memory->grow(array, ngrow, size_peratom_cols, "eph:array");
  array_atom = array;
}

double FixEPH::compute_vector(int i) {
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
  switch(state) {
    case FixState::RHO:
      for(int i = 0; i < n; ++i) {
        j = list[i];
        data[m++] = rho_i[j];
      }
      break;
    case FixState::BETA:
      for(int i = 0; i < n; ++i) {
        j = list[i];
        data[m++] = beta_i[j];
      }
      break;
    case FixState::XIX:
      for(int i = 0; i < n; ++i) {
        j = list[i];
        data[m++] = xi_i[j][0];
      }
      break;
    case FixState::XIY:
      for(int i = 0; i < n; ++i) {
        j = list[i];
        data[m++] = xi_i[j][1];
      }
      break;
    case FixState::XIZ:
      for(int i = 0; i < n; ++i) {
        j = list[i];
        data[m++] = xi_i[j][2];
      }
      break;
    case FixState::WX:
      for(int i = 0; i < n; ++i) {
        j = list[i];
        data[m++] = w_i[j][0];
      }
      break;
    case FixState::WY:
      for(int i = 0; i < n; ++i) {
        j = list[i];
        data[m++] = w_i[j][1];
      }
      break;
    case FixState::WZ:
      for(int i = 0; i < n; ++i) {
        j = list[i];
        data[m++] = w_i[j][2];
      }
      break;
    default:
      break;
  }
  
  /*
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
  */
  return m;
}

void FixEPH::unpack_forward_comm(int n, int first, double *data) {
  int m, last;
  m = 0;
  last = first + n;
  
  switch(state) {
    case FixState::RHO:
      for(int i = first; i < last; ++i) {
        rho_i[i] = data[m++];
      }
      break;
    case FixState::BETA:
      for(int i = first; i < last; ++i) {
        beta_i[i] = data[m++];
      }
      break;
    case FixState::XIX:
      for(int i = first; i < last; ++i) {
        xi_i[i][0] = data[m++];
      }
      break;
    case FixState::XIY:
      for(int i = first; i < last; ++i) {
        xi_i[i][1] = data[m++];
      }
      break;
    case FixState::XIZ:
      for(int i = first; i < last; ++i) {
        xi_i[i][2] = data[m++];
      }
      break;
    case FixState::WX:
      for(int i = first; i < last; ++i) {
        w_i[i][0] = data[m++];
      }
      break;
    case FixState::WY:
      for(int i = first; i < last; ++i) {
        w_i[i][1] = data[m++];
      }
      break;
    case FixState::WZ:
      for(int i = first; i < last; ++i) {
        w_i[i][2] = data[m++];
      }
      break;
    default:
      break;
  }
  
  /*
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
  */
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
