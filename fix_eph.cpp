/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: arturt@ut.ee
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

#include <cstring>
#include <limits>
#include <algorithm>

// internal headers
#include "fix_eph.h"
#include "eph_beta.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/**
   * FixEPH arguments
   * arg[ 0] <- fix ID
   * arg[ 1] <- group
   * arg[ 2] <- name
   * arg[ 3] <- rng seed
   * arg[ 4] <- eph parameter; 0 disable all terms; 1 enable friction term; 2 enable random force; 4 enable fde;
   * arg[ 5] <- eph model selection; 1 standard ttm with beta(rho); 2 PRB version; 3 PRB mod; 4 Mason like; 5 Testing
   * arg[ 6] <- electronic density; this might be changed with new fdm model
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
  
  vector_flag = 1; // fix is able to output a vector compute
  size_vector = 2; // 2 elements in the vector
  global_freq = 1; // frequency for vector data
  extvector = 1; // 
  nevery = 1; // call end_of_step every step
  peratom_flag = 1;
  size_peratom_cols = 5;
  peratom_freq = 1;
  
  /** comm test **/
  comm_forward = 1;
  comm_reverse = 1;
  
  // special needed for velocity things to work
  comm->ghost_velocity = 1;
  
  // initialise rng
  seed = atoi(arg[3]);
  random = new RanMars(lmp, seed + myID);
  
  // read model behaviour parameters
  eph_flag = atoi(arg[4]);
  eph_model = atoi(arg[5]);
  
  // TODO: magic parameters for passing values TEMPORARY
  v_alpha = atof(arg[6]);
  v_rho = atof(arg[7]);
  v_Ce = atof(arg[8]);
  v_Te = atof(arg[9]);
  v_struc = 1.0000; // remove structure factor for now
  // initialise beta(rho)
  types = atom->ntypes;
  
  if(types > (narg - 17)) {
    error->all(FLERR, "Fix eph: number of types larger than provided in fix");
  }
  
  typeMap = new unsigned int[types];
  beta = new EPH_Beta(arg[16]);
  
  if(beta->getElementsNumber() < 1) {
    error->all(FLERR, "Fix eph: no elements found in input file");
  }
  
  rcutoff = beta->getCutoff();
  
  // do element mapping
  for(unsigned int i = 0; i < types; ++i) {
    typeMap[i] = std::numeric_limits<unsigned int>::max();
    for(unsigned int j = 0; j < beta->getElementsNumber(); ++j) {
      if(strcmp(beta->getName(j), arg[17+i]) == 0) typeMap[i] = j;
    }
    if(typeMap[i] > types) {
      error->all(FLERR, "Fix eph: elements not found in input file");
    }
  }
  
  // set force prefactors
  beta_factor = 1.0 / force->ftm2v;
  
  // uniform distribution
  // eta_factor = sqrt(24.0 * force->boltz / update->dt / force->mvv2e) / force->ftm2v;
  
  // gaussian distribution
  eta_factor = sqrt(2.0 * force->boltz / update->dt / force->mvv2e) / force->ftm2v;
  
  // allocate arrays for fix_eph
  f_EPH = NULL;
  f_RNG = NULL;

  w_i = NULL;

  beta_i = NULL;
  rho_i = NULL;
  array = NULL;
  
  v_xi = NULL;
  
  rho_neigh = 512;
  rho_ij = NULL;
  rho_ji = NULL;
  
  grad_rho_i = NULL;
  
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  
  test_array = NULL;

  list = NULL;
  
  Ee = 0.0; // electronic energy is zero in the beginning
  
  // we ignore all other input parameters at the moment
}

// destructor
FixEPH::~FixEPH() {
  delete random;
  delete[] typeMap;
  delete beta;
  
  atom->delete_callback(id, 0);
  
  memory->destroy(beta_i);
  memory->destroy(rho_i);
  
  memory->destroy(array);
  
  memory->destroy(f_EPH);
  memory->destroy(f_RNG);
  
  memory->destroy(rho_ij);
  memory->destroy(rho_ji);
  
  memory->destroy(grad_rho_i);
  
  memory->destroy(v_xi);
  memory->destroy(w_i);

  memory->destroy(test_array);
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
}

void FixEPH::init_list(int id, NeighList *ptr) {
  this->list = ptr;
}

int FixEPH::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  
  return mask;
}

void FixEPH::end_of_step() {
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
        E_local -= f_EPH[i][0]*v[i][0]*update->dt;
        E_local -= f_EPH[i][1]*v[i][1]*update->dt;
        E_local -= f_EPH[i][2]*v[i][2]*update->dt;
      }
    }
  }
  
  if(eph_flag & Flag::RANDOM) {
    for(unsigned int i = 0; i < nlocal; ++i) {
      if(mask[i] & groupbit) {
        E_local -= f_RNG[i][0]*v[i][0]*update->dt;
        E_local -= f_RNG[i][1]*v[i][1]*update->dt;
        E_local -= f_RNG[i][2]*v[i][2]*update->dt;
      }
    }
  }
  
  MPI_Allreduce(MPI_IN_PLACE, &E_local, 1, MPI_DOUBLE, MPI_SUM, world);
  if(v_alpha > 0.0)
    v_Te += E_local / v_Ce / v_rho; // v_Ce * v_rho [K / eV]
  
  Ee += E_local;
  
  // store per atom values into array
  for(unsigned int i = 0; i < nlocal; ++i) {
    if(mask[i] & groupbit) {
      tagint itag = tag[i] - 1;
      
      array[i][ 0] = rho_i[itag];
      array[i][ 1] = beta_i[itag];
      array[i][ 2] = w_i[itag][0];
      array[i][ 3] = w_i[itag][1];
      array[i][ 4] = w_i[itag][2];
    }
    else {
      array[i][ 0] = 0.0;
      array[i][ 1] = 0.0;
      array[i][ 2] = 0.0;
      array[i][ 3] = 0.0;
      array[i][ 4] = 0.0;
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
    tagint itag = tag[i] - 1;
    
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
          
          rho_i[itag] += v_rho_ji;
          
          double v_grad_rho = -beta->getDRho(typeMap[jtype-1], r);
          
          grad_rho_i[i][0] += v_grad_rho * delx / r;
          grad_rho_i[i][1] += v_grad_rho * dely / r;
          grad_rho_i[i][2] += v_grad_rho * delz / r;
        }
      }
      
      beta_i[itag] = beta->getBeta(typeMap[itype-1], rho_i[itag]);
    }
  }
  
  MPI_Allreduce(MPI_IN_PLACE, &(rho_i[0]), atom->natoms, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &(beta_i[0]), atom->natoms, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
}

void FixEPH::force_ttm() {
  double **v = atom->v;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  
  // create friction forces
  if(eph_flag & Flag::FRICTION) {
    for(int i = 0; i < nlocal; ++i) {
      if(mask[i] & groupbit) {
        tagint itag = tag[i] - 1;
        
        double var = -beta_factor * beta_i[itag];
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
        tagint itag = tag[i] - 1;
        
        double var = eta_factor * sqrt(v_Te * beta_i[itag]);
        f_RNG[i][0] = var * v_xi[itag][0];
        f_RNG[i][1] = var * v_xi[itag][1];
        f_RNG[i][2] = var * v_xi[itag][2];
      }
    }
  }  
}

void FixEPH::force_prb() {
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
        
        tagint itag = tag[i] - 1;
        
        double var = beta_i[itag];
        f_EPH[i][0] = var * v[i][0];
        f_EPH[i][1] = var * v[i][1];
        f_EPH[i][2] = var * v[i][2];
        
        for(int j = 0; j < jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];
          
          tagint jtag = tag[jj] - 1;
          
          if(rho_ji[i][j] > 0.0D && rho_i[itag] > 0.0D) {
            double v_rho_ji = rho_ji[i][j];
            double var = beta_i[itag] * v_rho_ji / rho_i[itag];
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
        tagint itag = tag[i] - 1;
        
        double var = eta_factor * sqrt(v_Te * beta_i[itag]);
        f_RNG[i][0] = var * v_xi[itag][0];
        f_RNG[i][1] = var * v_xi[itag][1];
        f_RNG[i][2] = var * v_xi[itag][2];
      }
    }
  }
}

void FixEPH::force_prbmod() {
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
        
        tagint itag = tag[i] - 1;
        
        double alpha_i = sqrt(beta_i[itag]);
        w_i[itag][0] = alpha_i * v[i][0];
        w_i[itag][1] = alpha_i * v[i][1];
        w_i[itag][2] = alpha_i * v[i][2];
        
        for(int j = 0; j < jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];
          
          tagint jtag = tag[jj] - 1;
          
          if (rho_ji[i][j] > 0.0D && rho_i[itag] > 0.0D) {
            double v_rho_ji = rho_ji[i][j];
            double var = alpha_i * v_rho_ji / rho_i[itag];
            w_i[itag][0] -= var * v[jj][0];
            w_i[itag][1] -= var * v[jj][1];
            w_i[itag][2] -= var * v[jj][2];
          }
        }
      }
    }
    
    MPI_Allreduce(MPI_IN_PLACE, &(w_i[0][0]), 3 * atom->natoms, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    // now calculate the forces
    for(int i = 0; i < nlocal; ++i) {
      if(mask[i] & groupbit) {
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];
        
        tagint itag = tag[i] - 1;
        
        double alpha_i = sqrt(beta_i[itag]);
        f_EPH[i][0] = alpha_i * w_i[itag][0];
        f_EPH[i][1] = alpha_i * w_i[itag][1];
        f_EPH[i][2] = alpha_i * w_i[itag][2];
        
        for(int j = 0; j < jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];
          
          tagint jtag = tag[jj] - 1;
          
          if(rho_ij[i][j] > 0.0D &&rho_i[jtag] > 0.0D) {
            double alpha_j = sqrt(beta_i[jtag]);
            double v_rho_ij = rho_ij[i][j];
            double var = alpha_j * v_rho_ij / rho_i[jtag];
            f_EPH[i][0] -= var * w_i[jtag][0];
            f_EPH[i][1] -= var * w_i[jtag][1];
            f_EPH[i][2] -= var * w_i[jtag][2];
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
        
        tagint itag = tag[i] - 1;
        
        double var = sqrt(beta_i[itag]);
        
        f_RNG[i][0] = var * v_xi[itag][0];
        f_RNG[i][1] = var * v_xi[itag][1];
        f_RNG[i][2] = var * v_xi[itag][2];
        
        for(int j = 0; j < jnum; j++) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          
          tagint jtag = tag[jj] - 1;
          
          if(rho_ij[i][j] > 0.0D && rho_i[jtag] > 0.0) {          
            double v_rho_ij = rho_ij[i][j];
            var = sqrt(beta_i[jtag]);
            var *= v_rho_ij / rho_i[jtag];
            
            f_RNG[i][0] -= var * v_xi[jtag][0];
            f_RNG[i][1] -= var * v_xi[jtag][1];
            f_RNG[i][2] -= var * v_xi[jtag][2];
          }
        }
        
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
        
        tagint itag = tag[i] - 1;
        
        double alpha_i = sqrt(beta_i[itag]);
        
        for(int j = 0; j < jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];
          
          // calculate the e_ij vector
          double e_ij_x = x[jj][0] - x[i][0];
          double e_ij_y = x[jj][1] - x[i][1];
          double e_ij_z = x[jj][2] - x[i][2];
          
          double e_r_2 = e_ij_x * e_ij_x + e_ij_y * e_ij_y + e_ij_z * e_ij_z;
          
          tagint jtag = tag[jj] - 1;
          double alpha_j = sqrt(beta_i[jtag]);
          
          // first sum
          if (rho_ji[i][j] > 0.0D && rho_i[itag] > 0.0D && e_r_2 > 0.0D) {
            double v_rho_ji = rho_ji[i][j];
            double e_v_v = e_ij_x * v[i][0] + 
                          e_ij_y * v[i][1] + 
                          e_ij_z * v[i][2];
            
            double var = alpha_i * v_rho_ji / rho_i[itag] * e_v_v / e_r_2;
            
            w_i[itag][0] += var * e_ij_x;
            w_i[itag][1] += var * e_ij_y;
            w_i[itag][2] += var * e_ij_z;
          }
          
          // second sum
          if (rho_ji[i][j] > 0.0D && rho_i[itag] > 0.0D && e_r_2 > 0.0D) {
            double v_rho_ji = rho_ji[i][j];
            double e_v_v = e_ij_x * v[jj][0] + 
                          e_ij_y * v[jj][1] + 
                          e_ij_z * v[jj][2];
            
            double var = alpha_i * v_rho_ji / rho_i[itag] * e_v_v / e_r_2;
            
            w_i[itag][0] -= var * e_ij_x;
            w_i[itag][1] -= var * e_ij_y;
            w_i[itag][2] -= var * e_ij_z;
          }
        }
      }
    }
    
    MPI_Allreduce(MPI_IN_PLACE, &(w_i[0][0]), 3 * atom->natoms, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    // now calculate the forces
    // f_i = W_ij w_j
    for(int i = 0; i < nlocal; ++i) {
      if(mask[i] & groupbit) {
        int itype = type[i];
        int *jlist = firstneigh[i];
        int jnum = numneigh[i];
        
        tagint itag = tag[i] - 1;
        
        double alpha_i = sqrt(beta_i[itag]);
        
        for(int j = 0; j < jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];
          
          // calculate the e_ij vector
          double e_ij_x = x[jj][0] - x[i][0];
          double e_ij_y = x[jj][1] - x[i][1];
          double e_ij_z = x[jj][2] - x[i][2];
          
          double e_r_2 = e_ij_x * e_ij_x + e_ij_y * e_ij_y + e_ij_z * e_ij_z;
          
          tagint jtag = tag[jj] - 1;
          double alpha_j = sqrt(beta_i[jtag]);
          
          // first sum
          if (rho_ji[i][j] > 0.0D && rho_i[itag] > 0.0D && e_r_2 > 0.0D) {
            double v_rho_ji = rho_ji[i][j];
            double e_v_v = e_ij_x * w_i[itag][0] + 
                          e_ij_y * w_i[itag][1] + 
                          e_ij_z * w_i[itag][2];
            
            double var = alpha_i * v_rho_ji / rho_i[itag] * e_v_v / e_r_2;
            
            f_EPH[i][0] += var * e_ij_x;
            f_EPH[i][1] += var * e_ij_y;
            f_EPH[i][2] += var * e_ij_z;
          }
          
          // second sum
          if (rho_ij[i][j] > 0.0D && rho_i[jtag] > 0.0D && e_r_2 > 0.0D) {
            double v_rho_ij = rho_ij[i][j];
            double e_v_v = e_ij_x * w_i[jtag][0] + 
                          e_ij_y * w_i[jtag][1] + 
                          e_ij_z * w_i[jtag][2];
            
            double var = alpha_j * v_rho_ij / rho_i[jtag] * e_v_v / e_r_2;
            
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
        
        tagint itag = tag[i] - 1;
        
        double alpha_i = sqrt(beta_i[itag]);
        
        for(int j = 0; j < jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          
          // calculate the e_ij vector
          double e_ij_x = x[jj][0] - x[i][0];
          double e_ij_y = x[jj][1] - x[i][1];
          double e_ij_z = x[jj][2] - x[i][2];
          
          double e_r_2 = e_ij_x * e_ij_x + e_ij_y * e_ij_y + e_ij_z * e_ij_z;
          
          tagint jtag = tag[jj] - 1;
          double alpha_j = sqrt(beta_i[jtag]);
          
          // first sum
          if(rho_ji[i][j] > 0.0D && rho_i[itag] > 0.0D && e_r_2 > 0.0D) {
            double v_rho_ji = rho_ji[i][j];
            double e_v_xi = e_ij_x * v_xi[itag][0] + 
                            e_ij_y * v_xi[itag][1] + 
                            e_ij_z * v_xi[itag][2];
            
            double var = alpha_i * v_rho_ji / rho_i[itag] * e_v_xi / e_r_2;
            
            f_RNG[i][0] += var * e_ij_x;
            f_RNG[i][1] += var * e_ij_y;
            f_RNG[i][2] += var * e_ij_z;
          }
          
          // second sum
          if(rho_ij[i][j] > 0.0D && rho_i[jtag] > 0.0D && e_r_2 > 0.0D) {
            double v_rho_ij = rho_ij[i][j];
            double e_v_xi = e_ij_x * v_xi[jtag][0] + 
                            e_ij_y * v_xi[jtag][1] + 
                            e_ij_z * v_xi[jtag][2];
            
            double var = alpha_j * v_rho_ij / rho_i[jtag] * e_v_xi / e_r_2;
            
            f_RNG[i][0] -= var * e_ij_x;
            f_RNG[i][1] -= var * e_ij_y;
            f_RNG[i][2] -= var * e_ij_z;
          }
        }
        
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
        
        tagint itag = tag[i] - 1;
        
        for(int j = 0; j < jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];
          
          v_star_x -= v[jj][0] * rho_ji[i][j];
          v_star_y -= v[jj][1] * rho_ji[i][j];
          v_star_z -= v[jj][2] * rho_ji[i][j];
        }
        
        if(rho_i[itag] > 0.0) {
          v_star_x /= rho_i[itag];
          v_star_y /= rho_i[itag];
          v_star_z /= rho_i[itag];
          
          v_star_x += v[i][0];
          v_star_y += v[i][1];
          v_star_z += v[i][2];
          
          double v_rho = pow(rho_i[itag], 4.0/3.0);
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
          //printf("%.6f %.6f %.6f %.6f %.6f\n", v_alpha, grad_rho_i[i][0], v_rho, rho_i[itag], var);
          var = exp(var) * beta_i[itag];
          
          double v_th = grad_rho_i[i][0] * grad_rho_i[i][0] +
                        grad_rho_i[i][1] * grad_rho_i[i][1] +
                        grad_rho_i[i][2] * grad_rho_i[i][2];
          v_th *= v_alpha / rho_i[itag] / v_rho; // v_beta
                    
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
        
        tagint itag = tag[i] - 1;
        
        for(int j = 0; j < jnum; ++j) {
          int jj = jlist[j];
          jj &= NEIGHMASK;
          int jtype = type[jj];
          
          v_star_x -= v[jj][0] * rho_ji[i][j];
          v_star_y -= v[jj][1] * rho_ji[i][j];
          v_star_z -= v[jj][2] * rho_ji[i][j];
        }
        
        if(rho_i[itag] > 0.0) {
          v_star_x /= rho_i[itag];
          v_star_y /= rho_i[itag];
          v_star_z /= rho_i[itag];
          
          v_star_x += v[i][0];
          v_star_y += v[i][1];
          v_star_z += v[i][2];
          
          double v_rho = pow(rho_i[itag], 4.0/3.0);
          double v_v = sqrt(v_star_x * v_star_x + 
                            v_star_y * v_star_y + 
                            v_star_z * v_star_z);
          
          double v_e_x = v_star_x / v_v;
          double v_e_y = v_star_y / v_v;
          double v_e_z = v_star_z / v_v;
          
          double var = grad_rho_i[i][0] * v_e_x + 
                       grad_rho_i[i][1] * v_e_y +
                       grad_rho_i[i][2] * v_e_z;
          var *= v_alpha / rho_i[itag] / rho_i[itag];
          //printf("%.6f %.6f %.6f ", grad_rho_i[i][0], grad_rho_i[i][1], grad_rho_i[i][2]);
          //printf("%.6f %.6f %.6f ", v_star_x, v_star_y, v_star_z);
          //printf("%.6f\n", var);
          var = exp(var) * beta_i[itag];
          
          double v_th = grad_rho_i[i][0] * grad_rho_i[i][0] +
                        grad_rho_i[i][1] * grad_rho_i[i][1] +
                        grad_rho_i[i][2] * grad_rho_i[i][2];
          v_th *= v_alpha / rho_i[itag] / v_rho; // v_beta
          
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
  double **f = atom->f;
  int *mask = atom->mask;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  int *numneigh = list->numneigh;
  
  //zero all arrays
  std::fill_n(&(rho_i[0]), (atom->natoms), 0.0D);
  std::fill_n(&(beta_i[0]), (atom->natoms), 0.0D);
  std::fill_n(&(v_xi[0][0]), 3 * (atom->natoms), 0.0D);
  std::fill_n(&(w_i[0][0]), 3 * (atom->natoms), 0.0D);
  std::fill_n(&(f_EPH[0][0]), 3 * nlocal, 0.0D);
  std::fill_n(&(f_RNG[0][0]), 3 * nlocal, 0.0D);
  std::fill_n(&(grad_rho_i[0][0]), 3 * nlocal, 0.0D);
  
  // generate random forces and distribute them
  if(eph_flag & Flag::RANDOM) {
    for(int i = 0; i < nlocal; ++i) {
      tagint itag = tag[i] - 1;
      
      if(mask[i] & groupbit) {
        v_xi[itag][0] = random->gaussian();
        v_xi[itag][1] = random->gaussian();
        v_xi[itag][2] = random->gaussian();
      }
    }
    
    MPI_Allreduce(MPI_IN_PLACE, &(v_xi[0][0]), 3*atom->natoms, 
      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
  // This can be merged into previous loop
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
}

void FixEPH::reset_dt() {
  beta_factor = 1.0 / force->ftm2v;
  // this is true for uniform distribution
  // eta_factor = sqrt(24.0 * force->boltz / update->dt / force->mvv2e) / force->ftm2v;
  
  // this is true for a gaussian distribution
  eta_factor = sqrt(2.0 * force->boltz / update->dt / force->mvv2e) / force->ftm2v;
}

void FixEPH::grow_arrays(int ngrow) {
  memory->grow(f_EPH, ngrow, 3,"EPH:fEPH");
  memory->grow(f_RNG, ngrow, 3,"EPH:fRNG");
  
  memory->grow(beta_i, atom->natoms, "eph:beta_i");
  memory->grow(rho_i, atom->natoms, "eph:rho_i");
  memory->grow(w_i, atom->natoms, 3, "eph:w_i");
  memory->grow(v_xi, atom->natoms, 3, "eph:v_xi");
  
  memory->grow(rho_ij, ngrow, rho_neigh, "eph:rho_ij");
  memory->grow(rho_ji, ngrow, rho_neigh, "eph:rho_ji");
  
  memory->grow(grad_rho_i, ngrow, 3, "eph:grad_rho_i");
  
  // per atom values
  memory->grow(array, ngrow, 5, "eph:array");
  array_atom = array;
}

double FixEPH::compute_vector(int i) {
  if(i == 0)
    return Ee;
  else if(i == 1)
    return v_Te;
  
  return Ee;
}

int FixEPH::pack_forward_comm(int, int *, double *, int, int *) {
  
}

void FixEPH::unpack_forward_comm(int, int, double *) {
  
}

int FixEPH::pack_reverse_comm(int, int, double *) {
  
}

void FixEPH::unpack_reverse_comm(int, int *, double *) {
  
}


