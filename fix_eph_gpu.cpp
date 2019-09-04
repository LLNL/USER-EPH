/*
 * Authors of the extension Artur Tamm, Alfredo Correa
 * e-mail: artur.tamm.work@gmail.com
 */

#ifdef FIX_EPH_GPU

// external headers
#include <mpi.h>
#include <iostream>
#include <cstdio>

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
#include "fix_eph.h"
#include "fix_eph_gpu.h"
#include "eph_gpu.h"

using namespace LAMMPS_NS;
using namespace FixConst;

// constructor
FixEPHGPU::FixEPHGPU(LAMMPS *lmp, int narg, char **arg) :
  FixEPH(lmp, narg, arg) 
{
  eph_gpu = allocate_EPH_GPU(beta, types, type_map);
  eph_gpu.groupbit = groupbit;
  
  x_mirror = nullptr;
}

// destructor
FixEPHGPU::~FixEPHGPU() 
{
  deallocate_EPH_GPU(eph_gpu);
  delete[] x_mirror;
}

void FixEPHGPU::grow_arrays(int ngrow) 
{
  FixEPH::grow_arrays(ngrow);
  eph_gpu.grow(ngrow);
  
  if(x_mirror != nullptr) delete[] x_mirror;
  x_mirror = new double3d[ngrow];
}

void FixEPHGPU::post_force(int)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  
  //zero all arrays
  // TODO: remove these in the future
  std::fill_n(&(w_i[0][0]), 3 * nlocal, 0.0);
  std::fill_n(&(xi_i[0][0]), 3 * nlocal, 0.0);
  std::fill_n(&(f_EPH[0][0]), 3 * nlocal, 0.0);
  std::fill_n(&(f_RNG[0][0]), 3 * nlocal, 0.0);
  
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
  // push coordinates to GPU this is ugly
  eph_gpu.nlocal = nlocal;
  eph_gpu.nghost = nghost;
  
  int ntotal = nlocal + nghost;
  
  for(size_t i = 0; i < ntotal; ++i)
  {
    x_mirror[i][0] = x[i][0];
    x_mirror[i][1] = x[i][1];
    x_mirror[i][2] = x[i][2];
  }
  
  //cpu_to_device_EPH_GPU((void*) eph_gpu.x_gpu, (void*) x_mirror, 3*ntotal*sizeof(double));
  cpu_to_device_EPH_GPU((void*) eph_gpu.x_gpu, (void*) x[0], 3*ntotal*sizeof(double));
  cpu_to_device_EPH_GPU((void*) eph_gpu.type_gpu, (void*) type, ntotal*sizeof(int));
  cpu_to_device_EPH_GPU((void*) eph_gpu.mask_gpu, (void*) mask, ntotal*sizeof(int));
  
  calculate_environment();
  
  device_to_cpu_EPH_GPU((void*) rho_i, (void*) eph_gpu.rho_i_gpu, nlocal*sizeof(double));
  
  state = FixState::RHO;
  comm->forward_comm_fix(this);
  
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

void FixEPHGPU::calculate_environment()
{
  calculate_environment_gpu(eph_gpu);
}

#endif

