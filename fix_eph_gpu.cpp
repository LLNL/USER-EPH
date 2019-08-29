/*
 * Authors of the extension Artur Tamm, Alfredo Correa
 * e-mail: artur.tamm.work@gmail.com
 */

#ifdef FIX_EPH_GPU

// external headers
#include <mpi.h>

// lammps headers

// internal headers
#include "fix_eph_gpu.h"
#include "eph_gpu.h"

using namespace LAMMPS_NS;
using namespace FixConst;

// constructor
FixEPHGPU::FixEPHGPU(LAMMPS *lmp, int narg, char **arg) :
  FixEPH(lmp, narg, arg) 
{
  MPI_Comm_size(MPI_COMM_WORLD, &nrPS);
  MPI_Comm_rank(MPI_COMM_WORLD, &myID);

  
}

// destructor
FixEPHGPU::~FixEPHGPU() 
{
  
}

void FixEPHGPU::grow_arrays(int ngrow) 
{
  FixEPH::grow_arrays(ngrow);
}

void FixEPHGPU::post_force(int)
{
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
}

void FixEPHGPU::calculate_environment()
{
  
}

#endif

