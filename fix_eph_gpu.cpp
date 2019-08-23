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

  call_dummy(myID, nrPS);
}

// destructor
FixEPHGPU::~FixEPHGPU() 
{
}

#endif

