
#include "eph_gpu.h"

#include <mpi.h>

int main(int args, char **argv) 
{
  MPI_Init(&args, &argv);

  int myID;
  int nrPS;
  
  MPI_Comm_size(MPI_COMM_WORLD, &nrPS);
  MPI_Comm_rank(MPI_COMM_WORLD, &myID);

  call_dummy(myID, nrPS);

  MPI_Finalize();
  return 0;
}

