
#include <iostream>
#include <mpi.h>

#include "eph_fdm.h"

/* 
 * This example will create an electronic system and heat it from one end with
 * a laser source.
 */

// create an electronic system
EPH_FDM electrons {10000, 1, 1};

void laser_hit() {
  
}

int main(int args, char **argv) {
  // do the MPI_Init
  
  
  // do the MPI_Finalise
  
  return 0;
}

