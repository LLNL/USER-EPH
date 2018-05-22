
#include <iostream>
#include <mpi.h>

#include "eph_fdm.h"

/* 
 * This example will create an electronic system and heat it from one end with
 * a laser source.
 */

// create an electronic system
EPH_FDM electrons {1000000, 1, 1};

void laser_hit() {
  
}

int main(int args, char **argv) {
  // do the MPI_Init
  electrons.setBox( // in Ang
    -3000.0, 3000.0,
    0.0, 28.16, // 8 x 3.52 current shock simulation
    0.0, 28.16
    );
  electrons.setConstants(300.0, 3.5e-6, 1.0, 0.1248); // T_e, C_e, rho_e, kappa_e
  electrons.setSteps(1); // minimum nr. of steps
  electrons.setDt(0.0001); // in ps
  
  char a {'a'};
  while(a != 'e')
    std::cin >> a;
  
  // do the MPI_Finalise
  
  return 0;
}

