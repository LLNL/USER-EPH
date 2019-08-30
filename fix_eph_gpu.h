/**
 * This fix is a rewrite of our previous USER-EPH (a mod of fix_ttm) code and is based on our 
 * PRB 94, 024305 (2016) paper
 **/

/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

#ifdef FIX_EPH_GPU

#ifdef FIX_CLASS
FixStyle(eph/gpu,FixEPHGPU)
#else

#ifndef LMP_FIX_EPH_GPU_H
#define LMP_FIX_EPH_GPU_H

// external headers

// lammps headers

// internal headers
#include "fix_eph.h"
#include "eph_gpu.h"

namespace LAMMPS_NS {

// START GLOBAL THINGS 

// END GLOBAL THINGS

class FixEPHGPU : public FixEPH {
 public:
    FixEPHGPU(class LAMMPS *, int, char **); // constructor
    ~FixEPHGPU(); // destructor
    
    void grow_arrays(int) override; 
    void post_force(int) override;
    
  private:
    EPH_GPU eph_gpu;
    
    void calculate_environment(); 
    //void force_prl() override;
};

}
#endif
#endif
#endif
