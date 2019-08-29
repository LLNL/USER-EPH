
#ifdef FIX_EPH_GPU

#ifndef EPH_GPU_H
#define EPH_GPU_H

typedef double double3d[3];

// this thing will be passed to different kernels
struct EPH_GPU
{
  void* beta_gpu; // pointer to eph_beta_gpu object in gpu space
  // void* fdm_gpu; // pointer to the fdm thing in gpu space
  
  double r_cutoff;
  double r_cutoff_sq;
  double rho_cutoff;
  
  size_t nlocal;
  
  double3d *x_gpu; // atom positions [nlocal + nghost][3]
  double3d *v_gpu; // atom velocities [nlocal + nghost][3]
  double3d *f_gpu; // atom forces [nlocal + nghost][3]
  double3d *f_EPH_gpu; // friction force [nlocal][3]
  double3d *f_RNG_gpu; // random force [nlocal][3]
  
  double3d *xi_i_gpu; // random numbers [nlocal][3]
  double3d *w_i_gpu; // friction [nlocal][3]
  
  double *rho_i_gpu; // site density [nlocal]
};

EPH_GPU allocate_EPH_GPU(size_t n);
void deallocate_EPH_GPU(EPH_GPU& eph_gpu);

#endif
#endif
