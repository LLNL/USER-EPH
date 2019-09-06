
#ifdef FIX_EPH_GPU

#ifndef EPH_GPU_H
#define EPH_GPU_H

#include <vector>
#include <memory>

#include "eph_beta.h"

// Old-school object-oriented c with structures
typedef double double3d[3];

// this thing will be passed to different kernels
struct EPH_GPU
{
  size_t neighmask;
  
  double eta_factor;
  
  void* beta_gpu; // pointer to eph_beta_gpu object in gpu space
  // void* fdm_gpu; // pointer to the fdm thing in gpu space
  
  double r_cutoff;
  double r_cutoff_sq;
  double rho_cutoff;
  
  int threads;
  int blocks;

  int nlocal;
  int nghost;
  
  int groupbit;
  
  int n; // array sizes
  int n_neigh; // total number of elements in neighs_gpu
  
  int *number_neigh_gpu;
  int *index_neigh_gpu;
  int *neighs_gpu;
  
  int types;
  int* type_map_gpu;
  
  int *type_gpu;
  int *mask_gpu;
  double3d *x_gpu; // atom positions [n + nghost][3]
  double3d *v_gpu; // atom velocities [n + nghost][3]
  
  double3d *f_EPH_gpu; // friction force [n][3]
  double3d *f_RNG_gpu; // random force [n][3]
  
  double3d *xi_i_gpu; // random numbers [n][3]
  double3d *w_i_gpu; // friction [n][3]
  
  double3d *T_e_i_gpu; // electronic temperatures at each atomic size [n]
  double *rho_i_gpu; // site density [n]
  
  // instances in cpu memory
  void* beta;
  // void *fdm;
  
  // member functions  
  void grow(size_t ngrow);
  void grow_neigh(size_t ngrow);
};

// this is done on purpose
EPH_GPU allocate_EPH_GPU(Beta& eph_beta, int types, int* type_map);
void deallocate_EPH_GPU(EPH_GPU& eph_gpu);

void cpu_to_device_EPH_GPU(void* dst, void* src, size_t n);
void device_to_cpu_EPH_GPU(void* dst, void* src, size_t n);
 
void calculate_environment_gpu(EPH_GPU& eph_gpu);

void force_prl_stage1_gpu(EPH_GPU& eph_gpu);
void force_prl_stage2_gpu(EPH_GPU& eph_gpu);
void force_prl_stage3_gpu(EPH_GPU& eph_gpu);

#endif
#endif
