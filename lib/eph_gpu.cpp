
#include <cuda.h>
#include <cuda_runtime_api.h>

void run_dummy_test_cu(int, int);

void call_dummy(int myID, int nrPS) 
{
  double *f;
  int n = 1024;

  cudaMalloc((void**) &f, n*sizeof(double));
  run_dummy_test_cu(myID, nrPS);
  cudaFree(f);
}


