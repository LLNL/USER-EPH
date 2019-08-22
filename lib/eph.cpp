
#include <cuda.h>
#include <cuda_runtime_api.h>

void run_dummy_test_cu();

void call_dummy() 
{
  double *f;
  int n = 1024;

  cudaMalloc((void**) &f, n*sizeof(double));
  run_dummy_test_cu();
  cudaFree(f);
}


