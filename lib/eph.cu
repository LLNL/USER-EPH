
#include <cuda.h>
#include <cuda_runtime_api.h>

__global__
void dummy_test()
{
  int i = 0;
  while(1) 
  { 
    ++i; 
  }
}


void run_dummy_test_cu() 
{
  dummy_test<<<1, 256>>>();
}

