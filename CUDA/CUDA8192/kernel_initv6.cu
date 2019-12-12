#include "globals.h"
#include "cuda_runtime.h"
#include "curand_kernel.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <time.h>

extern curandState* devStates;
extern void free_GPU_memory ();


// INIT RANDOM VALUES, ONE PER EACH CELL:
__global__ void setup_kernel_v6 (curandState * state, unsigned long SEED) 
{        
    int i, idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;

    // Init Random vector
    for (i = idx; i<SIDE*SIDE; i += stride)
	    curand_init( (SEED << 20) + i, 0, 0, &state[i]);

}


int Launcher_init_v6 ()
{
    cudaError_t cudaStatus;

    // 0 - SETUP/INIT KERNEL 
    setup_kernel_v6 <<<num_blocks, num_threads>>> ( devStates, 1 );
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching setup_kernel_v6!\n", cudaStatus);
		goto Error;
    }
    
    return OK;

    Error:
        free_GPU_memory ();
        return (int) cudaStatus;
}


