#include "cuda_runtime.h"
#include "curand_kernel.h"
#include "device_launch_parameters.h"
#include <stdio.h>

cudaError_t GetMem(void ** devPtr, size_t size) {
	cudaError_t cudaStatus;

	cudaStatus = cudaMalloc(devPtr, size);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto GetMemError;
    }
	cudaStatus = cudaMemset(*devPtr, 0, size);
	if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemset failed!");
        goto GetMemError;
    }

GetMemError:
    return cudaStatus;
}

int ShowGPUInfo(unsigned int mem) {
	/////////////////////////////////////////////////////////////////
	// determine the number of CUDA capable GPUs
	//
	int num_gpus=0; 
	cudaGetDeviceCount(&num_gpus);
	if(num_gpus < 1)
	{
		printf("no CUDA capable devices were detected\n");
		return 1;
	}

	/////////////////////////////////////////////////////////////////
	// display GPU configuration
	//
	printf("- number of CUDA devices:\t%d\n", num_gpus);
	for(int i = 0; i < num_gpus; i++)
	{
		cudaDeviceProp dprop;
		cudaGetDeviceProperties(&dprop, i);
		printf("%d: %s\n", i, dprop.name);
		printf("\tMaximum size of Grid: %d\n", *(dprop.maxGridSize));
		printf("\tMaximum threads/block: %u\n", dprop.maxThreadsPerBlock);
		printf("\tMaximum threads/processor: %u\n", dprop.maxThreadsPerMultiProcessor);
		printf("\tNumber of processors: %u\n", dprop.multiProcessorCount);
		printf("\tMaximum sizes of each dimension of a block:    %d x %d x %d\n",
	                    dprop.maxThreadsDim[0],
	                    dprop.maxThreadsDim[1],
	                    dprop.maxThreadsDim[2]);
		printf("\tMemoria global: %zd\n", dprop.totalGlobalMem);

		if(mem>dprop.totalGlobalMem) {
			printf("ERROR: Memoria insuficiente\n");
		}
	}
	printf("**********************************************\n");

	return 0;
}
