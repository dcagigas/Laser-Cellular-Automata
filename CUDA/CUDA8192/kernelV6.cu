#include "globals.h"
#include "cuda_runtime.h"
#include "curand_kernel.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 


extern int numecaja[10000], numfcaja[10000];
extern int nemin;
extern int nemax;
extern int nfmin;
extern int nfmax;
extern int sumae, sumaf;
extern double frec, entropiae, entropiaf;
extern double muda1e, muda1f;
extern double anchocajae, anchocajaf, rncajae, rncajaf;
extern CELL_TYPE ncaja, ncajae, ncajaf;
extern double nemedio, nfmedio;
extern int sumanumecaja, sumanumfcaja;

// GPU extern variables:
extern	CELL_TYPE * GPU_f1;
extern	CELL_TYPE * GPU_f2;
extern	CELL_TYPE * GPU_f3;
extern	CELL_TYPE3 * GPU_tvf;
extern	CELL_TYPE2 * GPU_e;
extern int *GPU_ocurrencias_MAXFOT;
extern unsigned int * GPU_total_e;
extern unsigned int * GPU_total_f;
extern int  * GPU_nemin;
extern int  * GPU_nemax;
extern int  * GPU_nfmin;
extern int  * GPU_nfmax;
extern int  * GPU_sumae; 
extern int  * GPU_sumaf;  
extern int  * GPU_sumanumecaja; 
extern int  * GPU_sumanumfcaja; 
extern double *GPU_rncajae;
extern double *GPU_rncajaf;  
extern double *GPU_entropiae;
extern double *GPU_entropiaf;

extern int * GPU_numecaja;
extern int * GPU_numfcaja;
    
extern curandState* devStates;



// Functions and kernels definition:

    // Init variables
__global__ void init_time_step_v6 (CELL_TYPE2 * e, CELL_TYPE * f1, CELL_TYPE * f2, CELL_TYPE3 * tvf, int *ocurrencias_MAXFOT, int *nemin, int *nemax, int *nfmin, int *nfmax, int *sumae, int *sumaf, int *numecaja, int *numfcaja, int numdiv, int *sumanumecaja, int *sumanumfcaja, double *entropiae, double *entropiaf);

    // Shannon entropy kernels
__global__ void do_shannon_entropy_time_step_v6 (unsigned int * total_e, unsigned int * total_f, int *nemin, int *nemax, int *nfmin, int *nfmax, int *sumae, int *sumaf, double anchocajae, double anchocajaf, int *numecaja, int *numfcaja, double *rncajae, double *rncajaf);
__global__ void finish_shannon_entropy_time_step_v6 (SIMUL_DATA parametros, int *numecaja, int *numfcaja, int *sumanumecaja, int *sumanumfcaja, double *entropiae, double *entropiaf);

    // Main process kernels
__global__ void PhotonNoiseKernel_v6 (CELL_TYPE * f1, CELL_TYPE * f2, CELL_TYPE3 * tvf, SIMUL_DATA parametros, curandState * state, int *ocurrencias_MAXFOT);
__global__ void DecayKernel_v6 (CELL_TYPE2 * e, CELL_TYPE * f1,  CELL_TYPE * f2, CELL_TYPE3 * tvf);
__global__ void PumpingKernel_v6 (CELL_TYPE2 * e, CELL_TYPE * f1, CELL_TYPE * f2, CELL_TYPE3 * tvf, SIMUL_DATA parametros, unsigned int * total_e, unsigned int * total_f, curandState * state, int *ocurrencias_MAXFOT);

void print_parameters (SIMUL_DATA * parametros) ;



#ifdef __VIDEO
    #include "opencv2/core/core.hpp"
    #include "opencv2/highgui/highgui.hpp"
#endif

CELL_TYPE f[SIDE][SIDE];
CELL_TYPE e[SIDE][SIDE];

cudaError_t GetMem(void ** devPtr, size_t size); // Specified in cuda_functions.cu 
void free_GPU_memory ();                         // Specified in main.cpp


__global__ void init_time_step_v6 (unsigned int *total_e, unsigned int *total_f, CELL_TYPE2 *e, CELL_TYPE *f1, CELL_TYPE *f2, CELL_TYPE3 *tvf, int *ocurrencias_MAXFOT, int *nemin, int *nemax, int *nfmin, int *nfmax, int *sumae, int *sumaf, int *numecaja, int *numfcaja, int numdiv, int *sumanumecaja, int *sumanumfcaja, double *entropiae, double *entropiaf) 
{
    int i, k, idx = blockIdx.x * blockDim.x + threadIdx.x;  
    int stride = gridDim.x * blockDim.x;

    for (i = idx; i<numdiv; i += stride) {
        numecaja[i] = 0;
        numfcaja[i] = 0;
        total_e[i] = 0;
        total_f[i] = 0;
    }

    __syncthreads();

    // Init global variables:
    for (i = idx; i<SIDE*SIDE; i += stride) {
        e[i] = 0;
        f1[i] = 0;
        f2[i] = 0;
        for (k=0; k<MAXFOT; k++) {
            tvf[i * MAXFOT + k] = 0;
        }
    }
      
    if (idx == 0) {
        *ocurrencias_MAXFOT = 0;
        *entropiae = 0;
        *entropiaf = 0;
		*nemax = 0;
		*nemin = SIDE * SIDE;
		*nfmax = 0;
		*nfmin = MAXFOT * (*nemin);
		*sumae = 0;
		*sumaf = 0;
        *sumanumecaja = 0;
        *sumanumfcaja = 0;
    }
}


__global__ void do_shannon_entropy_time_step_v6 (unsigned int * total_e, unsigned int * total_f, int *nemin, int *nemax, int *nfmin, int *nfmax, int *sumae, int *sumaf, double anchocajae, double anchocajaf, int *numecaja, int *numfcaja, double *rncajae, double *rncajaf) 
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;  
    //double muda1e, muda1f;
    int ncajae, ncajaf;
      
    if (idx == 0) {
        if (*total_e < *nemin) *nemin = *total_e;
        if (*total_e > *nemax) *nemax = *total_e;
        if (*total_f < *nfmin) *nfmin = *total_f;
        if (*total_f > *nfmax) *nfmax = *total_f;

        // ADD TO TOTAL PHOTONS AND ELECTRONS 
        *sumae += *total_e;
        *sumaf += *total_f;
						
        // OCCURRENCE OF EACH BOX ('caja') FOR THE CALCULATION OF THE SHANNON ENTROPY
	    //muda1e = modf (*total_e / (*anchocajae), rncajae);
        modf ( (*total_e) / anchocajae, rncajae);
		ncajae = (int) (*rncajae);
		numecaja [ncajae] ++;

		//muda1f = modf (*total_f / (*anchocajaf), rncajaf);
        modf ((*total_f) / anchocajaf, rncajaf);
		ncajaf = (int) (*rncajaf);
		numfcaja [ncajaf] ++;
    }
}


__global__ void finish_shannon_entropy_time_step_v6 (SIMUL_DATA parametros, int *numecaja, int *numfcaja, int *sumanumecaja, int *sumanumfcaja, double *entropiae, double *entropiaf) {
    int i, idx = blockIdx.x * blockDim.x + threadIdx.x;  
    int stride = gridDim.x * blockDim.x;
//    int tid = threadIdx.x;
    double frec;

/**
	__shared__ int sumanumecaja_local[num_threads];
	__shared__ int sumanumfcaja_local[num_threads];


    if (idx < parametros.numdiv) {
	    sumanumecaja_local[tid] = numecaja[idx];
	    sumanumfcaja_local[tid] = numfcaja[idx];
    }
    __syncthreads();
    for (i = idx+stride; i<(parametros.numdiv); i += stride) {
	    sumanumecaja_local[tid] += numecaja[idx];
	    sumanumfcaja_local[tid] += numfcaja[idx];
    }
    __syncthreads();
        atomicAdd (sumanumecaja, sumanumecaja_local[tid]);
        atomicAdd (sumanumfcaja, sumanumfcaja_local[tid]);


/**
	for(int half=num_threads/2; half>0; half=half/2) {
		if (tid<half) {
			sumanumecaja_local[tid] += sumanumecaja_local[tid+half];
			sumanumfcaja_local[tid] += sumanumfcaja_local[tid+half];
		}
		__syncthreads();
	}

	if (tid==0)  {
        atomicAdd(sumanumecaja, sumanumecaja_local[0]);
		atomicAdd(sumanumfcaja, sumanumfcaja_local[0]);
		
	}
/**/


    for (i = idx; i<(parametros.numdiv); i += stride) {
        atomicAdd (sumanumecaja, numecaja [i]);
        atomicAdd (sumanumfcaja, numfcaja [i]);
	}



    __syncthreads();

    for (i = idx; i<(parametros.numdiv); i += stride) {
		if (numecaja [i] != 0){
			frec = (double) numecaja [i] / (*sumanumecaja);
            atomicAdd (entropiae, frec * log(frec));
			//entropiae = entropiae + frec * log(frec);
		}
		if (numfcaja [i] != 0){
			frec = (double) numfcaja [i] / (*sumanumfcaja);
            atomicAdd (entropiaf, frec * log(frec));
			//entropiaf = entropiaf + frec * log(frec);
		}
    }
}


// KERNEL LAUNCHER CALLED IN MAIN FUNCTION:
int Launcher_v6 (SIMUL_DATA * parametros, unsigned int * total_e, unsigned int * total_f) {
	cudaError_t cudaStatus;

#ifdef __VIDEO
	IplImage* img = 0;
	CvVideoWriter * writer;
    char video_file_name[80];

    sprintf (video_file_name, "Exit_video-%d.avi", SIDE);
	writer = cvCreateVideoWriter(video_file_name, CV_FOURCC('D', 'I', 'V', 'X') , 24 , cvSize(SIDE, SIDE), true);
	img = cvCreateImage(cvSize(SIDE, SIDE), IPL_DEPTH_8U, 3);
#endif


    init_time_step_v6 <<<num_blocks, num_threads>>> (GPU_total_e, GPU_total_f, GPU_e, GPU_f1, GPU_f2, GPU_tvf, GPU_ocurrencias_MAXFOT, GPU_nemin, GPU_nemax, GPU_nfmin, GPU_nfmax, GPU_sumae, GPU_sumaf, GPU_numecaja, GPU_numfcaja, parametros->numdiv, GPU_sumanumecaja, GPU_sumanumfcaja, GPU_entropiae, GPU_entropiaf);
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching shannon_entropy_time_step_v4!\n", cudaStatus);
		goto Error;
	}


	// CALLS TO MAIN KERNELS (1 TO 3):
    // TIME LOOP ('TMAX' TIME STEPS).
	for(int t=0; t<TMAX; t++) { 

        // 1 - PHOTON NOISE KERNEL:
		PhotonNoiseKernel_v6 <<<num_blocks, num_threads>>>(GPU_f1, GPU_f2, GPU_tvf, *parametros, devStates, GPU_ocurrencias_MAXFOT);
		// cudaDeviceSynchronize waits for the kernel to finish, and returns
		// any errors encountered during the launch.
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching fotones_ruido!\n", cudaStatus);
			goto Error;
		}

        // 2 - DECAY KERNEL:
		DecayKernel_v6 <<<num_blocks, num_threads>>> (GPU_e, GPU_f1, GPU_f2, GPU_tvf);
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching PhotonDecayKernel!\n", cudaStatus);
			goto Error;
		}

        // 3 - ELECTRON PUMPING AND STIMULATED EMISSION KERNEL:
		PumpingKernel_v6 <<<num_blocks, num_threads>>> (GPU_e, GPU_f1, GPU_f2, GPU_tvf, *parametros, GPU_total_e+t, GPU_total_f+t, devStates, GPU_ocurrencias_MAXFOT);
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching PumpingKernel!\n", cudaStatus);
			goto Error;
		}

        // 4 - CALCULATIONS FOR THE SHANNON ENTROPY 
        if (t > (parametros->ttrans) ) {
            do_shannon_entropy_time_step_v6 <<<1, 1>>> (GPU_total_e+t, GPU_total_f+t, GPU_nemin, GPU_nemax, GPU_nfmin, GPU_nfmax, GPU_sumae, GPU_sumaf, anchocajae, anchocajaf, GPU_numecaja, GPU_numfcaja, GPU_rncajae, GPU_rncajaf); 
		    cudaStatus = cudaDeviceSynchronize();
		    if (cudaStatus != cudaSuccess) {
			    fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching PumpingKernel!\n", cudaStatus);
			    goto Error;
		    }
        }

        // MATRIXES 'f1' AND 'f2' ARE EXCHANGED FOR THE NEXT LOOP ITERATION.
		GPU_f3 = GPU_f1;
		GPU_f1 = GPU_f2;
		GPU_f2 = GPU_f3;	

#ifdef __VIDEO
		cudaStatus = cudaMemcpy(e, GPU_e, SIDE * SIDE * sizeof(CELL_TYPE2), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}

		// VIDEO WRITE:
		cvRectangle(img, cvPoint(0,0),cvPoint(SIDE, SIDE),cvScalar(0),-1);
		for (int i=0; i < SIDE; i++) {
			for (int j = 0; j < SIDE; j++) {
				cvLine(img, cvPoint(j, i), cvPoint(j, i), cvScalar(255*(e[i][j]>0), 0, 0)); // (f[i][j]>0)*(127 + f[i][j]*128/MAXFOT)
			}
		}
		cvWriteFrame( writer, img );
#endif

	} // END TIME LOOP


    finish_shannon_entropy_time_step_v6 <<<num_blocks, num_threads>>> (*parametros, GPU_numecaja, GPU_numfcaja, GPU_sumanumecaja, GPU_sumanumfcaja, GPU_entropiae, GPU_entropiaf);
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching PumpingKernel!\n", cudaStatus);
		goto Error;
	}

	cudaMemcpy(total_e, GPU_total_e, TMAX*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(total_f, GPU_total_f, TMAX*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(&sumae, GPU_sumae, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&sumaf, GPU_sumaf, sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(&entropiae, GPU_entropiae, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&entropiaf, GPU_entropiaf, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(parametros->occurrencesMAXFOT), GPU_ocurrencias_MAXFOT, sizeof(int), cudaMemcpyDeviceToHost);


    // 'n' RELATIVE RANGE CALCULATION ("DELTAREL")
	nemedio = (double) sumae / ( (parametros->tmax) - (parametros->ttrans) );
	nfmedio = (double) sumaf / ( (parametros->tmax) - (parametros->ttrans) );

	entropiae *= -1.44269504;
	entropiaf *= -1.44269504;



Error:
#ifdef __VIDEO
	cvReleaseImage( &img );
	cvReleaseVideoWriter( &writer );
#endif

    free_GPU_memory ();

    return (int) cudaStatus;
}



//------------------------------------------------------------------
// NOISE PHOTONS FUNCTION:
__global__ void PhotonNoiseKernel_v6 (CELL_TYPE * f1, CELL_TYPE * f2, CELL_TYPE3 * tvf, SIMUL_DATA parametros, curandState * state, int *ocurrencias_MAXFOT) {
    int i, idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;
	int j;
	curandState localState;
	int x, y;
    CELL_TYPE f1_old;

    // SETS 'NFOTRUIDO' PHOTONS IN RANDOM CELL POSIITIONS
    for (i = idx; i<parametros.nfotruido; i += stride) {

        // IT GETS TWO RAMDON NUMBRES (curand_uniform(&localState)) BETWEEN 0 AND 1.
		localState = state[i];
		x = curand_uniform(&localState)*SIDE - 1;
		y = curand_uniform(&localState)*SIDE - 1;
		state[i] = localState;
	
		j = y*SIDE+x;

        // IT CAN HAPPENDS THAT TWO HW THREADS CAN TRY TO UPDATE THE SAME CELL (i.e. HAVE THE SAME 'j')
        // THIS MEANS THAT TWO PHOTONS CAN BE CREATED IN THE SAME CELL. THIS IS POSSIBLE.
        // TO PREVENT RACE CONDICITIONS, FIRST THE PHOTON COUNTER f1 MUST BE INCREMENTED IN AN ATOMIC WAY.
        // THIS FIRST ATOMIC INSTRUCTION MAKES ROOM FOT THE NEXT HW THREAD THAT WANTS TO UPDATE THE SAME CELL.
        f1_old = atomicAdd(&f1[j],1);
        
		//if (f1[j] < MAXFOT){
		if (f1_old < MAXFOT){
			// IT CREATES A NEW PHOTON IN THAT POSITION
			tvf[j*MAXFOT+f1[f1_old]] = parametros.tvfoton;
            f2[j] = f1[j];
		} else {
            atomicSub(&f1[j],1);
            atomicAdd(ocurrencias_MAXFOT,1);
        }
	}
}


// PHOTON DECAY FUNCTION:
__global__ void DecayKernel_v6 (CELL_TYPE2 *e, CELL_TYPE * f1, CELL_TYPE * f2, CELL_TYPE3 * tvf) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;
    int i;
	int it;
	int j, k, n;

    for (i = idx; i<SIDE*SIDE; i += stride) {
		k=0;
		j=0;
		n = f1[i];
		it = i * MAXFOT;

        // PHOTON DECAY:
        // 'n' HAS MAXFOT AS MAXIMUN VALUE (tipical 10).
        // THIS CODE MUST BE SEQUENTIAL BECAUSE IF A PHOTON DECAYS TO 0, 
        // ITS POSITION IN THE 'tvf' VECTOR MUST BE OCCUPIED BY THE NEXT PHOTON: THE PHOTON VECTOR IS SHIFTED.
        // THEREFORE, THERE ARE DEPENDENCIES BETWEEN CONSECUTIVE ELEMENTS IN THE 'tvf' VECTOR.
		while(j<n) {
			tvf[it+k] = tvf[it+j]-1;
			if(tvf[it+k]>0) k++;
			j++;
		}
		f1[i]=k;
        f2[i]=f1[i];

	    // ELECTRON DECAY:
	    if (e[i] > 0) {
            e[i] = e[i] - 1;
	    }
	}
}


// PUMPING AND STIMULATED EMISSION FUNCTIONS:
__global__ void PumpingKernel_v6 (CELL_TYPE2 * e, CELL_TYPE * f1, CELL_TYPE * f2, CELL_TYPE3 * tvf, SIMUL_DATA parametros, unsigned int * total_e, unsigned int * total_f, curandState * state, int *ocurrencias_MAXFOT) {
	int W, E, N, S, NE, SE, SW, NW;
    int neighbors = 0;
    int i, idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;
	int tid = threadIdx.x;
	curandState localState;
	CELL_TYPE partial_e=0;
	CELL_TYPE partial_f=0;

	__shared__ CELL_TYPE total_f_local[num_threads];
	__shared__ CELL_TYPE total_e_local[num_threads];

    for (i = idx; i<SIDE*SIDE; i += stride) {
		// PUMPING
		if (e[i] == 0) {
			localState = state[i];
		
            // IT GETS TWO RAMDON NUMBRES (curand_uniform(&localState)) BETWEEN 0 AND 1.
			if (curand_uniform(&localState) < parametros.pbombeo){
				e[i] = parametros.tvelectron;
			}

			state[i] = localState;
		} else {

		// STIMULATED EMISSION:
		// (e[i] > 0)
			int l, r;
			int u, d;
			int y = i / SIDE;
			int x = i % SIDE;

			l = x>0?x-1:SIDE-1;
			r = x<SIDE-1?x+1:0;
			u = y>0?y-1:SIDE-1;
			d = y<SIDE-1?y+1:0;
		
			E = y*SIDE + r;
			W = y*SIDE + l;
			N = u*SIDE + x;
			S = d*SIDE + x;
			NE = u*SIDE + r;
			NW = u*SIDE + l;
			SE = d*SIDE + r;
			SW = d*SIDE + l;
		
			neighbors =	f1[N] + f1[NE] + f1[E] + f1[SE] +
						f1[S] +	f1[SW] + f1[W] +f1[NW] + 
						f1[i];

			if (neighbors>parametros.threshold && f2[i]<MAXFOT) {
				tvf[i * MAXFOT + f2[i]] = parametros.tvfoton;
				//f2[i]++;
                f2[i] = f2[i] + 1;
				e[i] = 0;	
			} else if (f2[i]>=MAXFOT) {
                atomicAdd(ocurrencias_MAXFOT,1);
            }
		}
		
		partial_f += f2[i];
		partial_e += (e[i]>1);
	}

    // PHOTON AND ELECTRON COUNT:
    // REDUCTION: THIS IMROVES COMPUTING TIME TO LOG INSTEAD OF LINEAR.
/**/
	total_f_local[tid] = partial_f;
	total_e_local[tid] = partial_e;
    __syncthreads();

	for(int half=num_threads/2; half>0; half=half/2) {
		if (tid<half) {
			total_f_local[tid] += total_f_local[tid+half];
			total_e_local[tid] += total_e_local[tid+half];
		}
		__syncthreads();
	}

	if (tid==0)  {
        atomicAdd(total_f, total_f_local[0]);
		atomicAdd(total_e, total_e_local[0]);
	}
/**/



/** 
    __syncthreads();
    for (i = idx; i<SIDE*SIDE; i += stride) {
        atomicAdd (total_f, partial_f);
        atomicAdd (total_e, partial_e);
	}

/**/

}




void print_parameters (SIMUL_DATA * parametros) {
	printf ("\t pbombeomin: %f\n", (parametros->pbombeomin) );
	printf ("\t pbombeomax: %f\n", (parametros->pbombeomax) );
	printf ("\t pbombeopaso: %f\n", (parametros->pbombeopaso) );
	printf ("\t ----------- \n");
	printf ("\t tvelectronmin: %d\n", (parametros->tvelectronmin) );
	printf ("\t tvelectronmax: %d\n", (parametros->tvelectronmax) );
	printf ("\t tvelectronpaso: %d\n", (parametros->tvelectronpaso) );
	printf ("\t ----------- \n");	
	printf ("\t tvfotonmin: %d\n", (parametros->tvfotonmin) );
	printf ("\t tvfotonmax: %d\n", (parametros->tvfotonmax) );
	printf ("\t tvfotonpaso: %d\n", (parametros->tvfotonpaso) );
	printf ("\t ----------- \n");	
	printf ("\t tmax: %d\n", (parametros->tmax) );
	printf ("\t ttrans: %d\n", (parametros->ttrans) );
	printf ("\t numdiv: %d\n", (parametros->numdiv) );
	printf ("\t ----------- \n");
	printf ("\t nfotruido: %d\n", (parametros->nfotruido) );
	printf ("\t threshold: %d\n", (parametros->threshold) );
	printf ("\t seed: %d\n", (parametros->seed) );
	printf ("\t ----------- \n");	
}


