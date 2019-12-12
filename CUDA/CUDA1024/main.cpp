//----------------------------------------------------
// PROGRAM laserCAv6: SIMULATION OF LASER DYNAMICS THROUGH A CELLULAR AUTOMATA
//
// Version 6.0
// laserCAv6 --> Two output files
// "results-SIDE.txt" with experiment traces and
// "comport.txt" with data about the simulation.
//
// Calculate Shannon's entropy and range of values
// of n for a set of parameter values
// Daniel Cagigas-Muñiz, Manuel Ramón López-Torres, José Luis Guisado-Lizar, 
// University of Seville. 2019
//-----------------------------------------------------
// PROGRAMMED IN C
//-----------------------------------------------------

// INCLUSION OF COMPILER'S LIBRARY FILES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "globals.h"
#include "curand_kernel.h"
#include "cuda_runtime.h"

// If the video showing electron evolution in the "e" matriz should be created, the 
// "#define __VIDEO" must be uncommented.
// Please, take into account that execution time could increase in an exponential way
// in this case depending on the SIDE constant defined in "globals.h".
//#define __VIDEO

#ifdef __unix
	#define fopen_s(pFile,filename,mode) ((*(pFile))=fopen((filename),(mode)))==NULL
	#define fscanf_s fscanf 
#endif


//char input_dat_file[] = "input400.dat";
//char input_dat_file[] = "input512.dat";
char input_dat_file[] = "input1024.dat";
//char input_dat_file[] = "input2048.dat";
//char input_dat_file[] = "input4096.dat";
//char input_dat_file[] = "input8192.dat";



void allocate_GPU_memory (SIMUL_DATA parametros);
void free_GPU_memory ();

int LoadSimulationData (SIMUL_DATA * parameter, int argc, char *argv[]);
int SaveResults1 (FILE *resultados, SIMUL_DATA parametros );
int SaveResults2 (FILE *resultados, SIMUL_DATA parametros );

extern cudaError_t GetMem(void ** devPtr, size_t size);
extern int ShowGPUInfo(unsigned int mem);
extern void print_parameters (SIMUL_DATA * parametros) ;
extern __global__ void setup_kernel_v6 (curandState * state, unsigned long SEED);

// DECLARATION OF GLOBAL VARIABLES
    int numecaja[10000], numfcaja[10000];
	int nemin;
	int nemax;
	int nfmin;
	int nfmax;
	int sumae, sumaf;
	double frec, entropiae, entropiaf;
	double muda1e, muda1f;
	double anchocajae, anchocajaf, rncajae, rncajaf;
	CELL_TYPE ncaja, ncajae, ncajaf;
	double nemedio, nfmedio;
	int sumanumecaja, sumanumfcaja;
    unsigned int total_e[TMAX], total_f[TMAX];


// DECLARATION OF GPU GLOBAL VARIABLES
	CELL_TYPE * GPU_f1 = 0;
	CELL_TYPE * GPU_f2 = 0;
	CELL_TYPE * GPU_f3 = 0;
	CELL_TYPE3 * GPU_tvf = 0;
	CELL_TYPE2 * GPU_e = 0;

    int *GPU_ocurrencias_MAXFOT=0;

	unsigned int * GPU_total_e=0;
	unsigned int * GPU_total_f=0;
    int  * GPU_nemin=0;
    int  * GPU_nemax=0;
    int  * GPU_nfmin=0;
    int  * GPU_nfmax=0;
    int  * GPU_sumae=0; 
    int  * GPU_sumaf=0; 
    int  * GPU_sumanumecaja=0; 
    int  * GPU_sumanumfcaja=0; 
    double *GPU_rncajae=0;
    double *GPU_rncajaf=0;  
    double *GPU_entropiae=0;
    double *GPU_entropiaf=0;

    //unsigned int LOCAL_total_e=0;
    //unsigned int LOCAL_total_f=0;
    
    int * GPU_numecaja=0;
    int * GPU_numfcaja=0;
    
	curandState* devStates;



//------------------------------------------------------------------
// PROGRAM MAIN FUNCTION:
int main (int argc, char *argv[]) {

	SIMUL_DATA parametros;
	FILE *comport;
	char results_filename[50];
    char comport_filename[50];
    FILE *resultados;

	int result;

    cudaError_t cudaStatus;

    printf ("Total GPU memory used: %ld Bytes\n", TOTAL_MEMORY);


	// OPEN FILE comport.txt
    sprintf(comport_filename, "comport-%d.txt", SIDE);
	if (fopen_s(&comport, comport_filename,"w+")) {
		perror ("Error when opening comport.txt:\n");
		return ERROR;
	}

	if (ERROR == LoadSimulationData(&parametros, argc, argv)) {
        perror ("Error when loading parameters\n");
		return ERROR;
	}

	fprintf( comport, "# LASERAC V. 6\n");
	fprintf( comport, "#-------------------------------------------\n");
	fprintf( comport, "# SIDE = %d\n",SIDE );
	fprintf( comport, "# MAXFOT = %d\n",MAXFOT );
	fprintf( comport, "# NUMDIV = %d\n",parametros.numdiv );
	fprintf( comport, "# TMAX = %d\n",parametros.tmax );
	fprintf( comport, "# TTRANS = %d\n",parametros.ttrans );
	fprintf( comport, "# NFOTRUIDO = %d    TYPE(0=INIC, 1=CONT)=%d\n",parametros.nfotruido,parametros.ntiporuido );
	fprintf( comport, "# JUMP THRESHOLD 1->0 = %d\n",parametros.threshold );
	fprintf( comport, "# SEED = %d\n",parametros.seed );
	fprintf( comport, "#-------------------------------------------\n");
	fprintf( comport, "# PB\t   TVe   TVf   Se    Sf    Nemedio    Nfmedio   Ocurr.MAXFOT \n");
	fprintf( comport, "#-----------------------------------------------------------------------------------------------------------\n");

    sprintf(results_filename, "results-%d.txt", SIDE);
	if ( fopen_s(&resultados, results_filename,"w+")) {
		perror ("Error openning results.txt:\n");
		return ERROR;
	}
    SaveResults1 (resultados, parametros);

	anchocajae = SIDE * SIDE / parametros.numdiv;
	anchocajaf = MAXFOT * SIDE * SIDE / parametros.numdiv;

    allocate_GPU_memory (parametros);

    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }


    result = Launcher_init_v6 ();

	for (parametros.pbombeo=parametros.pbombeomin; parametros.pbombeo<parametros.pbombeomax + 0.00001 ; parametros.pbombeo += parametros.pbombeopaso){
		for (parametros.tvelectron=parametros.tvelectronmin; parametros.tvelectron<=parametros.tvelectronmax; parametros.tvelectron += parametros.tvelectronpaso){
			for (parametros.tvfoton=parametros.tvfotonmin; parametros.tvfoton<=parametros.tvfotonmax; parametros.tvfoton += parametros.tvfotonpaso){
				
                fprintf (resultados, "PARAMETERS: PB , TVE, TVF :  %8.4f  %d  %d\n",parametros.pbombeo,parametros.tvelectron,parametros.tvfoton);

				// VALUES INITIALIZATION: 				
				parametros.occurrencesMAXFOT = 0;
				parametros.numele = 0;
				parametros.numfot = 0;
				entropiae = 0.0;
				entropiaf = 0.0;
				sumae = 0;
				sumaf = 0;

				// LAUNCH SIMULATION:
				result = Launcher_v6(&parametros, total_e, total_f);
				

                // SAVE RESULTS TO FILE:
                SaveResults2 (resultados, parametros);		
				

				// WRITE DATA IN comport.txt FILE:
				fprintf( comport, "%8.4f   %d   %d", parametros.pbombeo, parametros.tvelectron, parametros.tvfoton);
				fprintf( comport, "    %9.4f  %9.4f", entropiae, entropiaf );
				fprintf( comport, "    %8.2f  %8.2f", nemedio, nfmedio );
				fprintf( comport, "    %d\n", parametros.occurrencesMAXFOT);

				if ( ferror( comport ) ) {
					//fclose ( comport );
					perror ("Error in comport.txt:\n");
					return ERROR;
				}

            // END PARAMETERS LOOP
			}
		}
	}


	free_GPU_memory ();
	fclose (comport);
    fclose (resultados);

	return OK;

    Error:
/**
#ifdef __VIDEO
	cvReleaseImage( &img );
	cvReleaseVideoWriter( &writer );
#endif
/**/
    free_GPU_memory ();

    return (int) cudaStatus;
}


void allocate_GPU_memory (SIMUL_DATA parametros) {
	cudaError_t cudaStatus;

	// Init random numbers
	cudaStatus = cudaMalloc (&devStates, SIDE*SIDE*sizeof( curandState ) );
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Allocate GPU buffers:
    
		// Allocate buffers and variables for Shannon Entropy:
	cudaStatus = GetMem((void**) &GPU_nemin, sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }  
	cudaStatus = GetMem((void**) &GPU_nemax, sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }  
	cudaStatus = GetMem((void**) &GPU_nfmin, sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }  
	cudaStatus = GetMem((void**) &GPU_nfmax, sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }  
	cudaStatus = GetMem((void**) &GPU_sumae, sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }  
	cudaStatus = GetMem((void**) &GPU_sumaf, sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }  
	cudaStatus = GetMem((void**) &GPU_sumanumecaja, sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }  
	cudaStatus = GetMem((void**) &GPU_sumanumfcaja, sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }  
	cudaStatus = GetMem((void**) &GPU_rncajae, sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }  
	cudaStatus = GetMem((void**) &GPU_rncajaf, sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }  
	cudaStatus = GetMem((void**) &GPU_numecaja, (parametros.numdiv)*sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }  
    cudaStatus = GetMem((void**) &GPU_numfcaja, (parametros.numdiv)*sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }  
    cudaStatus = GetMem((void**) &GPU_entropiae, (parametros.numdiv)*sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }     
    cudaStatus = GetMem((void**) &GPU_entropiaf, (parametros.numdiv)*sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }     
    
    
		// Allocate rest of GPU buffers:
	cudaStatus = GetMem((void**) &GPU_e, SIDE*SIDE*sizeof(CELL_TYPE2));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = GetMem((void**) &GPU_f1, SIDE*SIDE*sizeof(CELL_TYPE));
    if (cudaStatus != cudaSuccess) {
        goto Error;
    }

    cudaStatus = GetMem((void**) &GPU_f2, SIDE*SIDE*sizeof(CELL_TYPE));
    if (cudaStatus != cudaSuccess) {
        goto Error;
    }

    cudaStatus = GetMem((void**) &GPU_tvf, SIDE*SIDE*MAXFOT * sizeof(CELL_TYPE3));
    if (cudaStatus != cudaSuccess) {
        goto Error;
    }

	cudaStatus = GetMem((void**) &GPU_total_e, TMAX*sizeof(unsigned int));
    if (cudaStatus != cudaSuccess) {
        goto Error;
    }

	cudaStatus = GetMem((void**) &GPU_total_f, TMAX*sizeof(unsigned int));
    if (cudaStatus != cudaSuccess) {
        goto Error;
    }

	cudaStatus = GetMem((void**) &GPU_ocurrencias_MAXFOT, sizeof(int));
    if (cudaStatus != cudaSuccess) {
        goto Error;
    }

    return;

Error:
    free_GPU_memory ();
}


void free_GPU_memory () {
	// Free gpu random numbers memory
	cudaFree(devStates);

    // Free buffers and variables for Shannon Entropy:
    cudaFree(GPU_nemin);
    cudaFree(GPU_nemax);
    cudaFree(GPU_nfmin);
    cudaFree(GPU_nfmax);
    cudaFree(GPU_sumae);
    cudaFree(GPU_sumaf);
    cudaFree(GPU_sumanumecaja);
    cudaFree(GPU_sumanumfcaja);
    cudaFree(GPU_rncajae);
    cudaFree(GPU_rncajaf);
    cudaFree (GPU_numecaja);
    cudaFree (GPU_numfcaja);
    cudaFree (GPU_entropiae);
    cudaFree (GPU_entropiaf);

    // Free rest of GPU buffers:
	cudaFree(GPU_e);
    cudaFree(GPU_f1);
    cudaFree(GPU_f2);
    cudaFree(GPU_tvf);
	cudaFree(GPU_total_e);
	cudaFree(GPU_total_f);
    cudaFree(GPU_ocurrencias_MAXFOT);
}


int LoadSimulationData(SIMUL_DATA * parameter, int argc, char *argv[]) {
	FILE *entrada;
    int res;
	
    // OPENS INPUT DATA FILE
    if (argc == 1) {
	    if ( fopen_s(&entrada, input_dat_file,"r")) {
		    printf ("Error openning input data file : %s\n", input_dat_file);
		    return ERROR;
	    }
    } else if (argc == 2) {
	    if ( fopen_s(&entrada, argv[1],"r")) {
		    printf ("Error openning input data file: %s\n", argv[1]);
		    return ERROR;
	    }
    } else {
        printf ("Error,\n\t Usage 1: \'laserCAv6 \' %s or \n\t Usage 2: \'laserCAv6 \' file-path \n", input_dat_file);
		return ERROR;

    }

    // SYSTEM DATA INPUT ENTER FROM entrada.dat FILE:
	res=fscanf(entrada, "%f", &((*parameter).pbombeomin));		// 0.012500
	res=fscanf(entrada, "%f", &((*parameter).pbombeomax));		// 0.012500
	res=fscanf(entrada, "%f", &((*parameter).pbombeopaso));		// 1.000000
	
	res=fscanf(entrada, "%d", &((*parameter).tvelectronmin));		// 180
	res=fscanf(entrada, "%d", &((*parameter).tvelectronmax));		// 180
	res=fscanf(entrada, "%d", &((*parameter).tvelectronpaso));	// 1
	
	res=fscanf(entrada, "%d", &((*parameter).tvfotonmin));		// 10
	res=fscanf(entrada, "%d", &((*parameter).tvfotonmax));		// 10
	res=fscanf(entrada, "%d", &((*parameter).tvfotonpaso));		// 1
	
	res=fscanf(entrada, "%d", &((*parameter).tmax));				// 1000
	res=fscanf(entrada, "%d", &((*parameter).ttrans));			// ¿999?
	res=fscanf(entrada, "%d", &((*parameter).numdiv));			// 1000
	
	res=fscanf(entrada, "%d", &((*parameter).nfotruido));			// 55
	//fscanf_s(entrada, "%d", &((*parameter).ntiporuido));
	(*parameter).ntiporuido = 1;
	res=fscanf(entrada, "%d", &((*parameter).threshold));			// 1
	res=fscanf(entrada, "%d", &((*parameter).seed));			// 10003
	
	fclose ( entrada );

	return OK;
}


int SaveResults1 (FILE *resultados, SIMUL_DATA parametros ) {
	fprintf( resultados, "# LASERAC V. 6\n");
	fprintf( resultados, "#-------------------------------------------\n");
	fprintf( resultados, "# SIDE = %d\n",SIDE );
	fprintf( resultados, "# MAXFOT = %d\n",MAXFOT );
	fprintf( resultados, "# TMAX = %d\n",parametros.tmax );
	fprintf( resultados, "# NFOTRUIDO = %d    TYPE(0=INIC, 1=CONT)=%d\n",parametros.nfotruido,parametros.ntiporuido );
	fprintf( resultados, "# JUMP THRESHOLD 1->0 = %d\n",parametros.threshold );
	fprintf( resultados, "# SEED = %d\n",parametros.seed );
	fprintf( resultados, "# PBOMBEO = %f\n",parametros.pbombeo );
	fprintf( resultados, "# TVELECTRON = %d\n",parametros.tvelectronmax );
	fprintf( resultados, "# TVFOTON = %d\n",parametros.tvfotonmax );
	fprintf( resultados, "#-------------------------------------------\n");
	fprintf( resultados, "#  t    ne    nf\n");
	fprintf( resultados, "#-------------------------------------------\n");

	return OK;
}

int SaveResults2 (FILE *resultados, SIMUL_DATA parametros ) {
    // WRITE RESULTS DATA IN resultados FILE
	for(int t=0; t<TMAX; t++) {
		fprintf(resultados, "%d" ,t);
		fprintf(resultados, "\t%d", total_e[t]);
		fprintf(resultados, "\t%d\n", total_f[t]);
	}

	fprintf(resultados, "# MAXFOT OCCURRENCES= %d\n\n",parametros.occurrencesMAXFOT );

	return OK;
}
