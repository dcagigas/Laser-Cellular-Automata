//----------------------------------------------------
// LASER DYNAMICS SIMULATION USING CELLULAR AUTOMATAS:
// OPENMP AND SIMD OPTIMIZATIONS FOR THE PAPER: 
//"Efficiently Building Discrete Simulations on Multi-Core and Many-Core Architectures"
//
// THIS IS THE OPTIMIZED CODIFICATION.
// HOT ZONES WERE DETECTED IN THE BASIC VERSION laseracv10_omp_basic
// AND SEVERAL OPTIMIZATIONS WERE INSERTED. THey have been marked with symbols //@
//
// file "results.txt" stores photon and electron evolution for 
// a bunch of input parameters. 
// file "timing_results.txt" stores the timing of different pieces of the code. 
// random number generator is PCG 
// Several constants (#define) can be commented or not at the beginning of this file to
// use I/O functions for Visual Studio or GCC, to enable openmp, print inner function timing, etc.
//
// In addition, Shannon entropy is computed at the end of the simulation
//
// By Jose Luis Guisado Lizar, Daniel Cagigas-Mugniz, Ramon Lopez-Torres 
// and Fernando Diaz-del-Rio (Univ. Sevilla. 2019)

//---------------------------------------------------------------------------------
// CONFIGURATION FOR VARIOUS OPTIMIZATIONS 
//---------------------------------------------------------------------------------
#define ENABLE_OMP

// USE ONE OF THESE TWO STYLES:
#define GCC_STYLE
//#define VS_STYLE

#ifdef VS_STYLE
#define fopen_GCC_VS(a,b,c)  (fopen_s(&a, b, c) != 0)
#define fscanf_GCC_VS(a,b,c)  fscanf_s(a, b, c)
//#define PRINT_ENABLE
//#define OUTER_RULE_TIMING
#endif

#ifdef GCC_STYLE
#define fopen_GCC_VS(a,b,c)   ((a = fopen(b, c)) == NULL)
#define fscanf_GCC_VS(a,b,c)  fscanf(a, b, c) != 1
// #define PRINT_ENABLE
//#define OUTER_RULE_TIMING
#endif

#define MIN_NUM_THREADS 16
#define MAX_NUM_THREADS 16

// In addition a pregenerated vector of random binary numbers (instead of rand/PCG) 
//  can be used for functions Decay_v5() , Pumping_v5();   //@
//#define RANDOM_VECT_ENABLE

//@  SIZE OF THE VECTORS
typedef  short int  CELL_TYPE;
typedef CELL_TYPE  laser_e_t;
typedef CELL_TYPE  laser_lifet_f_t;
typedef unsigned char laser_f1_t;

//---------------------------------------------------------------------------------
// CONFIGURATION PARAMETERS  TO SIMULATE VARIOUS CELLULAR AUTOMATA 
//---------------------------------------------------------------------------------

#define SIDE 1024 // NUMBER OF CELLS PER SIDE 

//char input_dat_file[] = "input4096.dat";
//char input_dat_file[] = "input2048.dat";
char input_dat_file[] = "input1024.dat";
//char input_dat_file[] = "input512.dat";
//char input_dat_file[] = "input400.dat"; 

#define MAX_PHOT 10   // MAX NUM PHOTONS FOR EACH BOX 
#define TMAX 1000  

//---------------------------------------------------------------------------------
// ADDITIONAL PARAMETERS  AND INCLUDE FILES 
//---------------------------------------------------------------------------------
#define OK 0        // RETURNED VALUE IF NO ERRORS
#define ERROR 1     // RETURNED VALUE IF ERRORS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <omp.h>
#include "pcg_basic.h"

//---------------------------------------------------------------------------------
// GLOBAL VARIABLES 
//---------------------------------------------------------------------------------
// Using PCG Random Number Generation for C.
// (Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>)
// Using a local Random Number Generator
pcg32_random_t rng_Pumping[SIDE][SIDE];  //@ one seed stream per cell
pcg32_random_t rng_PhotonNoise;

#ifdef RANDOM_VECT_ENABLE 
// random vector to prevent calling rand() function inside the main loop:
bool  random_vect[TMAX][SIDE][SIDE]; //@
#endif

laser_e_t  e[SIDE][SIDE];
laser_f1_t f1[SIDE][SIDE], f2[SIDE][SIDE];
laser_lifet_f_t lifet_f[MAX_PHOT + 1][SIDE][SIDE]; //@optim_lifet_f

laser_f1_t (*pf1)[SIDE];
laser_f1_t (*pf2)[SIDE];

float ppumping;

int lifet_electron, lifet_photon;
int tmax, nphot_noise, threshold;
int events_max_phot;

int seed;

int tstep; 

FILE *input;
FILE *eph_results;
FILE *behaviour;
FILE *timing_results;

// function prototypes 
void PhotonNoise_v5(void), Decay_v5(void), Pumping_v5(void);
int neighboors(int, int);

//------------------------------------------------------------------
// MAIN
int main(void)
{
	float ppumpingmin, ppumpingmax, ppumping_step;
	int lifet_electronmin, lifet_electronmax, lifet_electron_step;
	int lifet_photonmin, lifet_photonmax, lifet_photon_step;
	int ntype_noise;
	int numdiv, ttrans;
	int numele, num_phot;
	double frec, entropye, entropyf;
	int nemin;
	int nemax;
	int nfmin;
	int nfmax;
	int sume, sumf;
	double ne_mean, nf_mean;
	double muda1e, muda1f;
	double box_width_e, box_width_f, rn_box_e, rn_box_f;
	int n_box_, n_box_e, n_box_f;
	int sumnume_box_, sumnumf_box_;
	int nume_box_[10000], numf_box_[10000];
	char fin = 'a';

	double t_sum_all_steps;
	double t_sum_all_param = 0.0;
	int nof_executions_param = 0;
	double t0_step, t_inc_step;
	
	double tcomplete_start = omp_get_wtime();

	pf1 = &f1[0];
	pf2 = &f2[0];

	//--------------------------------------------------
	// FILE DATA INPUT/OUTPUT 
	// 
	if ((fopen_GCC_VS(input, input_dat_file, "r"))) {  
		perror("Error opening input_dat_file:\n");
		return ERROR;
	}
	if ((fopen_GCC_VS(eph_results, "results.txt", "w+"))) {  
		perror("Error opening results.txt:\n");
		return ERROR;
	}
	if ((fopen_GCC_VS(behaviour, "behaviour.txt", "w+"))) {  
		perror("Error opening behaviour.txt:\n");
		return ERROR;
	}
	if ((fopen_GCC_VS(timing_results, "timing_results.txt", "a+"))) { 
		perror("Error opening timing_results.txt:\n");
		return ERROR;
	}

	fscanf_GCC_VS(input, "%f", &ppumpingmin);

	//printf("PPUMPING MAX= ");
	fscanf_GCC_VS(input, "%f", &ppumpingmax);

	//printf("PPUMPING STEP= ");
	fscanf_GCC_VS(input, "%f", &ppumping_step);

	//printf(" EXCITED ELECTRON LIFETIME : lifet_ELECTRON MIN= ");
	fscanf_GCC_VS(input, "%d", &lifet_electronmin);

	//printf("lifet_ELECTRON MAX= ");
	fscanf_GCC_VS(input, "%d", &lifet_electronmax);

	//printf("lifet_ELECTRON STEP= ");
	fscanf_GCC_VS(input, "%d", &lifet_electron_step);

	//printf("PHOTON LIFETIME : lifet_FOTON MIN = ");
	fscanf_GCC_VS(input, "%d", &lifet_photonmin);

	//printf("lifet_PHOTON MAX= ");
	fscanf_GCC_VS(input, "%d", &lifet_photonmax);

	//printf("lifet_PHOTON STEP= ");
	fscanf_GCC_VS(input, "%d", &lifet_photon_step);
	//printf("-------------------------------------------\n");

	//printf("MAX STEPS = ");
	fscanf_GCC_VS(input, "%d", &tmax);

	//printf("TRANSIENT = ");
	fscanf_GCC_VS(input, "%d", &ttrans);

	//printf("NOF DIVISIONS = ");
	fscanf_GCC_VS(input, "%d", &numdiv);

	//printf("NOF NOISE PHOTONES = ");
	fscanf_GCC_VS(input, "%d", &nphot_noise);

	//printf("threshold = ");
	fscanf_GCC_VS(input, "%d", &threshold);

	//printf("SEED FOR RANDOM NUMBERS ");
	fscanf_GCC_VS(input, "%d", &seed);
	//printf("-------------------------------------------\n");

	//NOISE TYPE: not implemented yet
	ntype_noise = 1;

	fprintf(eph_results, "# LASERAC V. 10 OPTIM\n");
	fprintf(eph_results, "#-------------------------------------------\n");
	fprintf(eph_results, "# SIDE = %d\n", SIDE);
	fprintf(eph_results, "# MAX_PHOT = %d\n", MAX_PHOT);
	fprintf(eph_results, "# TMAX = %d\n", tmax);
	fprintf(eph_results, "# nphot_noise = %d    type(0=INIC, 1=CONT)=%d\n", nphot_noise, ntype_noise);
	fprintf(eph_results, "# threshold 1->0 = %d\n", threshold);
	fprintf(eph_results, "# SEED = %d\n", seed);
	fprintf(eph_results, "# PUMPING PROBABILITY MIN = %f\n", ppumpingmin);
	fprintf(eph_results, "# PUMPING PROBABILITY MAX = %f\n", ppumpingmax);
	fprintf(eph_results, "# lifet_ELECTRON = %d\n", lifet_electron);
	fprintf(eph_results, "# lifet_PHOTON = %d\n", lifet_photon);
	fprintf(eph_results, "#-------------------------------------------\n");
	fprintf(eph_results, "#  t    ne    nf\n");
	fprintf(eph_results, "#-------------------------------------------\n");

	fprintf(behaviour, "# LASERAC V. 5\n");
	fprintf(behaviour, "#-------------------------------------------\n");
	fprintf(behaviour, "# SIDE = %d\n", SIDE);
	fprintf(behaviour, "# MAX_PHOT = %d\n", MAX_PHOT);
	fprintf(behaviour, "# NUMDIV = %d\n", numdiv);
	fprintf(behaviour, "# TMAX = %d\n", tmax);
	fprintf(behaviour, "# TTRANS = %d\n", ttrans);
	fprintf(behaviour, "# nphot_noise = %d    type(0=INIC, 1=CONT)=%d\n", nphot_noise, ntype_noise);
	fprintf(behaviour, "# threshold 1->0 = %d\n", threshold);
	fprintf(behaviour, "# SEED = %d\n", seed);
	fprintf(behaviour, "#-------------------------------------------\n");
	fprintf(behaviour, "# PB   lifet_e   lifet_f   Se    Sf    Ne_mean    Nf_mean   Events.MAX_PHOT \n");
	fprintf(behaviour, "#-----------------------------------------------------------------------------------------------------------\n");

	
	fprintf(timing_results, "# TMAX = %d, SIDE = %d, sizeof(lifet_f) = %zd, sizeof(e) = %zd, sizeof(f1) = %zd,  ",
		TMAX , SIDE, sizeof(lifet_f[0][0][0]), sizeof(e[0][0]), sizeof(f1[0][0]) );
#ifdef RANDOM_VECT_ENABLE
	fprintf(timing_results, " RANDOM_VECT_ENABLE = YES");
#endif
	fprintf(timing_results, "\n# threads\t time \n");
	fprintf(timing_results, "#-----------------------------------------------------------------------------------------------------------\n");

	//-----------------------------------------------------------------------------
	int num_thr=1;
#ifdef ENABLE_OMP  
	for (num_thr = MIN_NUM_THREADS; num_thr <= MAX_NUM_THREADS; num_thr++)
#endif
	{
#ifdef ENABLE_OMP  
		omp_set_num_threads(num_thr); 
		printf("\n# Num Threads: %d  \n", num_thr);
		fprintf(eph_results, "\n# NOF THREADS %d  \n#----------------------------------------------------------\n", num_thr);
		fprintf(behaviour, "\n# NOF THREADS %d  \n#----------------------------------------------------------\n", num_thr);
#endif
			// Using PCG Random Number Generation for C.
			// (Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>)
			// Using a local Random Number Generator
			// Seed with a fixed constant introduced in the parameters file

		// seed for the photon noise (photon creation )
		pcg32_srandom_r(&rng_PhotonNoise, seed, seed + 2);

#ifdef ENABLE_OMP  
#pragma omp parallel for default (none) shared (rng_Pumping, seed)
#endif
		// A seed for each cell. Used in the Pumping 
		// according to PGC:
		//     Seed the rng.  Specified in two parts, state initializer and a
		//     sequence selection constant (a.k.a. stream id)
		//		Controls which RNG sequence (stream) is
		//		selected. Must *always* be odd.

		for (int ii = 0; ii < SIDE; ii++) {
			// Seed with a fixed constant introduced in the parameters file and the cell linear address 
			for (int jj = 0; jj < SIDE; jj++) {
				int odd_stream_seed = (2 * ii * SIDE) + (2 * jj + 1);
				pcg32_srandom_r(&(rng_Pumping[ii][jj]), seed, odd_stream_seed);
			}
		}
			   		 
		box_width_e = SIDE * SIDE / numdiv;
		box_width_f = MAX_PHOT * SIDE * SIDE / numdiv;

#ifdef OUTER_RULE_TIMING
		double t0_total, t_total;
		double t0_param_previo, t_sum_param_previo = 0.0;

		double t_min_step, t_max_step; 

		t_min_step = 1.0e+15;
		t_max_step = 1.0e-15;
#endif
		t_sum_all_steps = 0.0; ;
		t_sum_all_param = 0.0; ;

		nof_executions_param = 0;;

		for (ppumping = ppumpingmin; ppumping < ppumpingmax + 0.00001; ppumping += ppumping_step) {
#ifdef RANDOM_VECT_ENABLE 
		//@ generating random_vect	
			uint32_t ppumping_integer_PCG_MAX = ppumping * ((double)(1ULL << 32ULL) - 1); // PCG RAND_MAX is (2^32)-1
			for (int tt = 0; tt < TMAX; tt++) { //@ big TMAX : enough memory 
// in order to have the same number for different numbers of threads, this CANNOT BE Executed in parallel: NO #pragma omp parallel for 
				for (int ii = 0; ii < SIDE; ii++) { 
					// in order to have the same number for different numbers of threads, this CANNOT BE Executed in parallel:  the same seed is always used 
					for (int jj = 0; jj < SIDE; jj++) { 
						uint64_t pp_rand_PCG = pcg32_random_r(&(rng_Pumping[ii][jj]));
						random_vect[tt][ii][jj] = (pp_rand_PCG < ppumping_integer_PCG_MAX);
					}
				}
			}
#endif
			for (lifet_electron = lifet_electronmin; lifet_electron <= lifet_electronmax; lifet_electron += lifet_electron_step) {
				for (lifet_photon = lifet_photonmin; lifet_photon <= lifet_photonmax; lifet_photon += lifet_photon_step)
				{
					nof_executions_param++;

					// printf(" PARAMETERS: PB , lifet_E, lifet_F :  %8.4f  %d  %d\n", ppumping, lifet_electron, lifet_photon);
					fprintf(eph_results, "PARAMETERS: PB , lifet_E, lifet_F :  %8.4f  %d  %d\n", ppumping, lifet_electron, lifet_photon);

					// INITIALIZATION 
					events_max_phot = 0;
					numele = 0;
					num_phot = 0;
					entropye = 0.0;
					entropyf = 0.0;
					nemax = 0;
					nemin = SIDE * SIDE;
					nfmax = 0;
					nfmin = MAX_PHOT * nemin;
					sume = 0;
					sumf = 0;
					sumnume_box_ = 0;
					sumnumf_box_ = 0;
					// nphot_noise = nphot_noise * (SIDE / 256)*(SIDE / 256);  //proportional to SIDE

#ifdef OUTER_RULE_TIMING
					t0_param_previo = t0_total = omp_get_wtime();
#endif
#ifdef ENABLE_OMP  
#pragma omp parallel for default (none) shared (e, f1, f2, lifet_f)
#endif
					for (int i = 0; i < SIDE; i++) {
						for (int j = 0; j < SIDE; j++) {  
							e[i][j] = 0;
							f1[i][j] = f2[i][j] = 0;
							for (int k = 0; k <= MAX_PHOT; k++) {  
								lifet_f[k][i][j] = 0; //@optim_lifet_f
							}
						}
					}
					
#ifdef ENABLE_OMP  
#pragma omp parallel for default (none) shared (nume_box_, numdiv, numf_box_)
#endif
					for (int n_box_ = 0; n_box_ < numdiv; n_box_++) {
						nume_box_[n_box_] = numf_box_[n_box_] = 0;
					}

#ifdef OUTER_RULE_TIMING
					t_sum_param_previo += omp_get_wtime() - t0_param_previo;  
					double t0, t_inc_Decay_v5, t_inc_Pumping_v5, t_inc_PhotonNoise_v5,
						funct_t_min_Decay_v5 = 1.0e+15, funct_t_min_Pumping_v5 = 1.0e+15, funct_t_min_PhotonNoise_v5 = 1.0e+15; 
					double funct_t_max_Decay_v5 = 1.0e-15, funct_t_max_Pumping_v5 = 1.0e-15, funct_t_max_PhotonNoise_v5 = 1.0e-15; 
					double funct_t_sum_Decay_v5 = 0, funct_t_sum_Pumping_v5 = 0, funct_t_sum_PhotonNoise_v5 = 0; 
#endif
#ifdef 	PRINT_ENABLE
					printf("Time init vectors = %lf \n", t_sum_param_previo); 
#endif

					tmax = TMAX;
					// SIMULATION TIMING LOOP 
					for (tstep = 0; tstep < tmax; ++tstep)
					{
						t0_step = omp_get_wtime();
#ifdef OUTER_RULE_TIMING
						t0 = t0_step;
#endif
						PhotonNoise_v5();

#ifdef OUTER_RULE_TIMING
						t_inc_PhotonNoise_v5 = omp_get_wtime() - t0;
						if (t_inc_PhotonNoise_v5 < funct_t_min_PhotonNoise_v5)
							funct_t_min_PhotonNoise_v5 = t_inc_PhotonNoise_v5; 
						if (t_inc_PhotonNoise_v5 > funct_t_max_PhotonNoise_v5)
							funct_t_max_PhotonNoise_v5 = t_inc_PhotonNoise_v5; 
						funct_t_sum_PhotonNoise_v5 += t_inc_PhotonNoise_v5; 

						t0 = omp_get_wtime();  
#endif
						Decay_v5();

#ifdef OUTER_RULE_TIMING
						t_inc_Decay_v5 = omp_get_wtime() - t0;
						if (t_inc_Decay_v5 < funct_t_min_Decay_v5) funct_t_min_Decay_v5 = t_inc_Decay_v5; 
						if (t_inc_Decay_v5 > funct_t_max_Decay_v5) funct_t_max_Decay_v5 = t_inc_Decay_v5; 
						funct_t_sum_Decay_v5 += t_inc_Decay_v5; 

						t0 = omp_get_wtime();
#endif
						Pumping_v5();

#ifdef OUTER_RULE_TIMING
						t_inc_Pumping_v5 = omp_get_wtime() - t0;
						if (t_inc_Pumping_v5 < funct_t_min_Pumping_v5) funct_t_min_Pumping_v5 = t_inc_Pumping_v5; 
						if (t_inc_Pumping_v5 > funct_t_max_Pumping_v5) funct_t_max_Pumping_v5 = t_inc_Pumping_v5; 
						funct_t_sum_Pumping_v5 += t_inc_Pumping_v5; 
#endif
						// SHANNON ENTROPY COMPTUATION 
						numele = 0;
						num_phot = 0;

#ifdef ENABLE_OMP
#pragma omp parallel for default (none) shared (e, pf1) reduction (+:numele ) reduction (+:num_phot )
#endif
						for (int i = 0; i < SIDE; i++) {
							for (int j = 0; j < SIDE; j++) {  
								numele = numele + (e[i][j] > 1);
								num_phot = num_phot + pf1[i][j];
							}
						}

						// COMPUTATION AFTERE TRANSIENT TIME 
						if (tstep > ttrans) {
							if (numele < nemin) nemin = numele;
							if (numele > nemax) nemax = numele;
							if (num_phot < nfmin) nfmin = num_phot;
							if (num_phot > nfmax) nfmax = num_phot;

							sume += numele;
							sumf += num_phot;

							// events in each _box_ for Shannon entropy 
							muda1e = modf(numele / box_width_e, &rn_box_e);
							n_box_e = (int)rn_box_e;
							nume_box_[n_box_e] ++;

							muda1f = modf(num_phot / box_width_f, &rn_box_f);
							n_box_f = (int)rn_box_f;
							numf_box_[n_box_f] ++;
						}

						// Total Execution time 
						t_inc_step = omp_get_wtime() - t0_step;
						t_sum_all_steps += t_inc_step; 
#ifdef OUTER_RULE_TIMING
						if (t_inc_step < t_min_step) t_min_step = t_inc_step;
						if (t_inc_step > t_max_step) t_max_step = t_inc_step;
#endif

						fprintf(eph_results, "%d", tstep);
						fprintf(eph_results, "\t%d", numele);
						fprintf(eph_results, "\t%d\n", num_phot);
						if (numele<0 || num_phot<0 || numele >(SIDE*SIDE) || num_phot >(SIDE*SIDE*MAX_PHOT)) {
							fclose(eph_results);
							fclose(behaviour);
							fclose(input);
							perror("Error counting num_phot, numele \n");
							return ERROR;
						}

						if (ferror(eph_results)) {
							//fclose ( eph_results );
							perror("Error accessing results.txt:\n");
							return ERROR;
						}

					} // end of for (tstep = 0; tstep < tmax; ++tstep)

#ifdef OUTER_RULE_TIMING
					t_total = omp_get_wtime() - t0_total;
#endif

					// END OF SIMULATION
					fprintf(eph_results, "# events MAX_PHOT= %d\n\n", events_max_phot);
					
					// n RELATIVE RANGES ("DELTAREL")
					ne_mean = (double)sume / (tmax - ttrans);
					nf_mean = (double)sumf / (tmax - ttrans);

					// SHANNON ENTROPY 
					//fprintf(behaviour, "# n_box_   nume   numf\n");
					//fprintf(behaviour, "#-------------------\n");
					for (n_box_ = 0; n_box_ < numdiv; n_box_++) {
						//fprintf(behaviour, "%d   %d   %d\n",n_box_, nume_box_ [n_box_], numf_box_ [n_box_]);
						sumnume_box_ += nume_box_[n_box_];
						sumnumf_box_ += numf_box_[n_box_];
					}

					for (n_box_ = 0; n_box_ < numdiv; n_box_++) {
						if (nume_box_[n_box_] != 0) {
							frec = (double)nume_box_[n_box_] / sumnume_box_;
							entropye = entropye + frec * log(frec);
						}
						if (numf_box_[n_box_] != 0) {
							frec = (double)numf_box_[n_box_] / sumnumf_box_;
							entropyf = entropyf + frec * log(frec);
						}
					}
					entropye *= -1.44269504;
					entropyf *= -1.44269504;

					fprintf(behaviour, "%8.4f   %d   %d", ppumping, lifet_electron, lifet_photon);
					fprintf(behaviour, "    %9.4f  %9.4f", entropye, entropyf);
					fprintf(behaviour, "    %8.2f  %8.2f", ne_mean, nf_mean);
					fprintf(behaviour, "    %d\n", events_max_phot);

					if (ferror(behaviour)) {
						//fclose ( behaviour );
						perror("Error accessing results.txt:\n");
						return ERROR;
					}

#ifdef OUTER_RULE_TIMING
					printf("\nfunct_t_min_Decay_v5 (Decay_v5) = %lf , funct_t_min_Pumping_v5 (Pumping_v5) = %lf , funct_t_min_PhotonNoise_v5 (noise photon ) = %lf \n", funct_t_min_Decay_v5, funct_t_min_Pumping_v5, funct_t_min_PhotonNoise_v5); 
					printf("funct_t_max_Decay_v5 (Decay_v5) = %lf , funct_t_max_Pumping_v5 (Pumping_v5) = %lf , funct_t_max_PhotonNoise_v5 (noise photon ) = %lf \n", funct_t_max_Decay_v5, funct_t_max_Pumping_v5, funct_t_max_PhotonNoise_v5); 
					printf("\nfunct_t_mean (Decay_v5) = %lf , funct_t_mean (Pumping_v5) = %lf , funct_t_mean2 (noise photon ) = %lf \n", funct_t_sum_Decay_v5 / tmax, funct_t_sum_Pumping_v5 / tmax, funct_t_sum_PhotonNoise_v5 / tmax); 
#endif
#ifdef 	PRINT_ENABLE
					printf("min time per step (omp_get_wtime)    = %lf \n", t_min_step); 
					printf("max time per step (omp_get_wtime)    = %lf \n", t_max_step); 
					printf("Mean time per step (omp_get_wtime)    = %lf \n", t_sum_all_steps / tmax); 
					printf("Total time all steps (omp_get_wtime) = %lf \n", t_sum_all_steps); 

					printf("Total time all steps PLUS previo (omp_get_wtime) = %lf \n", t_sum_all_steps + t_sum_param_previo); 
					printf("Total time INCLUDING printf (omp_get_wtime) = %lf \n", t_total); 
#endif
					fprintf(timing_results, "%d\t %8.4f\n", num_thr , t_sum_all_steps);
					
					t_sum_all_param += t_sum_all_steps;

				} // END OF  for (lifet_photon = lifet_photonmin; lifet_photon <= lifet_photonmax; lifet_photon += lifet_photon_step) 
			}
}

#ifdef 	PRINT_ENABLE
printf("\nTotal time (omp_get_wtime) = %lf ; mean time per execution = %lf \n", t_sum_all_param, t_sum_all_param / nof_executions_param);
#endif


	}//	for (num_thr = 1; num_thr < 10; num_thr++)

	fprintf(timing_results, "\n");

	fclose(eph_results);
	fclose(behaviour);
	fclose(timing_results);
	fclose(input);

	double tcomplete_end = omp_get_wtime();   
	printf("\nTOTAL PROGRAM TIME including FILE READING (omp_get_wtime) = %lf \n", tcomplete_end - tcomplete_start); 

	return OK;
}

//------------------------------------------------------------------
void Decay_v5(void)
{
#ifdef ENABLE_OMP  
#pragma omp parallel for default (none) shared (pf1, pf2, lifet_f )
#endif
	for (int i = 0; i < SIDE; i++)
	{   
		for (int j = 0; j < SIDE; j++)
		{
			// PHOTON DECAY and absortion:
			// 'n' HAS MAX_PHOT AS MAXIMUN VALUE (tipical 10).
			// THIS CODE MUST BE SEQUENTIAL BECAUSE IF A PHOTON DECAYS TO 0, 
			// ITS POSITION IN THE 'lifet_f' VECTOR MUST BE OCCUPIED BY THE NEXT PHOTON: THE PHOTON VECTOR IS SHIFTED.
			// THEREFORE, THERE ARE DEPENDENCIES BETWEEN CONSECUTIVE ELEMENTS IN THE 'lifet_f' VECTOR.

			int n = pf1[i][j];
			int j2 = 0;
			int k = 0;
			while (j2 < n) {
				lifet_f[k][i][j] = lifet_f[j2][i][j] - 1;  //@optim_lifet_f
				k += (lifet_f[k][i][j] > 0) ;  //@OPTIMIZATION
				j2++;
			}
			pf1[i][j] = k;
			pf2[i][j] = pf1[i][j];
		}
	}  //end of for (i = 0; i < SIDE; i++)

#ifdef ENABLE_OMP  
#pragma omp parallel for default (none) shared (e)
#endif
	for (int i = 0; i < SIDE; i++)
			{
				for (int j = 0; j < SIDE; j++)
				{
					// ELECTRON DECAY:
					/*	if (e[i][j] > 0) {
						e[i][j] = e[i][j] - 1;
						}
					*/
					//@ OPTIMIZATION
					e[i][j] = e[i][j] - (e[i][j] != 0);
		}
	}  //end of for (i = 0; i < SIDE; i++)

}
/////////////////////////////////////////////////////////////////////////////////

void Pumping_v5(void)
{
	// PUMPING AND STIMULATED EMISSION 
#ifdef ENABLE_OMP  
#pragma omp parallel for default (none) reduction (+:events_max_phot) shared(rng_Pumping, e, lifet_photon, ppumping, lifet_electron, lifet_f, pf2, threshold)
#endif

	for (int i = 0; i < SIDE; i++)
	{  
		for (int j = 0; j < SIDE; j++)
		{
			// PUMPING
			if (e[i][j] == 0)					
			{
#ifdef RANDOM_VECT_ENABLE 
					if (random_vect[tstep][i][j])   //@
					{ 
#endif
#ifndef RANDOM_VECT_ENABLE 
					double rand_numb3;
					rand_numb3 = ldexp(pcg32_random_r(&(rng_Pumping[i][j])), -32);
					if (rand_numb3 < ppumping)
					{ 
#endif
						e[i][j] = lifet_electron;
					}
				}
				// STIMULATED EMISSION 
				else
				{
						if (neighboors(i, j) > threshold)
						{
							if (pf2[i][j] < MAX_PHOT)
							{
								lifet_f[pf2[i][j]][i][j] = lifet_photon;  //@optim_lifet_f
								pf2[i][j] = pf2[i][j] + 1;
								e[i][j] = 0;
							}
							else
							{
								events_max_phot = events_max_phot + 1;
							}
						}
				}
					
			}  // END OF  for (int j 
		}  // END OF  for (int i = 0; i < SIDE; i++)

	// SWAPPING MATRIX POINTERS
		laser_f1_t  (*paux)[SIDE];
		paux = pf1; 	pf1 = pf2;  pf2 = paux;
		}

	//------------------------------------------------------------------
	// FUNCITON FOR INTRODUCING A CONTINUOUS PHOTON NOISE LEVEL (IN RANDOM POSITIONS )
	void PhotonNoise_v5(void)
	{
		// @ no omp is introduced here because this loop is not time significant
		// @ is omp were introduced, be careful because two photons could fall in the same (x,y) (atomic add would be needed)
		for (int n = 0; n < nphot_noise; n++) {
			int x, y;
			double rand_numb1, rand_numb2, xx, yy, mudo1, mudo2;
			rand_numb1 = ldexp(pcg32_random_r(&rng_PhotonNoise), -32);
			rand_numb2 = ldexp(pcg32_random_r(&rng_PhotonNoise), -32);

			mudo1 = modf(rand_numb1*SIDE, &xx);
			mudo2 = modf(rand_numb2*SIDE, &yy);
			x = (int)xx;	// RANDOM POSITION X 
			y = (int)yy;	// RANDOM POSITION Y 

			int f1_old = pf1[x][y];
			pf1[x][y] = f1_old + 1;

			if (f1_old < MAX_PHOT) {
				// IT CREATES A NEW PHOTON IN THAT POSITION
				lifet_f[f1_old][x][y] = lifet_photon; //@optim_lifet_f
				pf2[x][y] = pf1[x][y];
			}
			else {
				pf1[x][y] = f1_old - 1;
				events_max_phot = events_max_phot + 1;
			}

		} // end of for (n = 0; n < nphot_noise; n++) {
	}

	//------------------------------------------------------------------
	// COMPUTING PHOTON SUM IN THE NEIGHBORHOOD  (MOORE) AND THE OWN CELL 
	int neighboors(int i, int j)
	{
		int a1, a2, a3, a4;
		int sum = 0;

		a1 = i - 1;
		if (a1 < 0)
			a1 = SIDE + a1;
		a2 = i + 1;
		if (a2 >= SIDE)
			a2 = a2 - SIDE;
		a3 = j - 1;
		if (a3 < 0)
			a3 = SIDE + a3;
		a4 = j + 1;
		if (a4 >= SIDE)
			a4 = a4 - SIDE;

		sum = pf1[a1][j] + pf1[a2][j] + pf1[i][a3] + pf1[i][a4];
		sum = sum + pf1[a1][a3] + pf1[a1][a4] + pf1[a2][a3] + pf1[a2][a4] + pf1[i][j];
		return sum;
	}
