// CONSTANT DEFINITION:
#define SIDE 8192	// MATRIX SIDE
#define MAXFOT 10   // MAX NUMBER OF PHOTONS IN EACH CELL (X,Y)
#define OK 0        // RETURNED VALUE IF NO ERRORS
#define ERROR 1     // RETURNED VALUE IF ERRORS
#define TMAX 1000

#define num_threads 256
#define num_blocks (SIDE*SIDE/num_threads)

                                // Change data size if convenient:
#define CELL_TYPE int           // for GPU_f1 and GPU_f2 matrixes
#define CELL_TYPE2 short int    // for GPU_e matrix
#define CELL_TYPE3 short int    // for GPU_tvf matrix

                        // GPU memory allocated:
                        // 1) Buffers and variables for Shannon Entropy:
                        // GPU_nemin + GPU_nemax + GPU_nfmin + GPU_nfmax + GPU_sumae + GPU_sumaf + GPU_sumanumecaja + GPU_sumanumfcaja + 
                        // GPU_rncajae + GPU_rncajaf +
                        // GPU_numecaja + GPU_numfcaja + 
                        // GPU_entropiae + GPU_entropiaf
                        // = 8*sizeof(int) + 2*sizeof(double) + 2*(parametros.numdiv)*sizeof(int) + (parametros.numdiv)*sizeof(double)
                        // 2) Main GPU buffers:
                        // GPU_e + GPU_f1 + GPU_f2 + GPU_tvf + GPU_total_e + GPU_total_f + GPU_ocurrencias_MAXFOT
                        // = SIDE*SIDE*sizeof(CELL_TYPE2) + 2*SIDE*SIDE*sizeof(CELL_TYPE) + SIDE*SIDE*MAXFOT*sizeof(CELL_TYPE2) + 
                        // 2*TMAX*sizeof(unsigned int) + sizeof(int)
#define TOTAL_MEMORY (8*sizeof(int) + 2*sizeof(double) + 2*TMAX*sizeof(int) + TMAX*sizeof(double)) + (SIDE*SIDE*sizeof(CELL_TYPE2)+ 2*SIDE*SIDE*sizeof(CELL_TYPE) + SIDE*SIDE*MAXFOT*sizeof(CELL_TYPE3) + 2*TMAX*sizeof(unsigned int) + sizeof(int))


typedef struct {
	float pbombeo, pbombeomin, pbombeomax, pbombeopaso;
	CELL_TYPE tvelectronmin, tvelectronmax, tvelectronpaso;
	CELL_TYPE tvfotonmin, tvfotonmax, tvfotonpaso;
	CELL_TYPE ntiporuido;
	CELL_TYPE numdiv, ttrans;
	CELL_TYPE numele, numfot;
	CELL_TYPE tvelectron, tvfoton;
	CELL_TYPE tmax, nfotruido, threshold;
	int occurrencesMAXFOT;		
	int seed;
} SIMUL_DATA;

int Launcher_init_v6 ();
int Launcher_v6(SIMUL_DATA * parametros, unsigned int * total_e, unsigned int * total_f);

int ShowGPUInfo(unsigned int mem);
void PrintgGrid(void);


