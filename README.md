# Laser-Cellular-Automata
Source code for the simulations of laser dynamics

This code is related to paper:

Cagigas-Muñiz, D.; Diaz-del-Rio, F.; López-Torres, M.R.; Jiménez-Morales, F.; Guisado, J.L. 
Developing Efficient Discrete Simulations on Multicore and GPU Architectures. 
Electronics 2020, 9, 189. 

Two solutions are provided:

* CUDA: source code that can be compiled in a Linux machine with a NVIDIA compatible card. CUDA libraries and compiler are neccesary. OpenCV libraries are also needed if any video output is required. 
  - Graphics folder has the gnuplots scripts for generating the output plots and the plots for experiments with cellular automatas with 400, 512, 1024, 2048 and 4096 sides.
  - Input folder has the data input files (.dat) for different sides of the cellular automata. These files can be used as templates for running new simulations.
  - Output folder has the output pair of files for each side experiment. The comport files (.txt) provide data about the behavior of simulations. The results files (.txt) show electron and photon productions for each time step. The results files are the input files for plot generation in Graphics folder.
  - CUDA400, CUDA512, ... CUDA8192 folders have the source, makefile and input files for each experiment.
  - Makefiles need cuda_functions.cu, kernel_initv6.cu, kernelV6.cu, globals.h and main.cpp files to create the binary output file. Makefiless can be adapted to any NVIDIA architecture card. See comments inside. The header file globals.h contains the constant parameters of the simulation like for example the grid side, or the number of threads and blocks of threads used by the CUDA hardware. 
  - Scripts run_compilations.sh, run_erase_compilations.sh and run_simulations.sh automate the compilation and launch of all the experiments.

* OpenMP:source code that can be compiled in a Linux machine. OpenMP libraries and g++ compiler are neccesary.
  - Graphics has the same CUDA experiments plus 8192 side results and graphics.
  - Input files: same as CUDA. There are also .dat files for 16384, 32768, 65536 and 131072 sides.
  - Output: contains only file time_results.txt. It has the execution times obtained for sides 400 to 8192.  
  - OPENMP400, OPENMP512, ... OPENMP8192 folders have the source, makefile and input files for each experiment. Here is where result timing and simulation behaviour files are created. This is a difference with the CUDA solution.
  - Makefiles need the laseracv10_omp_optim.cpp and pcg_basic.cpp source files, and the pcg_basic.h header file. 
  - Scripts run_compilations.sh, run_erase_compilations.sh and run_simulations.sh automate the compilation and launch of all the experiments. 

* OpenMP_NOT_OPTIMIZED: source code that can be compiled in a Linux machine. OpenMP libraries and g++ compiler are neccesary. It is a previous version of OpenMP source code that is not optimized. Execution times are higher and the 8192 SIDE experiment can not be compiled due to integer overflow.
