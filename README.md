# Laser-Cellular-Automata
Source code for the simulations of laser dynamics

Two solutions are provided:

* CUDA: source code that can be compiled in a Linux machine with a NVIDIA compatible card.
  - Graphics folder has the gnuplots scripts for generating the output plots and the plots for experiments with cellular automatas with 400, 512, 1024, 2048 and 4096 sides.
  - Input folder has the data input files (.dat) for different sides of the cellular automata. These files can be used as templates for running new simulations.
  - Output folder has the output pair of files for each side experiment. The comport files (.txt) provide data about the simulations. The results files (.txt) show electron and photon productions for each time step. The results files are the input files for plot generation in Graphics folder.
  - Makefile need cuda_functions.cu, kernel_initv6.cu, kernelV6.cu and main.cpp to create the binary output file. This Makefile can be adapted to the NVIDIA architecture of your NVIDIA card. See comments inside.

* OpenMP:source code that can be compiled in a Linux machine. OpenMP libraries and g++ compiler are neccesary.
