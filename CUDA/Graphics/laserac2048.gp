# LASERAC PROGRAM: SIMULATION OF LASER DYNAMICS USING A CELLULAR AUTOMATA
#------------------------------------------------------------------------
# GRAPHICS RESULTS USING GNUPLOT

set terminal postscript eps color

set xlabel "Time Steps"
set ylabel "Population"
set key tmargin

set output "CUDA_parallel_application2048.eps"
set title "CUDA parallel application. Lattice size = 2048x2048"
plot 'results-2048.txt' using 1:2 title 'Population Inversion',  'results-2048.txt' using 1:3 title 'Laser Photons'

reset
# pause -1 "Press 'Return' to continue"
