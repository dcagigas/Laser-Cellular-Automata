# LASERAC PROGRAM: SIMULATION OF LASER DYNAMICS USING A CELLULAR AUTOMATA
#------------------------------------------------------------------------
# GRAPHICS RESULTS USING GNUPLOT

set terminal postscript eps color

set xlabel "Time Steps"
set ylabel "Population"
set key tmargin

set output "CUDA_parallel_application4096.eps"
set title "CUDA parallel application. Lattice size = 4096x4096"
plot 'results-4096.txt' using 1:2 title 'Population Inversion',  'results-4096.txt' using 1:3 title 'Laser Photons'

reset
# pause -1 "Press 'Return' to continue"
