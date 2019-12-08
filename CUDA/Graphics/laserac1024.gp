# LASERAC PROGRAM: SIMULATION OF LASER DYNAMICS USING A CELLULAR AUTOMATA
#------------------------------------------------------------------------
# GRAPHICS RESULTS USING GNUPLOT

set terminal postscript eps color

set xlabel "Time Steps"
set ylabel "Population"
set key tmargin

set output "CUDA_parallel_application1024.eps"
set title "CUDA parallel application. Lattice size = 1024x1024"
plot 'results-1024.txt' using 1:2 title 'Population Inversion',  'results-1024.txt' using 1:3 title 'Laser Photons'

reset
# pause -1 "Press 'Return' to continue"
