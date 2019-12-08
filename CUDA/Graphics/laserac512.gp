# LASERAC PROGRAM: SIMULATION OF LASER DYNAMICS USING A CELLULAR AUTOMATA
#------------------------------------------------------------------------
# GRAPHICS RESULTS USING GNUPLOT

set terminal postscript eps color

set xlabel "Time Steps"
set ylabel "Population"
set key tmargin

set output "CUDA_parallel_application512.eps"
set title "CUDA parallel application. Lattice size = 512x512"
plot 'results-512.txt' using 1:2 title 'Population Inversion',  'results-512.txt' using 1:3 title 'Laser Photons'

reset
# pause -1 "Press 'Return' to continue"
