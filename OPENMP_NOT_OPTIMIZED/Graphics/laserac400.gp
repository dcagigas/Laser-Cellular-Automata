# PROGRAMA LASERAC: SIMULACION DE LA DINAMICA DEL
# LASER MEDIANTE UN AUTOMATA CELULAR
#--------------------------------------------------------
# GRAFICA DE RESULTADOS USANDO GNUPLOT

set terminal postscript eps color

set xlabel "Time Steps"
set ylabel "Population"
set key tmargin

set output "OMP_parallel_application400.eps"
set title "OMP parallel application. Lattice size = 400x400"
plot 'results-400.txt' using 1:2 title 'Population Inversion', 'results-400.txt' using 1:3 title 'Laser Photons'

reset
# pause -1 "Press 'Return' to continue"
