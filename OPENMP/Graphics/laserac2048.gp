# PROGRAMA LASERAC: SIMULACION DE LA DINAMICA DEL
# LASER MEDIANTE UN AUTOMATA CELULAR
#--------------------------------------------------------
# GRAFICA DE RESULTADOS USANDO GNUPLOT

set terminal postscript eps color

set xlabel "Time Steps"
set ylabel "Population"
set key tmargin


set output "OpenMP_parallel_application2048.eps"
set title "OpenMP application. Lattice size = 2048x2048"
plot 'results-2048.txt' using 1:2 title 'Population Inversion',  'results-2048.txt' using 1:3 title 'Laser Photons'

reset
# pause -1 "Pulsa 'Return' para seguir"
