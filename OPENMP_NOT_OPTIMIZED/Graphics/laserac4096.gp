# PROGRAMA LASERAC: SIMULACION DE LA DINAMICA DEL
# LASER MEDIANTE UN AUTOMATA CELULAR
#--------------------------------------------------------
# GRAFICA DE RESULTADOS USANDO GNUPLOT

set terminal postscript eps color

set xlabel "Time Steps"
set ylabel "Population"
set key tmargin


set output "OpenMP_parallel_application4096.eps"
set title "OpenMP application. Lattice size = 4096x4096"
plot 'results-4096.txt' using 1:2 title 'Population Inversion',  'results-4096.txt' using 1:3 title 'Laser Photons'

reset
# pause -1 "Pulsa 'Return' para seguir"
