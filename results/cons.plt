set terminal postscript enhanced eps color
set output 'cons_abs_p.eps'

set xlabel 'time [fm/c]'
set ylabel 'potential energy per particle [MeV]'

plot 'cons_abs.dat' u 1:4 w li

