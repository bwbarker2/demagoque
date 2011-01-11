set terminal postscript enhanced eps color
set output "plotcons.eps"

set xlabel "time [fm/c]"
set ylabel "energy [MeV]"

plot './results/cons_abs.dat' u 1:(398.3484*$2+$4-113.564):(398.3484*$3+$5) w yerr tit "total energy with correction to kinetic" lt 1, \
     ''  u 1:(398.3484*$2+$4-113.564) w li notit lt 1
