set terminal postscript enhanced eps color
set output 'results/1dsymw.eps'

set xlabel 'time [2*fm/c]'
set ylabel 'sum_i |d_i - d_{ave}| / N

plot \
     'results/diagwreudsym.dat' u 1:($2) w li tit 'wre ud', \
     'results/diagwrelrsym.dat' u 1:($2) w li tit 'wre lr', \
     'results/diagwimudsym.dat' u 1:2 w li tit 'wim ud', \
     'results/diagwimlrsym.dat' u 1:2 w li tit 'wim lr'
#     'results/diagkreudsym.dat' u 1:2 w li tit 'kre ud', \
#     'results/diagkrelrsym.dat' u 1:2 w li tit 'kre lr', \
#     'results/diagkimudsym.dat' u 1:2 w li tit 'kim ud', \
#     'results/diagkimlrsym.dat' u 1:2 w li tit 'kim lr'
#     'results/diagxreudsym.dat' u 1:2 w li tit 'xre ud', \
#     'results/diagxrelrsym.dat' u 1:($2) w li tit 'xre lr', \
#     'results/diagximudsym.dat' u 1:($2) w li tit 'xim ud', \
#     'results/diagximlrsym.dat' u 1:($2) w li tit 'xim lr'
