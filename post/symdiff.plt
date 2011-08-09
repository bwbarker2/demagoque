#plot 'diagxreudsym.dat' u 1:($2/5e-9) w li, \
#     'diagxrelrsym.dat' u 1:($2/2.5e-10) w li, \
#     'diagximudsym.dat' u 1:($2/1.4e-16) w li, \
#     'diagximlrsym.dat' u 1:($2/0.008) w li, \
#     'diagwreudsym.dat' u 1:($2/0.0035) w li, \
#     'diagwrelrsym.dat' u 1:($2/0.0035) w li
plot 'diagwimudsym.dat' u 1:($2/9e-10) w li, \
     'diagwimlrsym.dat' u 1:($2/7e-17) w li, \
     'diagkreudsym.dat' u 1:($2/1e-16) w li, \
     'diagkrelrsym.dat' u 1:($2/9e-17) w li, \
     'diagkimudsym.dat' u 1:($2/0.006) w li, \
     'diagkimlrsym.dat' u 1:($2/0.006) w li
