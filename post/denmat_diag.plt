# to be used with plot2d.sh
#
# plots x=x' density

set terminal postscript enhanced eps color
set output 'OUTPUTFILE.eps'

set object 1 rectangle from screen 0, screen 0 to screen 1, screen 1 behind fc rgbcolor "white"

unset key

set title 'probability density ("diagonal"), time = TIME fm/c'

set size 0.7

set xlabel "x [fm]" 
set ylabel "density [fm^{-1}]"

#set xrange [-10:10]
set yrange [-0.025:0.5]

plot 'INPUTFILE.dat' i INDEX u 1:NUM w lines lw 2 tit "time = TIME fm/c"
#     'denmatan_x_t.dat' i INDEX u 1:NUM w li tit "analytical"
