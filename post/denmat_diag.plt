# to be used with plot2d.sh
#
# plots x=x' density

set terminal postscript enhanced eps color
set output 'OUTPUTFILE.eps'

set object 1 rectangle from screen 0, screen 0 to screen 1, screen 1 behind fc rgbcolor "white"

#set xrange [-10:10]
set yrange [-0.025:0.125]

plot 'INPUTFILE.dat' i INDEX u 1:NUM w lines tit "time = TIME"
#     'denmatan_x_t.dat' i INDEX u 1:NUM w li tit "analytical"
