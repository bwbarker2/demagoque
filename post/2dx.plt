# to be used with plot2d.sh , only physical matrices

set terminal postscript enhanced eps color
#set terminal gif
set output "OUTPUTFILE.eps"

set pm3d map

set xrange [-20:20]
set yrange [-20:20]
#set zrange [-0.3:0.3]
set cbrange [-0.2:0.3]

set nokey
#unset colorbox

set xlabel 'x_{average} [fm]'
set xtics -20,10,20

set ylabel 'x_{relative} [fm]'
set ytics -20,10,20

set cblabel 'density [fm^{-1}]'

set title 'density matrix, position, t=TIME fm/c'

#set size 0.5
set size ratio -1 0.5,0.5

#set palette defined (-0.3 "black", 0 "blue", 0.3 "white")

set palette rgbformulae 30,31,32

#set multiplot

splot "INPUTFILE.dat" i INDEX u 1:2:3 #every 2

#unset pm3d
#
#set xrange [-35:35]
#set yrange [-35:35]
#set size 0.7,0.7
#set origin 0.1,0.1
#
#plot x
