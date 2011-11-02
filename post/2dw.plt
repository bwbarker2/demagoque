# to be used with plot2d.sh , only physical matrices

set terminal postscript enhanced eps color
#set terminal gif
#set terminal jpeg large crop
set output "OUTPUTFILE.eps"

set pm3d
set pm3d map
set hidden3d

set xrange [10:90]
set yrange [10.4:93.6]
#set yrange [6.4:57.6]
#set yrange [20:180]
#set zrange [-0.3:0.3]
set cbrange [-0.8:0.8]

set nokey
#unset colorbox

set xlabel 'x_{average} [fm]'
set xtics ("-20" 10, "-10" 30, "0" 50, "10" 70, "20" 90)

set ylabel 'k_{average} [fm^{-1}]'
set ytics ("-5" 10.4, "-2.5" 31.2, "0" 52, "2.5" 72.8, "5" 93.6)
#set ytics ("-5" 6.4, "-2.5" 19.2, "0" 32, "2.5" 44.7, "5" 57.4)
#set ytics ("-5" 20, "-2.5" 60, "0" 100, "2.5" 140, "5" 180)

set cblabel 'density [fm^{-1}]'
set cbtics -0.8,0.4,1.2

set title 'Wigner distribution, t=TIME fm/c'

#set size 0.5
set size ratio -0.96 0.5,0.5

#set palette defined (-0.3 "black", 0 "blue", 0.3 "white")

set palette rgbformulae 30,31,32

#set multiplot

splot "INPUTFILE.dat" i INDEX matrix w li

#unset pm3d
#
#set xrange [-35:35]
#set yrange [-35:35]
#set size 0.7,0.7
#set origin 0.1,0.1
#
#plot x
