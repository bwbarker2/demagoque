# to be used with plot2d.sh , only physical matrices

set terminal postscript enhanced eps color
#set terminal gif
set output "OUTPUTFILE.eps"

set pm3d map

set xrange [0:100]
set yrange [0:200]
#set zrange [-0.3:0.3]
set cbrange [-0.2:0.2]

#set nokey
#unset colorbox

set size ratio -1

#set palette defined (-0.3 "black", 0 "blue", 0.3 "white")

set palette rgbformulae 30,31,32

#set multiplot

splot "INPUTFILE.dat" i INDEX matrix

#unset pm3d
#
#set xrange [-35:35]
#set yrange [-35:35]
#set size 0.7,0.7
#set origin 0.1,0.1
#
#plot x
