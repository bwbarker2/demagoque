# to be used with plot2d.sh , only spectral matrices

#set terminal gif
set terminal postscript enhanced eps color
set output "OUTPUTFILE.eps"

set pm3d map

#set xrange [0:100]
#set yrange [0:100]
#set zrange [-0.3:0.3]
#set cbrange [-0.001:0.001]

set size ratio -1

#set palette defined (-0.3 "black", 0 "blue", 0.3 "white")

set palette rgbformulae 30,31,32

splot "INPUTFILE.dat" i INDEX matrix
