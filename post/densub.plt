#set terminal postscript enhanced color eps
#set output 'densub_0NUMBER.eps'

set pm3d map
set palette rgbformulae 30,31,32
set hidden3d

set size square

set xlabel 'x_{average}'
set ylabel 'x_{relative}'
set title '{/Symbol r}_{arnau}-{/Symbol r}_{brent}, time=NUMBER0 fm/c'

splot 'densub.dat' i NUMBER matrix

