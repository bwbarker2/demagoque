#!/bin/bash
#
# plots 2d density matrices

# test file to get timestep information
filetest=2dxre

filelistx="2dxre 2dxim"
filelistk="2dkre 2dkim"

# get number of timesteps in file
nindex=`grep time ${filetest}.dat | wc -l`

# get timestep information
for ((index=1; index <= nindex; index++))
  do
  times[index]=`grep time ${filetest}.dat \
               | head -n${index} | tail -n1 \
               | awk '{print $3}' \
               | sed 's/\.\([0-9]\{3\}\)0*/\_\1/'`
done

# find maximum character length for all timesteps
let "nin0 = nindex - 1"
maxlen=${#times[nin0]}
    
# pad times with zeroes for correct filename scheme
for ((index=0; index < nindex; index++))
  do
  while [ ${#times[index]} -lt $maxlen ]
    do
    times[index]=0${times[index]}
  done
  
  echo index = $index, time = ${times[index]}
done

for ((index=0; index < nindex; index++)) ; do

    echo time = ${times[index]}
  
  for i in $filelistx ;  do
      
      # create subdirectory if none exists
      mkdir -p $i
      
      echo $i

    sed -e "s/INPUTFILE/$i/" -e "s/OUTPUTFILE/${i}\/${i}${times[index]}/" -e "s/INDEX/$index/" \
	2dx.plt | gnuplot
  done

  for i in $filelistk ;  do
      
      # create subdirectory if none exists
      mkdir -p $i

      echo $i
      
    sed -e "s/INPUTFILE/$i/" -e "s/OUTPUTFILE/${i}\/${i}${times[index]}/" -e "s/INDEX/$index/" \
	2dk.plt | gnuplot
  done

done

