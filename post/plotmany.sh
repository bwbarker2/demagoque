#!/bin/bash
#
# plots 2d density matrices
# WARNING: Incomplete, not ready to run yet.

path_to_results="./"
path_from_results="../post/"

if [$# -ne 1]
then
 echo "Usage: $0 [basefname]"
 exit -1
fi

outmode=$1

case "${outmode}" in
 2dxre)		pltfile="2dx.plt"
		basefname="2dxre"
		options="matrix"
		;;
 2dx_re)	pltfile="2dx.plt"
                basefname="2dx"
		options="u 1:2:3"
		;;
 *)		echo "plot2d.sh: Not valid outmode: ${outmode}"
		;;
esac

# test file to get timestep information
filetest="${path_to_results}/2dxre"

#filelistx="2dxre 2dxim"
filelistx="2dx"
#filelistk="2dkre 2dkim"
filelistk=""  # only printing k right now

# get number of timesteps in file
nindex=`grep 'time' ${basefname}.dat | wc -l`

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
for ((index=1; index <= nindex; index++))
  do
  while [ ${#times[index]} -lt $maxlen ]
    do
    times[index]=0${times[index]}
  done
  
  echo index = $index, time = ${times[index]}
done

for ((index=1; index <= nindex; index++)) ; do

    echo time = ${times[index]}
  
  cd ${path_to_results}

  for i in $filelistx ;  do
      
      # create subdirectory if none exists
      mkdir -p $i
      
      echo $i

      i=${i}

    let "ind0 = index - 1"
    sed -e "s/INPUTFILE/$i/" -e "s/OUTPUTFILE/${i}\/${i}${times[index]}/" -e "s/INDEX/$ind0/" \
	${path_from_results}/2dx.plt | gnuplot
  done

  for i in $filelistk ;  do
      
      # create subdirectory if none exists
      mkdir -p $i

      echo $i
      
      i=${i}
      let "ind0 = index - 1"
      sed -e "s/INPUTFILE/$i/" -e "s/OUTPUTFILE/${i}\/${i}${times[index]}/" -e "s/INDEX/$ind0/" \
	${path_from_results}/2dk.plt | gnuplot
  done

done

