#!/bin/bash
#
# plots 1d diagonal time evolution in x,k space

path_to_results="./"
path_from_results="../post/"

# test file to get timestep information
filetest=${path_to_results}denmat_x_t

filelist="denmat_x_t denmat_k_t"

# get number of timesteps in file
nindex=`grep time ${filetest}.dat | wc -l`
echo $nindex
# get timestep information
for ((index=1; index <= nindex; index++))
do
  times[index]=`grep time ${filetest}.dat \
               | head -n${index} | tail -n1 \
               | awk '{print $3}' \
               | sed 's/\.\([0-9]\{3\}\)0*/\_\1/g'`
done

# find maximum character length for all timesteps
let "nin0 = nindex - 1"
maxlen=${#times[nindex]}
    
# pad times with zeroes for correct filename scheme
for ((index=1; index <= nindex; index++))
  do
  while [ ${#times[index]} -lt $maxlen ]
    do
    times[index]=0${times[index]}
  done
  
  echo index = $index, time = ${times[index]}
done

# iterate over files in list(x,k)
for i in $filelist ;  do

  # iterate over real, imaginary columns in data
  for ((reim=2; reim <= 3; reim++)) ; do


    # set column to read from in data file
    if [ "$reim" = "2" ]; then
      fili="${i}_re"
    else
      fili="${i}_im"
    fi

    cd ${path_to_results}

    # iterate over time indices
    for ((index=1; index <= nindex; index++)) ; do

#      echo time = ${times[index]}

      # make index formatted for gnuplot INDEX
      let "ind1 = index - 1"

      # create subdirectory if none exists
      mkdir -p ${fili}
      
#      echo ${fili}

      sed -e "s/INPUTFILE/$i/g" \
          -e "s/OUTPUTFILE/${fili}\/${fili}${times[index]}/g" \
          -e "s/INDEX/$ind1/g" \
          -e "s/TIME/${times[index]}/g" \
	  -e "s/NUM/$reim/g" \
	  ${path_from_results}denmat_diag.plt | gnuplot
    done #time
  done #real,imag
done #files

