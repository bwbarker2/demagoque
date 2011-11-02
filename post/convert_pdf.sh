#!/bin/bash
#
# Converts all images of given basename and extension to given format with 
# simply ordered sequential numbers
#
# Requires imagemagick, ext and new_ext must be formats supported by imagemagick

if [ $# -ne "3" ]
then
 echo "Usage: $0 [base-filename] [old-extension] [new-extension]"
 exit -1
fi


basefname=$1

ext=$2

new_ext=$3

num=0

for fname in $(ls ${basefname}*.${ext})
do
 echo ${fname}
# echo ${num}

# rasterizes at high resolution
 convert -density 300 ${fname} ${basefname}${num}.pdf

# converts, keeping vector format
# epstopdf --outfile=${basefname}${num}.pdf ${fname}
 num=$(echo "$num+1" | bc)
done
