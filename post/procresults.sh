#!/bin/bash
#

for folder in `ls | grep results.25`
do
 cd ${folder}
 ../post/plot1d.sh
 cd ..
done
