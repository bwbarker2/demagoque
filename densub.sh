#!/bin/bash

for i in `seq 0 4`;
do
 sed "s/NUMBER/$i/" ../post/densub.plt | gnuplot
done
