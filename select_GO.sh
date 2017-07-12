#!/bin/bash

input_tsv_file=$1
output_file=$2
uniq_out=$3

cut -f9,14 $1 | sed -r -n -e '/E/p' > temp
sed -r -n -e '/GO/p' temp > temp2
cut -f2 temp2 > temp3
perl -p -e 's/\|/\n/g' temp3 | sort > temp4
uniq temp4 > temp5

mv temp4 $2
mv temp5 $3

rm temp*
