#!/bin/bash

# pear has been installed locally. So, no need to load module here
fq_R1=$1
fq_R2=$2
pear_out=$3

mkdir -p pear_output

pear \
   -o $pear_out \
   -f $fq_R1 \
   -r $fq_R2 \
   --min-overlap 1 \
   --min-assembly-length 1 \
   --min-trim-length 1 \
   --quality-threshold 30 \
   -j 16 
       
## the -k command only matters for the discarded sequences: 
# -k  "Do not reverse and complement the reverse reads when writing
#          the unassembled and discarded reads output."
