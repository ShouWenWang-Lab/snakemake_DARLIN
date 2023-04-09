#!/bin/bash

CARLIN_dir=$1
input_dir=$2
SampleList=$3
template=$4

module load matlab

echo "Running interactive mode for merging samples and allele bank generation"
#command_str="make_allele_bank('$SampleList','$input_dir','$template')"
command_str="merge_samples('$SampleList','$input_dir','$template')"
echo $command_str
cur_dir=$(pwd)
cd $CARLIN_dir
matlab -nodisplay -nosplash -nodesktop -r "$command_str; exit"
cd $cur_dir


