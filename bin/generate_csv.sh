#!/bin/bash

CARLIN_dir=$1
input_dir=$2
SampleList=$3
template=$4


module load matlab

echo "Running interactive mode for csv generation"
command_str="csv_reports('$SampleList','$input_dir','$template')"
echo $command_str
cur_dir=$(pwd)
cd $CARLIN_dir
matlab -nodisplay -nosplash -nodesktop -r "$command_str; exit"
cd $cur_dir



