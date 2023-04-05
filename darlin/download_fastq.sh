#!/bin/bash
project_name=$1
project_ID=$2
data_dir=$3
mkdir -p $data_dir
bs download project $project_name --id $project_ID -o $data_dir
mkdir -p $data_dir/raw_fastq
rm -r -f $data_dir/raw_fastq/*
mv $data_dir/*/*.fastq.gz $data_dir/raw_fastq
rm -r -f $data_dir/*_ds.*
