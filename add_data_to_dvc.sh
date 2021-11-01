#!/bin/bash

folder_name=$1
root_path=/n/data1/bch/hemonc/camargo/li/DATA
cd $root_path

input_dir=$folder_name
git add $input_dir/*.yaml
dvc add $input_dir/CARLIN
mkdir -p $input_dir/QC
cp $input_dir/fastqc_after_pear/multiqc_report.html $input_dir/QC/multiqc_report_after_pear.html
cp $input_dir/fastqc_before_pear/multiqc_report.html $input_dir/QC/multiqc_report_before_pear.html
dvc add $input_dir/QC
#dvc add $input_dir/raw_fastq # This will waste too much resource 
dvc push