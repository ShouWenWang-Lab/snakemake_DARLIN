#!/bin/bash

folder_name=$1
#root_folder_name=$2
root_path=/n/data1/bch/hemonc/camargo/li
input_dir=$root_path/DATA/$folder_name
output_dir=$root_path/temp_data/$folder_name
mkdir -p $output_dir
rsync -avP $input_dir/raw_fastq $output_dir
cp $input_dir/fastqc_after_pear/multiqc_report.html $output_dir/multiqc_report_after_pear.html
cp $input_dir/fastqc_before_pear/multiqc_report.html $output_dir/multiqc_report_before_pear.html
cp $input_dir/*.yaml $output_dir
#cp $root_path/$root_folder_name/*.json $root_path/o2_data/$root_folder_name
