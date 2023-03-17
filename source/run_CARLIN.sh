#!/bin/bash

CARLIN_dir=$1
input_dir=$2
output_dir=$3
sample=$4
cfg_type=$5
template=$6
read_cutoff_UMI_override=$7
read_cutoff_CB_override=$8
requested_memory=$9
sbatch_job=${10}
max_run_time=${11}

# # This is used to test the pipeline
# CARLIN_dir='/n/data1/bch/hemonc/camargo/li/CARLIN_pipeline/Custom_CARLIN'
# input_dir='/n/data1/bch/hemonc/camargo/li/CARLIN_pipeline/Custom_CARLIN/tests/Tigre_Carlin_test'
# output_dir='/n/data1/bch/hemonc/camargo/li/CARLIN_pipeline/Custom_CARLIN/tests/Tigre_Carlin_test/output2'
# sample='test_Tigre_Carlin'
# cfg_type='BulkDNA_Tigre'
# template='Tigre'
# read_cutoff_CB_override=3
# read_cutoff_UMI_override=8
# requested_memory=1024 # the unit here is M. By default, we request 1G
# sbatch_job=1 # 1 to run sbatch, 0 for not
# max_run_time=1

#echo sbatch_job=${sbatch_job}
# file_size_command="$input_dir/$sample.trimmed.pear.assembled.fastq";file_size="$(du $file_size_command)";requested_memory="$($CARLIN_memory_factor*$file_size)"

module load matlab
mkdir -p log

if [[ $sbatch_job -eq 0 ]]
then
    echo "Running interactive mode for CARLIN analysis"
command_str="my_CARLIN_pipeline('$sample','$cfg_type','$input_dir','$output_dir','$template','read_cutoff_UMI_override',$read_cutoff_UMI_override,'read_cutoff_CB_override',$read_cutoff_CB_override,'CARLIN_dir','$CARLIN_dir')"
    #echo $command_str
    cur_dir=$(pwd)
    cd $CARLIN_dir
    matlab -nodisplay -nosplash -nodesktop -r "$command_str; exit"
    cd $cur_dir
else
    echo "Running batch jobs for CARLIN analysis"
    sbatch -p intel-sc3 -c 6 -t $max_run_time:00:00 --mem=${requested_memory}G --job-name $sample --output=log/${sample}-%j.o  --error=log/${sample}-%j.e --mail-type=TIME_LIMIT_90,FAIL,END --wrap="cd ${CARLIN_dir}; matlab -batch \"my_CARLIN_pipeline(\\\"${sample}\\\", \\\"${cfg_type}\\\", \\\"${input_dir}\\\", \\\"${output_dir}\\\", \\\"${template}\\\", \\\"read_cutoff_CB_override\\\", $read_cutoff_CB_override, \\\"read_cutoff_UMI_override\\\", $read_cutoff_UMI_override, \\\"CARLIN_dir\\\", \\\"${CARLIN_dir}\\\")\""
    
fi
