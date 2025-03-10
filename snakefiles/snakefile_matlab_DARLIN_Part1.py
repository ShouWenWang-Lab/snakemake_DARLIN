import os
import sys
from pathlib import Path

import pandas as pd

from darlin import help_functions as hf
from darlin.settings import script_dir, QC_dir, ref_dir, CARLIN_dir
#configfile: "config.yaml"  # command line way to set it: --configfile 'path/to/config'
#workdir: config['data_dir'] # set working directory, a command-line way to set it: --directory 'path/to/your/dir'
config['data_dir']=str(os.getcwd())

##################
## preprocessing
################## 
CARLIN_dir=hf.update_CARLIN_dir(CARLIN_dir,config['template'])
print("Updated CARLIN_dir:"+ str(CARLIN_dir))

cfg_type=config['cfg_type']
if 'sbatch_mode' in config.keys():
    sbatch_mode=config['sbatch_mode']
else:
    sbatch_mode='intel-sc3'

if config['template'] == 'Tigre':
    print("------------Warn: remember that the Tigre template is inversed-------------")

if len(config['SampleList'])==0: 
    df=pd.read_csv('raw_fastq/sample_info.csv')
    SampleList=list(df['Sample_long'])
else:
    SampleList=config['SampleList']

    
# this is to make the pipeline compatible with earlier bulk config files
if cfg_type.startswith('Bulk') and ('read_cutoff_UMI_override' not in config.keys()) and ('read_cutoff_override' in config.keys()):
    config['read_cutoff_UMI_override']=config['read_cutoff_override']
    
DARLIN_sub_dir=[f"results_cutoff_override_{xx}" for xx in config['read_cutoff_UMI_override']]

    

# remove the flag file of the workflow if the sbatch is not actually run to finish
for sample in SampleList:
    if not os.path.exists(f'DARLIN/{DARLIN_sub_dir}/{sample}/DARLIN_analysis_actually.done'):
        if os.path.exists(f'DARLIN/{DARLIN_sub_dir}/{sample}/DARLIN_analysis.done'):
            os.remove(f'DARLIN/{DARLIN_sub_dir}/{sample}/DARLIN_analysis.done') 
        
##################
## start the rules
################## 
rule all:
    input: 
        expand("DARLIN/{sub_dir}/{sample}/DARLIN_analysis.done",sample=SampleList,sub_dir=DARLIN_sub_dir)

rule DARLIN:        
    input:
        fq_R1="raw_fastq/{sample}_R1.fastq.gz",
        fq_R2="raw_fastq/{sample}_R2.fastq.gz"
    output:
        touch("DARLIN/{sub_dir}/{sample}/DARLIN_analysis.done")
    run:
        output_dir=config['data_dir']+f'/DARLIN/{wildcards.sub_dir}'
        
        if cfg_type.startswith('Bulk'):
            input_dir=config['data_dir']+'/pear_output'
            if (not os.path.exists(f'{output_dir}/{wildcards.sample}/Summary.mat')):
                # part 1: run pear
                out_dir=f"pear_output"
                os.makedirs(out_dir,exist_ok=True)
                command_1=f"bash {script_dir}/run_pear.sh {input.fq_R1} {input.fq_R2} pear_output/{wildcards.sample}.trimmed.pear"
                command_2=f"bash {script_dir}/run_fastqc.sh {input.fq_R1} fastqc_before_pear; sh {script_dir}/run_fastqc.sh {input.fq_R2} fastqc_before_pear"
                command_3=f"bash {script_dir}/run_fastqc.sh pear_output/{wildcards.sample}.trimmed.pear.assembled.fastq  fastqc_after_pear"
                #command_4=f"bash {script_dir}/run_multiqc.sh fastqc_before_pear; sh {script_dir}/run_multiqc.sh fastqc_after_pear"
            else:
                command_1="echo command_1"
                command_2="echo command_2"
                command_3="echo command_3"
                
            
        elif cfg_type.startswith('sc'):
            input_dir=config['data_dir']+'/raw_fastq'
            if (not os.path.exists(f'{output_dir}/{wildcards.sample}/Summary.mat')):
                command_1=f"bash {script_dir}/run_fastqc.sh {input.fq_R1} fastqc_before_pear; sh {script_dir}/run_fastqc.sh {input.fq_R2} fastqc_before_pear"
                command_2="echo command_2"
                command_3="echo command_3"
                #command_4="echo command_2"
            else:
                command_1="echo command_1"
                command_2="echo command_2"
                command_3="echo command_3"

            
        print(input_dir)
        template=config['template']
        
        CARLIN_memory_factor=config['CARLIN_memory_factor']
        sbatch=config['sbatch']
        CARLIN_max_run_time=config['CARLIN_max_run_time']
        read_cutoff_UMI_override=int(wildcards.sub_dir.split('_')[-1])
        read_cutoff_CB_override=read_cutoff_UMI_override
        
        file_size = os.path.getsize(f'{input.fq_R1}')*5/1000000000
        print(f"{wildcards.sample}:   FileSize {file_size} G")
        requested_memory=int(file_size*CARLIN_memory_factor)
        if requested_memory<20:
            requested_memory=20 # at least request 20G memory
        if requested_memory>250:
            requested_memory=250 # do not request more than 200G memory
        print(f"{wildcards.sample}:   Requested memory {requested_memory} G")
        os.makedirs(f'{output_dir}/{wildcards.sample}',exist_ok=True)
        
        # do not run sbatch within this command
        command_4=f"bash {script_dir}/run_CARLIN.sh {CARLIN_dir} {input_dir} {output_dir} {wildcards.sample} {cfg_type} {template} {read_cutoff_UMI_override} {read_cutoff_CB_override} {requested_memory} {0} {CARLIN_max_run_time}"

        combined_command=f"{command_1}; {command_2}; {command_3}; {command_4}"
        
        
        job_name=f'Car_{wildcards.sample}'
        if config['sbatch']==0:
            print("Run on terminal directly")
            os.system(combined_command)
        else:
            print("Run on sbatch")
            os.system(f"python {script_dir}/run_sbatch.py  --sbatch_mode {sbatch_mode}  --job_name {job_name} --cores 1 --mem {requested_memory}G --time {CARLIN_max_run_time} --command '{combined_command}' ") # we use ' in '{command}' to avoid bash expansion
