import os
import sys
from pathlib import Path

import pandas as pd

from source import help_functions as hf

#configfile: "config.yaml"  # command line way to set it: --configfile 'path/to/config'
#workdir: config['data_dir'] # set working directory, a command-line way to set it: --directory 'path/to/your/dir'
config['data_dir']=str(os.getcwd())

##################
## preprocessing
################## 

config['CARLIN_dir']=hf.update_CARLIN_dir(config['CARLIN_dir'],config['template'])
print("Updated CARLIN_dir:"+ str(config['CARLIN_dir']))
if 'sbatch_mode' in config.keys():
    sbatch_mode=config['sbatch_mode']
else:
    sbatch_mode='intel-sc3'

print(f'Current work dir: {os.getcwd()}')
if config['template'] == 'Tigre':
    print("------------Warn: remember that the Tigre template is inversed-------------")

if len(config['SampleList'])==0: 
    df=pd.read_csv('raw_fastq/sample_info.csv')
    SampleList=list(df['Sample_long'])
else:
    SampleList=config['SampleList']
print(f'SampleList: {SampleList}')
    
CARLIN_sub_dir=[f"results_cutoff_override_{xx}" for xx in config['read_cutoff_override']]
    

# remove the flag file of the workflow if the sbatch is not actually run to finish
for sample in SampleList:
    if not os.path.exists(f'CARLIN/{CARLIN_sub_dir}/{sample}/CARLIN_analysis_actually.done'):
        if os.path.exists(f'CARLIN/{CARLIN_sub_dir}/{sample}/CARLIN_analysis.done'):
            os.remove(f'CARLIN/{CARLIN_sub_dir}/{sample}/CARLIN_analysis.done') 
        
##################
## start the rules
################## 
rule all:
    input: 
        expand("CARLIN/{sub_dir}/{sample}/CARLIN_analysis.done",sample=SampleList,sub_dir=CARLIN_sub_dir)
 
        
rule CARLIN:
    input:
        "pear_output/{sample}.trimmed.pear.assembled.fastq"
    output:
        touch("CARLIN/{sub_dir}/{sample}/CARLIN_analysis.done")
    run:
        script_dir=config['script_dir']
        CARLIN_dir=config['CARLIN_dir']
        input_dir=config['data_dir']+'/pear_output'
        cfg_type=config['cfg_type']
        template=config['template']
        CARLIN_memory_factor=config['CARLIN_memory_factor']
        sbatch=config['sbatch']
        CARLIN_max_run_time=config['CARLIN_max_run_time']
        output_dir=config['data_dir']+f'/CARLIN/{wildcards.sub_dir}'
        read_cutoff_override=int(wildcards.sub_dir.split('_')[-1])
        
        file_size = os.path.getsize(f'{input_dir}/{wildcards.sample}.trimmed.pear.assembled.fastq')/1000000000
        print(f"{wildcards.sample}:   FileSize {file_size} G")
        requested_memory=int(file_size*CARLIN_memory_factor)
        if requested_memory<10:
            requested_memory=10 # at least request 10G memory
        if requested_memory>250:
            requested_memory=250 # do not request more than 200G memory
        print(f"{wildcards.sample}:   Requested memory {requested_memory} G")
        os.makedirs(f'{output_dir}/{wildcards.sample}',exist_ok=True)
        

        command=f"sh {script_dir}/run_CARLIN.sh {CARLIN_dir} {input_dir} {output_dir} {wildcards.sample} {cfg_type} {template} {read_cutoff_override} {read_cutoff_override} {requested_memory} {0} {CARLIN_max_run_time}"

        
        job_name=f'Car_{wildcards.sample}'
        if config['sbatch']==0:
            print("Run on terminal directly")
            os.system(command)
        else:
            print("Run on sbatch")
            os.system(f"python {script_dir}/run_sbatch.py  --sbatch_mode {sbatch_mode} --job_name {job_name} --cores 1 --mem {requested_memory}G --time {CARLIN_max_run_time} --command '{command}' ") # we use ' in '{command}' to avoid bash expansion
