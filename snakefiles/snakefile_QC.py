import pandas as pd
from pathlib import Path
import os
import sys
from darlin import help_functions as hf
from darlin.settings import script_dir, QC_dir, ref_dir, CARLIN_dir
#configfile: "config.yaml"  # command line way to set it: --configfile 'path/to/config'
#workdir: config['data_dir'] # set working directory, a command-line way to set it: --directory 'path/to/your/dir'
config['data_dir']=str(os.getcwd())

##################
## preprocessing
################## 

print(f'Current work dir: {os.getcwd()}')

if len(config['SampleList'])==0: 
    df=pd.read_csv('raw_fastq/sample_info.csv')
    SampleList=list(df['Sample_long'])
else:
    SampleList=config['SampleList']
    
# this is to make the pipeline compatible with earlier bulk config files
if cfg_type.startswith('Bulk') and ('read_cutoff_UMI_override' not in config.keys()) and ('read_cutoff_override' in config.keys()):
    config['read_cutoff_UMI_override']=config['read_cutoff_override']
    config['read_cutoff_CB_override']=10 
    
CARLIN_sub_dir=[f"results_cutoff_override_{xx}" for xx in config['read_cutoff_UMI_override']]
    
        
##################
## start the rules
################## 
rule all:
    input: 
        expand("fastqc_before_pear/{sample}.done",sample=SampleList),
        expand("fastqc_after_pear/{sample}.done",sample=SampleList,sub_dir=CARLIN_sub_dir)
    
rule fastqc_before_pear:
    input:
        fq_R1="raw_fastq/{sample}_L001_R1_001.fastq.gz",
        fq_R2="raw_fastq/{sample}_L001_R2_001.fastq.gz"
    output:
        touch("fastqc_before_pear/{sample}.done")
    run:
        commands=[f"sh {script_dir}/run_fastqc.sh {input.fq_R1} fastqc_before_pear", f"sh {script_dir}/run_fastqc.sh {input.fq_R2} fastqc_before_pear"]

        file_size = os.path.getsize(f'{input.fq_R1}')/1000000000
        print(f"{wildcards.sample}:   FileSize {file_size} G")
        requested_memory=int(file_size*10)
        if requested_memory<1:
            requested_memory=1 # at least request 10G memory
        if requested_memory>200:
            requested_memory=200 # do not request more than 20G memory
        print(f"{wildcards.sample}:   Requested memory {requested_memory} G")
        
        job_name=f'QC_{wildcards.sample}'
        if config['sbatch']==0:
            print("Run on terminal directly")
            for command in commands:
                os.system(command)
        else:
            for command in commands:
                os.system(f"python {script_dir}/run_sbatch.py --job_name {job_name} --cores 1 --mem {requested_memory}G --time 1 --command '{command}' ") # we use ' in '{command}' to avoid bash expansion
        
        
rule fastqc_after_pear:
    input:
        "pear_output/{sample}.trimmed.pear.assembled.fastq"
    output:
        touch("fastqc_after_pear/{sample}.done")
    run:
        command=f"sh {script_dir}/run_fastqc.sh {input} fastqc_after_pear"
        
        file_size = os.path.getsize(f'{input}')/1000000000
        print(f"{wildcards.sample}:   FileSize {file_size} G")
        requested_memory=int(file_size*10)
        if requested_memory<1:
            requested_memory=1 # at least request 10G memory
        if requested_memory>200:
            requested_memory=200 # do not request more than 20G memory
        print(f"{wildcards.sample}:   Requested memory {requested_memory} G")
        
        job_name=f'QP_{wildcards.sample}'
        if config['sbatch']==0:
            print("Run on terminal directly")
            os.system(command)
        else:
            os.system(f"python {script_dir}/run_sbatch.py --job_name {job_name} --cores 1 --mem {requested_memory}G --time 1 --command '{command}' ") # we use ' in '{command}' to avoid bash expansion
            