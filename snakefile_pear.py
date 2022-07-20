import pandas as pd
from pathlib import Path
import os
import sys
from source import help_functions as hf
#configfile: "config.yaml"  # command line way to set it: --configfile 'path/to/config'
#workdir: config['data_dir'] # set working directory, a command-line way to set it: --directory 'path/to/your/dir'
config['data_dir']=str(os.getcwd())

##################
## preprocessing
################## 

config['CARLIN_dir']=hf.update_CARLIN_dir(config['CARLIN_dir'],config['template'])
print("Updated CARLIN_dir:"+ str(config['CARLIN_dir']))


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
        expand("pear_output/{sample}.done",sample=SampleList),
    
        
rule pear:
    input:
        fq_R1="raw_fastq/{sample}_L001_R1_001.fastq.gz",
        fq_R2="raw_fastq/{sample}_L001_R2_001.fastq.gz"
    output:
        touch("pear_output/{sample}.done")
    run:
        script_dir=config['script_dir']
        out_dir=f"pear_output"
        os.makedirs(out_dir,exist_ok=True)
        command=f"sh {script_dir}/run_pear.sh {input.fq_R1} {input.fq_R2} pear_output/{wildcards.sample}.trimmed.pear"
        
        file_size = os.path.getsize(f'raw_fastq/{wildcards.sample}_L001_R1_001.fastq.gz')/1000000000
        print(f"{wildcards.sample}:   FileSize {file_size} G")
        requested_memory=int(file_size*10)
        if requested_memory<1:
            requested_memory=1 # at least request 10G memory
        if requested_memory>200:
            requested_memory=200 # do not request more than 20G memory
        print(f"{wildcards.sample}:   Requested memory {requested_memory} G")
        
        job_name=f'P_{wildcards.sample}'
        if config['sbatch']==0:
            print("Run on terminal directly")
            os.system(command)
        else:
            os.system(f"python {script_dir}/run_sbatch.py --job_name {job_name} --cores 1 --mem {requested_memory}G --time 6 --command '{command}' ") # we use ' in '{command}' to avoid bash expansion