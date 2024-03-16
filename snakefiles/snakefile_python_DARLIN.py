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

CARLIN_dir=hf.update_CARLIN_dir(CARLIN_dir,config['template'])
print("Updated CARLIN_dir:"+ str(CARLIN_dir))
cfg_type=config['cfg_type']

if len(config['SampleList'])==0: 
    df=pd.read_csv('raw_fastq/sample_info.csv')
    SampleList=list(df['Sample_long'])
else:
    SampleList=config['SampleList']

    
# this is to make the pipeline compatible with earlier bulk config files
if cfg_type.startswith('Bulk') and ('read_cutoff_UMI_override' not in config.keys()) and ('read_cutoff_override' in config.keys()):
    config['read_cutoff_UMI_override']=config['read_cutoff_override']
    config['read_cutoff_CB_override']=10 
    
    
# parameters
coarse_grained_readcutoff_floor=config['python_DARLIN_pipeline']['coarse_grained_readcutoff_floor']
distance_relative_threshold=config['python_DARLIN_pipeline']['distance_relative_threshold']
read_ratio_threshold=config['python_DARLIN_pipeline']['read_ratio_threshold']
seq_3prime_upper_N=config['python_DARLIN_pipeline']['seq_3prime_upper_N']
output_folder=config['python_DARLIN_pipeline']['output_folder']
kernel=config['python_DARLIN_pipeline']['kernel']
DARLIN_sub_dir=[output_folder]

    
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
        fq_R1="raw_fastq/{sample}_L001_R1_001.fastq.gz",
        fq_R2="raw_fastq/{sample}_L001_R2_001.fastq.gz"
    output:
        touch("DARLIN/{sub_dir}/{sample}/DARLIN_analysis.done")
    run:
        output_dir=config['data_dir']+f'/DARLIN/{wildcards.sub_dir}'
        
        if (cfg_type not in ['scCamellia','sc10xV3']):
            raise ValueError("This pipeline is only intended for scCamellia or sc10xV3")
            
        else:
            input_dir=config['data_dir']+'/raw_fastq'
            
        print(input_dir)
        template=config['template']
        CARLIN_memory_factor=config['CARLIN_memory_factor']
        sbatch=config['sbatch']
        CARLIN_max_run_time=config['CARLIN_max_run_time']
        file_size = os.path.getsize(f'{input.fq_R1}')*5/1000000000
        print(f"{wildcards.sample}:   FileSize {file_size} G")
        requested_memory=int(file_size*CARLIN_memory_factor)
        if requested_memory<20:
            requested_memory=20 # at least request 20G memory
        if requested_memory>250:
            requested_memory=250 # do not request more than 200G memory
        print(f"{wildcards.sample}:   Requested memory {requested_memory} G")
        os.makedirs(f'{output_dir}/{wildcards.sample}',exist_ok=True)
        
        print("----generate report -----")
        data_dir=config['data_dir']
        
        if cfg_type=='scCamellia':
            command=f"""
            papermill  {QC_dir}/single_cell_DARLIN-scCamellia.ipynb -k {kernel} {output_dir}/{wildcards.sample}/single_cell_DARLIN-scCamellia.ipynb  -p sample {wildcards.sample} -p template {template} -p data_path {data_dir} -p ref_dir {ref_dir} -p output_dir {output_dir}/{wildcards.sample} -p cfg {cfg_type} -p coarse_grained_readcutoff_floor {coarse_grained_readcutoff_floor} -p distance_relative_threshold {distance_relative_threshold} -p read_ratio_threshold {read_ratio_threshold} -p seq_3prime_upper_N {seq_3prime_upper_N}
            jupyter nbconvert --to html {output_dir}/{wildcards.sample}/single_cell_DARLIN-scCamellia.ipynb
            """
        elif cfg_type=='sc10xV3':
            command=f"""
            papermill  {QC_dir}/single_cell_DARLIN-10x.ipynb -k {kernel} {output_dir}/{wildcards.sample}/single_cell_DARLIN-10x.ipynb  -p sample {wildcards.sample} -p template {template} -p data_path {data_dir} -p output_dir {output_dir}/{wildcards.sample} -p cfg {cfg_type} -p ref_dir {ref_dir} -p coarse_grained_readcutoff_floor {coarse_grained_readcutoff_floor} -p distance_relative_threshold {distance_relative_threshold} -p read_ratio_threshold {read_ratio_threshold} -p seq_3prime_upper_N {seq_3prime_upper_N}
            jupyter nbconvert --to html {output_dir}/{wildcards.sample}/single_cell_DARLIN-10x.ipynb
            """
            
        
        job_name=f'scL_{wildcards.sample}'
        if config['sbatch']==0:
            print("Run on terminal directly")
            os.system(command)
        else:
            print("Run on sbatch")
            os.system(f"python {script_dir}/run_sbatch.py --job_name {job_name} --cores 1 --mem {requested_memory}G --time 1 --command '{command}' ") # we use ' in '{command}' to avoid bash expansion
