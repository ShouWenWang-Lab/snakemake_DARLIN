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
        "fastqc_after_pear/multiqc_report.html",
        expand("fastqc_after_pear/{sample}.trimmed.pear.assembled_fastqc.html",sample=SampleList),
        expand("CARLIN/{sub_dir}/{sample}/CARLIN_analysis.done",sample=SampleList,sub_dir=CARLIN_sub_dir)
 
        
# rule fastqc_after_pear:
#     input:
#         "pear_output/{sample}.trimmed.pear.assembled.fastq"
#     output:
#         "fastqc_after_pear/{sample}.trimmed.pear.assembled_fastqc.html"
#     params:
#         script_dir=config['script_dir']
#     shell:
#         "sh {params.script_dir}/run_fastqc.sh {input} fastqc_after_pear"

        
        
# rule multiqc_after_pear:
#     input:
#         expand("fastqc_after_pear/{sample}.trimmed.pear.assembled_fastqc.html",sample=SampleList)
#     output:
#         "fastqc_after_pear/multiqc_report.html"
#     params:
#         script_dir=config['script_dir']
#     run: 
#         shell("sh {params.script_dir}/run_multiqc.sh fastqc_after_pear")
        
        
        
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
        read_cutoff_floor=config['read_cutoff_floor']
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
        if requested_memory>200:
            requested_memory=200 # do not request more than 20G memory
        print(f"{wildcards.sample}:   Requested memory {requested_memory} G")
        os.makedirs(f'{output_dir}/{wildcards.sample}',exist_ok=True)
        
        shell(f"sh {script_dir}/run_CARLIN.sh {CARLIN_dir} {input_dir} {output_dir} {wildcards.sample} {cfg_type} {template} {read_cutoff_override} {read_cutoff_floor} {requested_memory} {sbatch} {CARLIN_max_run_time}")

