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
        "fastqc_before_pear/multiqc_report.html",
        "fastqc_after_pear/multiqc_report.html",
        expand("fastqc_after_pear/{sample}.trimmed.pear.assembled_fastqc.html",sample=SampleList),
        expand("DARLIN/{sub_dir}/{sample}/DARLIN_analysis.done",sample=SampleList,sub_dir=DARLIN_sub_dir)
    
rule fastqc_before_pear:
    input:
        fq_R1="raw_fastq/{sample}_L001_R1_001.fastq.gz",
        fq_R2="raw_fastq/{sample}_L001_R2_001.fastq.gz"
    output:
        "fastqc_before_pear/{sample}_L001_R1_001_fastqc.html",
        "fastqc_before_pear/{sample}_L001_R2_001_fastqc.html"
    run:
        shell(f"bash {script_dir}/run_fastqc.sh {input.fq_R1} fastqc_before_pear"),
        shell(f"bash {script_dir}/run_fastqc.sh {input.fq_R2} fastqc_before_pear")
        
rule pear:
    input:
        fq_R1="raw_fastq/{sample}_L001_R1_001.fastq.gz",
        fq_R2="raw_fastq/{sample}_L001_R2_001.fastq.gz"
    output:
        "pear_output/{sample}.trimmed.pear.assembled.fastq"
    run:
        shell(f"bash {script_dir}/run_pear.sh {input.fq_R1} {input.fq_R2} pear_output/{wildcards.sample}.trimmed.pear")
        
        
        
rule fastqc_after_pear:
    input:
        "pear_output/{sample}.trimmed.pear.assembled.fastq"
    output:
        "fastqc_after_pear/{sample}.trimmed.pear.assembled_fastqc.html"
    run:
        shell(f"bash {script_dir}/run_fastqc.sh {input} fastqc_after_pear")


rule multiqc_before_pear:
    input:
        expand("fastqc_before_pear/{sample}_L001_R1_001_fastqc.html",sample=SampleList)
    output:
        "fastqc_before_pear/multiqc_report.html"
    run: 
        shell(f"bash {script_dir}/run_multiqc.sh fastqc_before_pear")
        
        
rule multiqc_after_pear:
    input:
        expand("fastqc_after_pear/{sample}.trimmed.pear.assembled_fastqc.html",sample=SampleList)
    output:
        "fastqc_after_pear/multiqc_report.html"
    run: 
        shell(f"bash {script_dir}/run_multiqc.sh fastqc_after_pear")
        
        
rule DARLIN:
    input:
        "pear_output/{sample}.trimmed.pear.assembled.fastq"
    output:
        touch("DARLIN/{sub_dir}/{sample}/DARLIN_analysis.done")
    run:
        input_dir=config['data_dir']+'/pear_output'
        cfg_type=config['cfg_type']
        template=config['template']
        read_cutoff_floor=config['read_cutoff_floor']
        CARLIN_memory_factor=config['CARLIN_memory_factor']
        sbatch=config['sbatch']
        CARLIN_max_run_time=config['CARLIN_max_run_time']
        output_dir=config['data_dir']+f'/DARLIN/{wildcards.sub_dir}'
        read_cutoff_override=int(wildcards.sub_dir.split('_')[-1])
        
        file_size = os.path.getsize(f'{input_dir}/{wildcards.sample}.trimmed.pear.assembled.fastq')/1000000000
        print(f"{wildcards.sample}:   FileSize {file_size} G")
        requested_memory=int(file_size*CARLIN_memory_factor)
        if requested_memory<10:
            requested_memory=10 # at least request 10G memory
        if requested_memory>250:
            requested_memory=250 # do not request more than 250G memory
        print(f"{wildcards.sample}:   Requested memory {requested_memory} G")
        os.makedirs(f'{output_dir}/{wildcards.sample}',exist_ok=True)
        
        shell(f"bash {script_dir}/run_CARLIN.sh {CARLIN_dir} {input_dir} {output_dir} {wildcards.sample} {cfg_type} {template} {read_cutoff_override} {read_cutoff_floor} {requested_memory} {sbatch} {CARLIN_max_run_time}")
