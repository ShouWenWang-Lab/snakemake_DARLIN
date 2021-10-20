import pandas as pd
from pathlib import Path
#configfile: "config.yaml"

df=pd.read_csv('raw_fastq/sample_info.csv')
SampleList=list(df['Sample_long'])
#CARLIN_analysis=''
#SampleList='LL605-LK_S3'


rule all:
    input: 
        "fastqc_before_pear/multiqc_report.html",
        "fastqc_after_pear/multiqc_report.html",
        expand("fastqc_after_pear/{sample}.trimmed.pear.assembled_fastqc.html",sample=SampleList),
        #expand("CARLIN/{sample}/Summary.mat",sample=SampleList)
        expand("CARLIN/{sample}/CARLIN_analysis.done",sample=SampleList)
        
# rule download_fastq:
#     params:
#         project_name=config["project_name"],
#         project_ID=config["project_ID"],
#         root_path=config['data_dir']
#     output: 
#         touch("raw_fastq/download_fastq.done")
#     run:
#         #shell("./download_fastq.sh {params.project_name} {params.project_ID} {params.root_path}")
#         import glob
#         import re
#         SampleList = []
#         locations=[]
#         for file in glob.glob(config['data_dir'] + "/raw_fastq", recursive=True):
#             if '_L001_R1_001.fastq.gz' in file:
#                 sample_name = file.split('/')[-1].split('_L001_R1')[0]
#                 SampleList.append(sample_name)

#         df=pd.DataFrame({'Sample':SampleList})     
#         df.to_csv('raw_fastq/sample_info.csv')
    
rule fastqc_before_pear:
    input:
        fq_R1="raw_fastq/{sample}_L001_R1_001.fastq.gz",
        fq_R2="raw_fastq/{sample}_L001_R2_001.fastq.gz"
        #flag_file="raw_fastq/download_fastq.done"
    output:
        "fastqc_before_pear/{sample}_L001_R1_001_fastqc.html",
        "fastqc_before_pear/{sample}_L001_R2_001_fastqc.html"
    params:
        script_dir=config['script_dir']
    run:
        shell("sh {params.script_dir}/run_fastqc.sh {input.fq_R1} fastqc_before_pear"),
        shell("sh {params.script_dir}/run_fastqc.sh {input.fq_R2} fastqc_before_pear")
        
rule pear:
    input:
        fq_R1="raw_fastq/{sample}_L001_R1_001.fastq.gz",
        fq_R2="raw_fastq/{sample}_L001_R2_001.fastq.gz"
        #flag_file="raw_fastq/download_fastq.done"
    output:
        "pear_output/{sample}.trimmed.pear.assembled.fastq"
    params:
        script_dir=config['script_dir'],
        is_RNA=config['is_RNA']
    shell:
        "sh {params.script_dir}/run_pear.sh {input.fq_R1} {input.fq_R2} pear_output/{wildcards.sample}.trimmed.pear {params.is_RNA}"
        
        
        
rule fastqc_after_pear:
    input:
        "pear_output/{sample}.trimmed.pear.assembled.fastq"
    output:
        "fastqc_after_pear/{sample}.trimmed.pear.assembled_fastqc.html"
    params:
        script_dir=config['script_dir']
    shell:
        "sh {params.script_dir}/run_fastqc.sh {input} fastqc_after_pear"


rule multiqc_before_pear:
    input:
        expand("fastqc_before_pear/{sample}_L001_R1_001_fastqc.html",sample=SampleList)
    output:
        "fastqc_before_pear/multiqc_report.html"
    params:
        script_dir=config['script_dir']
    shell:
        "sh {params.script_dir}/run_multiqc.sh fastqc_before_pear"
        
        
rule multiqc_after_pear:
    input:
        expand("fastqc_after_pear/{sample}.trimmed.pear.assembled_fastqc.html",sample=SampleList)
    output:
        "fastqc_after_pear/multiqc_report.html"
    params:
        script_dir=config['script_dir']
    shell:
        "sh {params.script_dir}/run_multiqc.sh fastqc_after_pear"
        
        
rule CARLIN:
    input:
        "pear_output/{sample}.trimmed.pear.assembled.fastq"
    output:
        #"CARLIN/{sample}/Summary.mat" # if not succeeded, this file will be deleted
        touch("CARLIN/{sample}/CARLIN_analysis.done")
    params:
        script_dir=config['script_dir'],
        CARLIN_dir=config['CARLIN_dir'],
        input_dir=config['data_dir']+'/pear_output',
        output_dir=config['data_dir']+'/CARLIN',
        cfg_type=config['cfg_type'],
        template=config['template'],
        read_cutoff_override=config['read_cutoff_override'],
        read_cutoff_floor=config['read_cutoff_floor'],
        CARLIN_memory_factor=config['CARLIN_memory_factor'],
        sbatch=config['sbatch'],
        CARLIN_max_run_time=config['CARLIN_max_run_time']
    run:
        import os 
        file_size = os.path.getsize(f'{params.input_dir}/{wildcards.sample}.trimmed.pear.assembled.fastq')/1000000
        print(f"{wildcards.sample}:   FileSize {file_size} M")
        requested_memory=int(file_size*params.CARLIN_memory_factor)
        print(f"{wildcards.sample}:   Requested memory {requested_memory} M")
        shell("sh {params.script_dir}/run_CARLIN.sh {params.CARLIN_dir} {params.input_dir} {params.output_dir} {wildcards.sample} {params.cfg_type} {params.template} {params.read_cutoff_override} {params.read_cutoff_floor} {requested_memory} {params.sbatch} {params.CARLIN_max_run_time}")