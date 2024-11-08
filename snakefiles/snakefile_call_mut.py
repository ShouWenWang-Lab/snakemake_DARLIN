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

# parameters
python_sub_dir=config['template'] + '_' + config['python_DARLIN_pipeline']['output_folder']
kernel=config['python_DARLIN_pipeline']['kernel']
template=config['template']

# remove the flag file of the workflow if the sbatch is not actually run to finish
for sample in SampleList:
    if not os.path.exists('DARLIN/'+config['template']+f'_cutoff_override_1/{sample}/DARLIN_analysis_actually.done'):
        if os.path.exists('DARLIN/'+config['template']+f'results_cutoff_override_1/{sample}/DARLIN_analysis.done'):
            os.remove('DARLIN/'+config['template']+f'results_cutoff_override_1/{sample}/DARLIN_analysis.done') 

os.path.exists('slim_fastq') or os.makedirs('slim_fastq')

##################
## start the rules
################## 
rule all:
    input: 
        expand(f"DARLIN/{python_sub_dir}/{{sample}}/final.done",sample=SampleList)

rule gen_slim_fastq:
    input:
        called_bc_file=f"DARLIN/{python_sub_dir}/{{sample}}/called_barcodes_by_SW_method.csv"
    output:
        fq_R1="slim_fastq_"+config['template']+"/{sample}_R1.fastq.gz",
        fq_R2="slim_fastq_"+config['template']+"/{sample}_R2.fastq.gz"
    run:
        job_name=f'Car_{wildcards.sample}'
        if cfg_type=='sc10xV3':
            command=f"papermill {QC_dir}/extra-gen_10xv3_fastq.ipynb -k {kernel} DARLIN/{python_sub_dir}/{wildcards.sample}/extra-gen_10xv3_fastq.ipynb -p called_bc_file {input.called_bc_file} -p locus {template} -p R1_file {output.fq_R1} -p R2_file {output.fq_R2}"
        elif cfg_type=='scCamellia':
            command=f"papermill {QC_dir}/extra-gen_scCamellia_fastq.ipynb -k {kernel} DARLIN/{python_sub_dir}/{wildcards.sample}/extra-gen_scCamellia_fastq.ipynb -p called_bc_file {input.called_bc_file} -p locus {template} -p R1_file {output.fq_R1} -p R2_file {output.fq_R2}"
        if config['sbatch']==0:
            print("Run on terminal directly")
            os.system(command)
        else:
            print("Run on sbatch")
            os.system(f"python {script_dir}/run_sbatch.py  --sbatch_mode {sbatch_mode}  --job_name {job_name} --cores 1 --mem {requested_memory}G --time {CARLIN_max_run_time} --command '{command}' ") # we use ' in '{command}' to avoid bash expansion

rule DARLIN:        
    input:
        # fq_R1="slim_fastq/{sample}_L001_R1_001.fastq.gz",
        # fq_R2="slim_fastq/{sample}_L001_R2_001.fastq.gz"
        fq_R1="slim_fastq_"+config['template']+"/{sample}_R1.fastq.gz",
        fq_R2="slim_fastq_"+config['template']+"/{sample}_R2.fastq.gz"
    output:
        "DARLIN/"+config['template']+"_cutoff_override_1/{sample}/Actaul_CARLIN_seq.txt",
        "DARLIN/"+config['template']+"_cutoff_override_1/{sample}/AlleleAnnotations.txt"
    run:
        output_dir=config['data_dir']+'/DARLIN/'+config['template']+'_cutoff_override_1'
        input_dir=config['data_dir']+'/slim_fastq_'+config['template']
        print(input_dir)

        CARLIN_memory_factor=config['CARLIN_memory_factor']
        sbatch=config['sbatch']
        CARLIN_max_run_time=config['CARLIN_max_run_time']
        read_cutoff_UMI_override=1
        read_cutoff_CB_override=read_cutoff_UMI_override
        
        file_size = os.path.getsize(f'{input.fq_R1}')*5/1000000000
        print(f"{wildcards.sample}:   FileSize {file_size} G")
        requested_memory=5
        if requested_memory<20:
            requested_memory=20 # at least request 20G memory
        if requested_memory>250:
            requested_memory=250 # do not request more than 200G memory
        print(f"{wildcards.sample}:   Requested memory {requested_memory} G")
        os.makedirs(f'{output_dir}/{wildcards.sample}',exist_ok=True)
        
        # do not run sbatch within this command
        command_1=f"bash {script_dir}/run_CARLIN.sh {CARLIN_dir} {input_dir} {output_dir} {wildcards.sample} {cfg_type} {template} {read_cutoff_UMI_override} {read_cutoff_CB_override} {requested_memory} {0} {CARLIN_max_run_time}"

        combined_command=f"{command_1}"
        
        
        job_name=f'Car_{wildcards.sample}'
        if config['sbatch']==0:
            print("Run on terminal directly")
            os.system(combined_command)
        else:
            print("Run on sbatch")
            os.system(f"python {script_dir}/run_sbatch.py  --sbatch_mode {sbatch_mode}  --job_name {job_name} --cores 1 --mem {requested_memory}G --time {CARLIN_max_run_time} --command '{combined_command}' ") # we use ' in '{command}' to avoid bash expansion

rule join_results:
    input:
        "DARLIN/"+config['template']+"_cutoff_override_1/{sample}/Actaul_CARLIN_seq.txt",
        "DARLIN/"+config['template']+"_cutoff_override_1/{sample}/AlleleAnnotations.txt"
    output:
        touch(f"DARLIN/{python_sub_dir}/{{sample}}/final.done")
    run:
        called_bc_file = f"DARLIN/{python_sub_dir}/{wildcards.sample}/called_barcodes_by_SW_method.csv"
        carlin_seq_file = "DARLIN/"+config['template']+f"_cutoff_override_1/{wildcards.sample}/Actaul_CARLIN_seq.txt"
        allele_annot_file = "DARLIN/"+config['template']+f"_cutoff_override_1/{wildcards.sample}/AlleleAnnotations.txt"
        final_file = f"DARLIN/{python_sub_dir}/{wildcards.sample}/called_barcodes_by_SW_method_allele.csv"
        command_1=f"python {script_dir}/extra-join_results.py {called_bc_file} {carlin_seq_file} {allele_annot_file} {final_file}"
        combined_command=f"{command_1}"
        job_name=f'Car_{wildcards.sample}'
        if config['sbatch']==0:
            print("Run on terminal directly")
            os.system(combined_command)
        else:
            print("Run on sbatch")
            os.system(f"python {script_dir}/run_sbatch.py  --sbatch_mode {sbatch_mode}  --job_name {job_name} --cores 1 --mem {requested_memory}G --time {CARLIN_max_run_time} --command '{combined_command}' ") # we use ' in '{command}' to avoid bash expansion
