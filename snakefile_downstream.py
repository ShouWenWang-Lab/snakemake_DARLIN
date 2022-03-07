import pandas as pd
from pathlib import Path
from source import help_functions as hf
import os
import sys
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
    print("Warn: remeber that the Tigre template is inversed")

if len(config['SampleList'])==0: 
    df=pd.read_csv('raw_fastq/sample_info.csv')
    SampleList=list(df['Sample_long'])
else:
    SampleList=config['SampleList']
print(f'SampleList: {SampleList}')
    
CARLIN_sub_dir=[f"results_cutoff_override_{xx}" for xx in config['read_cutoff_override']]
print(f"Subdir: {CARLIN_sub_dir}")
#CARLIN_sub_dir="results_cutoff_override_"+str(config['read_cutoff_override'])


##################
## start the rules
################## 


rule all:
    input: 
        expand("CARLIN/{sub_dir}/merge_all/refined_results.csv",sub_dir=CARLIN_sub_dir),
        expand("CARLIN/{sub_dir}/merge_all/all_insertion_pattern.png",sub_dir=CARLIN_sub_dir),
        expand("CARLIN/{sub_dir}/merge_all/allele_annotation.mat",sub_dir=CARLIN_sub_dir),
        expand("CARLIN/{sub_dir}/merge_all/CARLIN_report.html",sub_dir=CARLIN_sub_dir),
        
rule merge_all_sample:
    output:
        "CARLIN/{sub_dir}/merge_all/Bank.mat",
        "CARLIN/{sub_dir}/merge_all/allele_breakdown_by_sample.mat",
        "CARLIN/{sub_dir}/merge_all/allele_annotation.mat"
    params:
        script_dir=config['script_dir'],
        CARLIN_dir=config['CARLIN_dir'],
        Samples=','.join(SampleList),
        template=config['template']
    run:
        input_dir=config['data_dir']+f'/CARLIN/{wildcards.sub_dir}'
        shell("sh {params.script_dir}/merge_sample.sh  {params.CARLIN_dir} {input_dir} {params.Samples}  {params.template}")
        
        
        
rule CARLIN_csv:
    input:
        "CARLIN/{sub_dir}/merge_all/Bank.mat"
    output:
        "CARLIN/{sub_dir}/merge_all/refined_results.csv"
    params:
        script_dir=config['script_dir'],
        CARLIN_dir=config['CARLIN_dir'],
        template=config['template'],
        Samples=','.join(SampleList+["merge_all"])
    run:
        input_dir=config['data_dir']+f'/CARLIN/{wildcards.sub_dir}'
        print(f"SampleList:   {params.Samples}")
        shell("sh {params.script_dir}/generate_csv.sh {params.CARLIN_dir} {input_dir} {params.Samples}  {params.template}")
        
        

rule more_plots:
    input:
        "CARLIN/{sub_dir}/merge_all/refined_results.csv"
    output:
        "CARLIN/{sub_dir}/merge_all/all_insertion_pattern.png"
    params:
        script_dir=config['script_dir'],
        Samples=','.join(SampleList)
    run:
        input_dir=config['data_dir']+f'/CARLIN/{wildcards.sub_dir}'
        output_dir=config['data_dir']+f'/CARLIN/{wildcards.sub_dir}'
        shell("python {params.script_dir}/make_more_plots.py --input_dir {input_dir} --SampleList {params.Samples} --output_dir {output_dir}")
        #shell("python {params.script_dir}/clonal_analysis.py --data_path {input_dir} --SampleList {params.Samples}") #The computationally poor
        
        for sample in SampleList:
            new_input_dir=f'{output_dir}/{sample}'
            shell(f"python {params.script_dir}/plot_cumulative_insert_del_freq.py --input_dir {new_input_dir}")
        
        
# make sure you have a kernel that could run this notebook, and set it as the default kernel for this notebook
rule generate_report:
    input:
        "CARLIN/{sub_dir}/merge_all/all_insertion_pattern.png"
    output:
        "CARLIN/{sub_dir}/merge_all/CARLIN_report.html"
    run:
        script_dir=config['script_dir']
        Samples=','.join(SampleList)
        data_dir=os.path.join(config['data_dir'],'CARLIN', wildcards.sub_dir)
        shell(f"papermill  {script_dir}/CARLIN_report.ipynb  {data_dir}/merge_all/CARLIN_report.ipynb  -p data_dir {data_dir} -p Samples {Samples}")
        shell("jupyter nbconvert --to html {data_dir}/merge_all/CARLIN_report.ipynb")
        
        
        