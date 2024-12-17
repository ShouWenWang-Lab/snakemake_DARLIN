import pandas as pd
from pathlib import Path
from darlin import help_functions as hf
import os
import sys
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
if config['template'] == 'Tigre':
    print("Warn: remeber that the Tigre template is inversed")

if len(config['SampleList'])==0: 
    df=pd.read_csv('raw_fastq/sample_info.csv')
    SampleList=list(df['Sample_long'])
else:
    SampleList=config['SampleList']
print(f'SampleList: {SampleList}')
    
DARLIN_sub_dir=[f"results_cutoff_override_{xx}" for xx in config['read_cutoff_override']]
print(f"Subdir: {DARLIN_sub_dir}")
#DARLIN_sub_dir="results_cutoff_override_"+str(config['read_cutoff_override'])


##################
## start the rules
################## 


rule all:
    input: 
        expand("DARLIN/{sub_dir}/merge_all/refined_results.csv",sub_dir=DARLIN_sub_dir),
        expand("DARLIN/{sub_dir}/merge_all/all_insertion_pattern.png",sub_dir=DARLIN_sub_dir),
        expand("DARLIN/{sub_dir}/merge_all/allele_annotation.mat",sub_dir=DARLIN_sub_dir),
        expand("DARLIN/{sub_dir}/merge_all/DARLIN_report.html",sub_dir=DARLIN_sub_dir),
        
        
rule merge_all_sample:
    output:
        "DARLIN/{sub_dir}/merge_all/allele_annotation.mat"
    params:
        script_dir=script_dir,
        CARLIN_dir=CARLIN_dir,
        Samples=','.join(SampleList),
        template=config['template']
    run:
        print("----Merge all samples -----")
        input_dir=config['data_dir']+f'/DARLIN/{wildcards.sub_dir}'
        shell("bash {params.script_dir}/merge_sample.sh  {params.CARLIN_dir} {input_dir} {params.Samples}  {params.template}")
        
        
        
rule DARLIN_csv:
    input:
        "DARLIN/{sub_dir}/merge_all/allele_annotation.mat"
    output:
        "DARLIN/{sub_dir}/merge_all/refined_results.csv"
    params:
        script_dir=script_dir,
        CARLIN_dir=CARLIN_dir,
        template=config['template'],
        Samples=','.join(SampleList+["merge_all"])
    run:
        print("----generate csv results -----")
        input_dir=config['data_dir']+f'/DARLIN/{wildcards.sub_dir}'
        print(f"SampleList:   {params.Samples}")
        shell("bash {params.script_dir}/generate_csv.sh {params.CARLIN_dir} {input_dir} {params.Samples}  {params.template}")
        
        

rule more_plots:
    input:
        "DARLIN/{sub_dir}/merge_all/refined_results.csv"
    output:
        "DARLIN/{sub_dir}/merge_all/all_insertion_pattern.png"
    params:
        script_dir=script_dir,
        Samples=','.join(SampleList)
    run:
        print("----more plots -----")
        input_dir=config['data_dir']+f'/DARLIN/{wildcards.sub_dir}'
        output_dir=config['data_dir']+f'/DARLIN/{wildcards.sub_dir}'
        shell("python {params.script_dir}/make_more_plots.py --input_dir {input_dir} --SampleList {params.Samples} --output_dir {output_dir}")
        #shell("python {params.script_dir}/clonal_analysis.py --data_path {input_dir} --SampleList {params.Samples}") #The computationally poor
        
        for sample in SampleList:
            new_input_dir=f'{output_dir}/{sample}'
            shell(f"python {params.script_dir}/plot_cumulative_insert_del_freq.py --input_dir {new_input_dir}")
        
        
# make sure you have a kernel that could run this notebook, and set it as the default kernel for this notebook
rule generate_report:
    input:
        "DARLIN/{sub_dir}/merge_all/all_insertion_pattern.png"
    output:
        "DARLIN/{sub_dir}/merge_all/DARLIN_report.html"
    run:
        print("----generate report -----")
        Samples=','.join(SampleList)
        data_dir=os.path.join(config['data_dir'],'CARLIN', wildcards.sub_dir)
        shell(f"papermill  {QC_dir}/DARLIN_report.ipynb  {data_dir}/merge_all/DARLIN_report.ipynb  -p data_dir {data_dir} -p Samples {Samples}")
        shell("jupyter nbconvert --to html {data_dir}/merge_all/DARLIN_report.ipynb")
        
        
        