import os
import sys
from pathlib import Path

import pandas as pd
from scipy.io import loadmat
from tqdm import tqdm

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

print(f'Current work dir: {os.getcwd()}')
if config['template'] == 'Tigre':
    print("Warn: remeber that the Tigre template is inversed")

if len(config['SampleList'])==0: 
    df=pd.read_csv('raw_fastq/sample_info.csv')
    SampleList=list(df['Sample_long'])
else:
    SampleList=config['SampleList']
print(f'SampleList: {SampleList}')
    
# this is to make the pipeline compatible with earlier bulk config files
if cfg_type.startswith('Bulk') and ('read_cutoff_UMI_override' not in config.keys()) and ('read_cutoff_override' in config.keys()):
    config['read_cutoff_UMI_override']=config['read_cutoff_override']
    config['read_cutoff_CB_override']=10 
    
DARLIN_sub_dir=[f"results_cutoff_override_{xx}" for xx in config['read_cutoff_UMI_override']]
print(f"Subdir: {DARLIN_sub_dir}")
#DARLIN_sub_dir="results_cutoff_override_"+str(config['read_cutoff_override'])


##################
## start the rules
################## 


rule all:
    input: 
        expand("DARLIN/{sub_dir}/merge_all/refined_results.csv",sub_dir=DARLIN_sub_dir),
        expand("DARLIN/{sub_dir}/merge_all/DARLIN_report.html",sub_dir=DARLIN_sub_dir),
        
        

rule plots:
    output:
        "DARLIN/{sub_dir}/merge_all/refined_results.csv"
    params:
        script_dir=script_dir,
        CARLIN_dir=CARLIN_dir,
        template=config['template'],
    run:
        hf.set_rcParams()
        input_dir=config['data_dir']+f'/DARLIN/{wildcards.sub_dir}'
        # print("---- Allele analysis -----")
        # hf.analyze_allele_frequency_count(input_dir,SampleList)
        print("---- Sample statistics csv (also do allele analysis) -----")
        hf.generate_csv(input_dir,SampleList,cfg_type=cfg_type)

        file_path=f'DARLIN/{wildcards.sub_dir}/merge_all/refined_results.csv'
        df_results=pd.read_csv(file_path,index_col=0).sort_values('sample')
        df_results.to_csv(file_path)

        print("---- Plot sample statistics -----")
        hf.plot_data_statistics_across_samples(input_dir)
        
        print("---- Cumulative insertion/deletion plot per sample -----")

        
        for sample in SampleList:
            data = loadmat(f"{input_dir}/{sample}/allele_annotation.mat")
            freq = data["allele_freqs"].flatten()[1:]
            allele_annot = data["AlleleAnnotation"].flatten()[1:]
            allele_annot=[x[0] for x in allele_annot]
            df_input=pd.DataFrame({'allele':allele_annot,'UMI_count':freq})
            df_input.to_csv(f"{input_dir}/{sample}/allele_UMI_count.csv",index=0)
            hf.plot_cumulative_insert_del_freq(df_input, input_dir+'/'+sample)
        df_input=pd.read_csv(f"{input_dir}/merge_all/allele_UMI_count.csv")
        hf.plot_cumulative_insert_del_freq(df_input, input_dir+'/merge_all')

        print("---- Insertion pattern -----")
        hf.plot_insertion_patterns(input_dir,SampleList)
        
# make sure you have a kernel that could run this notebook, and set it as the default kernel for this notebook
rule generate_report:
    input:
        "DARLIN/{sub_dir}/merge_all/refined_results.csv"
    output:
        "DARLIN/{sub_dir}/merge_all/DARLIN_report.html"
    run:
        print("----generate report -----")
        Samples=','.join(SampleList)
        data_dir=os.path.join(config['data_dir'],'DARLIN', wildcards.sub_dir)
        shell(f"papermill  {QC_dir}/DARLIN_report.ipynb  {data_dir}/merge_all/DARLIN_report.ipynb  -p data_dir {data_dir} -p Samples {Samples}")
        shell(f"jupyter nbconvert --to html {data_dir}/merge_all/DARLIN_report.ipynb")
        
        
        