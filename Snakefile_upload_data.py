import pandas as pd
from pathlib import Path
from source import help_functions as hf
import os
import sys
#configfile: "config.yaml"  # command line way to set it: --configfile 'path/to/config'
#workdir: config['data_dir'] # set working directory, a command-line way to set it: --directory 'path/to/your/dir'
config['data_dir']=str(os.getcwd())

print(f"Current working dir: {os.getcwd()}")

if len(config['SampleList'])==0: 
    df=pd.read_csv('raw_fastq/sample_info.csv')
    SampleList=list(df['Sample_long'])
else:
    SampleList=config['SampleList']
    
SampleList.append('merge_all')
print(f'SampleList: {SampleList}')
    
CARLIN_sub_dir=[f"results_cutoff_override_{xx}" for xx in config['read_cutoff_override']]
##################
## start the rules
################## 


rule all:
    input: 
        "transfer_data.done",
        "dvc_upload.done"

        
        
rule transfer_data:
    output:
        touch("transfer_data.done")
    params:
        script_dir=config['script_dir']
    run:
        root_sub_dir=config['data_dir'].split('/DATA/')[1]
        shell("sh {params.script_dir}/transfer_data.sh {root_sub_dir}")
        
        
rule dvc:
    output:
        touch("dvc_upload.done")
    params:
        script_dir=config['script_dir']
    run:
        #os.system(f'dvc remove CARLIN.dvc') # to be removed later
        data_dir=config['data_dir']
        # add the CARLIN folder
        for cur_dir in CARLIN_sub_dir:
            for sample in SampleList:
                os.system(f"dvc add {data_dir}/CARLIN/{cur_dir}/{sample}")
        
        os.system(f"sh {params.script_dir}/add_data_to_dvc.sh {data_dir}") # skip adding the CARLIN folder, includes "dvc push"
        
        
        

