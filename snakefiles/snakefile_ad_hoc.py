import os
import sys
from pathlib import Path

import pandas as pd

from darlin import help_functions as hf
from darlin import settings
from darlin.settings import script_dir, QC_dir, ref_dir, CARLIN_dir
#configfile: "config.yaml"  # command line way to set it: --configfile 'path/to/config'
#workdir: config['data_dir'] # set working directory, a command-line way to set it: --directory 'path/to/your/dir'
config['data_dir']=str(os.getcwd())
global_file='/Users/shouwenwang/Dropbox (HMS)/shared_folder_with_Li/DATA/scripts/running_script.m'
if not os.path.exists(global_file):
    os.system(f"touch '{global_file}'")

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
    
CARLIN_sub_dir=[f"results_cutoff_override_{xx}" for xx in config['read_cutoff_override']]
print(f"Subdir: {CARLIN_sub_dir}")
#CARLIN_sub_dir="results_cutoff_override_"+str(config['read_cutoff_override'])

##################
## start the rules
################## 


rule all:
    input: 
        expand("CARLIN/{sub_dir}/{sample}/ad_hoc.done",sub_dir=CARLIN_sub_dir,sample=SampleList),
        
rule analysis:
    output:
        touch("CARLIN/{sub_dir}/{sample}/ad_hoc.done")
    params:
        script_dir=script_dir,
        CARLIN_dir=CARLIN_dir,
        Samples=','.join(SampleList),
        template=config['template']
    run:
        input_dir=config['data_dir']+f'/CARLIN/{wildcards.sub_dir}'
        
        # command=f"$matlab_terminal  -nodisplay -nosplash -nodesktop -r \"cd \'{params.CARLIN_dir}\'; output_selected_from_summary(\'{wildcards.sample}\',\'{input_dir}\',\'{params.template}\'); exit;\""
        command=f"cd '{params.CARLIN_dir}'; output_selected_from_summary('{wildcards.sample}','{input_dir}','{params.template}');"
        with open(global_file,'a') as f:
            f.write(command)
            f.write('\n\n')

