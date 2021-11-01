import pandas as pd
from pathlib import Path
import help_functions as hf
import os
import sys
sys.path.append('/n/data1/bch/hemonc/camargo/li/snakeTemplate')
#configfile: "config.yaml"  # command line way to set it: --configfile 'path/to/config'
#workdir: config['data_dir'] # set working directory, a command-line way to set it: --directory 'path/to/your/dir'
config['data_dir']=str(os.getcwd())

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
        root_sub_dir=config['data_dir'].split('/DATA/')[1]
        shell("sh {params.script_dir}/add_data_to_dvc.sh {root_sub_dir}")
        
        

