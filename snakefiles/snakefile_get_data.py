import pandas as pd
from pathlib import Path
import os
from darlin.settings import script_dir, QC_dir, ref_dir, CARLIN_dir
#configfile: "config.yaml"  # command line way to set it: --configfile 'path/to/config'
#workdir: config['data_dir'] # set working directory, a command-line way to set it: --directory 'path/to/your/dir'
config['data_dir']=str(os.getcwd())

rule all:
    input: 
        "raw_fastq/download_fastq.done"
        
rule download_fastq:
    params:
        script_dir=script_dir,
        project_name=config["project_name"],
        project_ID=config["project_ID"],
        data_path=config['data_dir']
    output: 
        touch("raw_fastq/download_fastq.done")
    run:
        shell("{params.script_dir}/download_fastq.sh '{params.project_name}' '{params.project_ID}' '{params.data_path}'")
        shell("python {params.script_dir}/get_sample_info.py --data_path {params.data_path}/raw_fastq")

    
