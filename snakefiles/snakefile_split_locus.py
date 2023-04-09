import pandas as pd
from pathlib import Path
import seaborn as sns
from matplotlib import pyplot as plt 
import os
import sys
from darlin import help_functions as hf
from darlin.settings import script_dir, QC_dir, ref_dir, CARLIN_dir
#configfile: "config.yaml"  # command line way to set it: --configfile 'path/to/config'
#workdir: config['data_dir'] # set working directory, a command-line way to set it: --directory 'path/to/your/dir'
config['data_dir']=str(os.getcwd())
script_dir=script_dir
##################
## preprocessing
################## 

print(f'Current work dir: {os.getcwd()}')
if len(config['SampleList'])==0: 
    df=pd.read_csv('raw_fastq/sample_info.csv')
    SampleList=list(df['Sample_long'])
else:
    SampleList=config['SampleList']

        
##################
## start the rules
################## 
rule all:
    run:
        #wk_dir=config['data_dir']+'/pear_output'
        df_list=[]
        for sample in SampleList:
            os.system(f"sh {script_dir}/split_locus.sh  {sample} pear_output")
            df_temp=pd.read_csv(f'pear_output/{sample}.split_locus.tsv',header=None,names=['locus','data','sample'])
            df_list.append(df_temp)
        df_data=pd.concat(df_list)
        df_data['reads']=df_data['data'].apply(lambda x:x.split(' ')[0])
        df_data=df_data.drop(columns='data')
        df_data['reads']=df_data['reads'].astype(int)
        df_data.to_csv('pear_output/all.split_locus.csv')
        sns.factorplot(data=df_data,x='sample',y='reads',hue='locus',kind='bar')
        plt.xticks(rotation=90)
        plt.xlabel('')
        plt.savefig(f'split_locus.png')

        SampleString="["
        for xx in SampleList:
            SampleString += "\'"+xx+f".Tigre\',"
        SampleString=SampleString[:-1]+"]"
        print(f'SampleList: {SampleString}')

        SampleString="["
        for xx in SampleList:
            SampleString += "\'"+xx+f".cCARLIN\',"
        SampleString=SampleString[:-1]+"]"
        print(f'SampleList: {SampleString}')

        SampleString="["
        for xx in SampleList:
            SampleString += "\'"+xx+f".Rosa\',"
        SampleString=SampleString[:-1]+"]"
        print(f'SampleList: {SampleString}')
            
    