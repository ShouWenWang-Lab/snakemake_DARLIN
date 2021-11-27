import os
import pandas as pd
import argparse
from os import listdir
from os.path import isfile, join
# parse cmd line arguments
parser = argparse.ArgumentParser(description="Download data")

parser.add_argument(
    "--data_path",
    type=str,
    default=".",
    help="Path of the raw fastq files",
)

data_path = parser.parse_args().data_path


## step two, get sample information
SampleList = []
SampleID = []
SampleShortName = []
onlyfiles = [f for f in listdir(data_path) if isfile(join(data_path, f))]
for file in onlyfiles:
    if '_L001_R1_001.fastq.gz' in file:
        sample_name = file.split('/')[-1].split('_L001_R1')[0]
        SampleList.append(sample_name)
        SampleShortName.append(sample_name.split('_')[0])
        SampleID.append(sample_name.split('_')[1])
        
df=pd.DataFrame({'Sample_long':SampleList,'Sample_ID':SampleID,'Sample_short':SampleShortName})     
df.to_csv(f'{data_path}/sample_info.csv')
SampleString="["
for xx in SampleList:
    SampleString += "\'"+xx+"\',"
SampleString=SampleString[:-1]+"]"
print(f'SampleList: {SampleString}')


