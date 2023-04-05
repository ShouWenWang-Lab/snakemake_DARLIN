import os
from pathlib import Path
import os
import pandas as pd
import argparse
from os import listdir
from os.path import isfile, join
# parse cmd line arguments
parser = argparse.ArgumentParser(description="rename files")

######
# Warining: the data_path should be at the path/to/raw_fastq
# We assume that under raw_fastq, there are a list fo sample folder
# under each sample folder there will be R1 and R2 files
# we will rename files

# To put them at the level of raw_fastq, and remove the folders
#         
    # os.system(f'mv */*.fastq.gz temp_dir; rm -r {sample}')    
    # os.system('mv temp_dir/* .; rm -r temp_dir')
######



parser.add_argument(
    "--data_path",
    type=str,
    default=".",
    help="Path/to/raw_fastq",
)
parser.add_argument(
    "--suffix",
    type=str,
    default=".fq.gz",
    help="suffix of the file",
)


args=parser.parse_args()

args.data_path=os.path.abspath(args.data_path)
all_files = sorted(
    Path(args.data_path).glob(
        os.path.join(
            "*",
            f"*{args.suffix}",
        )
    )
)


os.chdir(args.data_path)
os.system(f'mkdir -p temp_dir')
print(f'Current path: {args.data_path}')
for file in all_files:
    print(file)
    file=str(file)
    print(file)
    sample=file.split('/')[-2]
    print(f'Sample: {sample}')
    if (f'_1{args.suffix}' in file) or (f'_R1{args.suffix}' in file):
        os.system(f'mv {file} {sample}/{sample}_L001_R1_001.fastq.gz')
    if (f'_2{args.suffix}' in file) or (f'_R2{args.suffix}' in file):
        os.system(f'mv {file} {sample}/{sample}_L001_R2_001.fastq.gz')
