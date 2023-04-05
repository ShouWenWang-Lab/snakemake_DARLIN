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
# Warining: the data_path should be at the path/to/raw_fastq (i.e. under raw_fastq)
# We assume that under raw_fastq, there are a list fo sample folder
# under each sample folder there will be R1 and R2 files
# we will rename files

# To put them at the level of raw_fastq, and remove the folders
#
# os.system(f'mv */*.fastq.gz temp_dir; rm -r {sample}')
# os.system('mv temp_dir/* .; rm -r temp_dir')
######

parser.add_argument(
    "--suffix",
    type=str,
    default=".fastq.gz",
    help="suffix of the file",
)


args = parser.parse_args()


all_files = sorted(
    Path(".").glob(
        os.path.join(
            f"*{args.suffix}",
        )
    )
)


for file in all_files:
    file = str(file)
    print(file)
    base_file = file
    condition_1 =  (
        f"_R1_001{args.suffix}" in base_file
    )
    condition_2 =  (
        f"_R2_001{args.suffix}" in base_file
    )
    
    if condition_1 and ('_S' in base_file):
        new_sample=base_file.split('_S')[0]
        os.system(f"mv {base_file} {new_sample}_L001_R1_001.fastq.gz")
    elif condition_2 and ('_S' in base_file):
        new_sample=base_file.split('_S')[0]
        os.system(f"mv {base_file} {new_sample}_L001_R2_001.fastq.gz")
    else:
        print(
            "Error. Please place the command at the one level upper directory above samples (typically under raw_fastq/"
        )
