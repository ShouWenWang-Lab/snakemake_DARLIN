#!/bin/bash

cd bulk
# run matlab-based DARLIN pipeline 
snakemake -s ../../snakefiles/snakefile_matlab_DARLIN_Part1.py  --configfile config.yaml --core 10 --config sbatch=0 -R DARLIN
# merge results from all samples (we have only one sample 'test' here. So, there is not much difference)
snakemake -s ../../snakefiles/snakefile_matlab_DARLIN_Part2.py --configfile config.yaml --core 5 --ri -R generate_report -R plots 

cd ../sc-10X
# run matlab-based DARLIN pipeline 
snakemake -s ../../snakefiles/snakefile_matlab_DARLIN_Part1.py  --configfile config.yaml --core 10 --config sbatch=0 -R DARLIN

# run customized python pipeline on single-cell data
snakemake -s ../../snakefiles/snakefile_python_DARLIN.py  --configfile config.yaml --core 5 --config sbatch=0 -R DARLIN

# call alleles via Matlab-based DARLIN pipeline (fast)
snakemake -s ../../snakefiles/snakefile_call_mut.py --configfile config.yaml --core 5 --config sbatch=0

cd ../scCamellia
# run matlab-based DARLIN pipeline 
snakemake -s ../../snakefiles/snakefile_matlab_DARLIN_Part1.py  --configfile config.yaml --core 10 --config sbatch=0 -R DARLIN

# call alleles via Matlab-based DARLIN pipeline (fast)
snakemake -s ../../snakefiles/snakefile_call_mut.py --configfile config.yaml --core 5 --config sbatch=0

# run customized python pipeline on single-cell data
snakemake -s ../../snakefiles/snakefile_python_DARLIN.py  --configfile config.yaml --core 5 --config sbatch=0 -R DARLIN
