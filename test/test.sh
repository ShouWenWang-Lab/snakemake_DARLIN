#!/bin/bash

cd bulk
# run matlab-based DARLIN pipeline 
snakemake -s ../../snakefiles/snakefile_matlab_DARLIN_Part1.py  --configfile config.yaml --core 10 --config sbatch=0 -R DARLIN
# merge results from all samples (we have only one sample 'test' here. So, there is not much difference)
snakemake -s ../../snakefiles/snakefile_matlab_DARLIN_Part2.py --configfile config.yaml --core 5 --ri -R generate_report -R plots 

cd ../sc-10X
# run matlab-based DARLIN pipeline 
snakemake -s ../../snakefiles/snakefile_matlab_DARLIN_Part1.py  --configfile config.yaml --core 10 --config sbatch=0 -R DARLIN

cd ../scCamellia
# run matlab-based DARLIN pipeline 
snakemake -s ../../snakefiles/snakefile_matlab_DARLIN_Part1.py  --configfile config.yaml --core 10 --config sbatch=0 -R DARLIN
