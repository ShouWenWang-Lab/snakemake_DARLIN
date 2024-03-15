#!/bin/bash

cd bulk_cCARLIN
# run CARLIN pipeline (based on matlab)
snakemake -s ../../snakefiles/snakefile_integrate_CARLIN.py  --configfile config.yaml --core 10 --config sbatch=0 -R CARLIN

cd ../sc-10X
# run CARLIN pipeline (based on matlab)
snakemake -s ../../snakefiles/snakefile_integrate_CARLIN.py  --configfile config.yaml --core 10 --config sbatch=0 -R CARLIN
# run customized python pipeline on single-cell data
snakemake -s ../../snakefiles/snakefile_single_cell_DARLIN.py  --configfile config.yaml --core 5 --config sbatch=0 -R CARLIN


cd ../scCamellia
# run CARLIN pipeline (based on matlab)
snakemake -s ../../snakefiles/snakefile_integrate_CARLIN.py  --configfile config.yaml --core 10 --config sbatch=0 -R CARLIN
# run customized python pipeline on single-cell data
snakemake -s ../../snakefiles/snakefile_single_cell_DARLIN.py  --configfile config.yaml --core 5 --config sbatch=0 -R CARLIN
