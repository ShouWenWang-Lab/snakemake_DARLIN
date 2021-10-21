In order to run Snakemake pipeline, you can do the following.

1, To grab the new data from basespace, first view sample information

```bash
 bs list project
```

2, Make a folder under `DATA` with a proper name for this dataset, and create a YAML file uner this folder with the following sample information

```yaml
project_name : '3_pulses_LL605_LL607'
project_ID : '295673483'
script_dir: '/n/data1/bch/hemonc/camargo/li/my_scripts/snakeTemplate'
CARLIN_dir : '/n/data1/bch/hemonc/camargo/li/CARLIN_pipeline'
SampleList : [] # [] means no selection, and all will be selected;
cfg_type : 'BulkRNA' # BulkDNA_Tigre (with a UMI QC threshold 25), BulkDNA (UMI QC threshold 30);  
template : 'cCARLIN' # (for Tigre, the template is switched)
read_cutoff_override : 3 
read_cutoff_floor : 1  # This variable is overrided by read_cutoff_override
CARLIN_memory_factor : 100 # request memory at X times the size of the pear fastq file.
sbatch : 1 # 1, run sbatch job;  0, run in the interactive mode.  This only affects the CARLIN analysis and csv generation. 
CARLIN_max_run_time : 12 # hours
```

3, Run snakemake pipeline sequentially `within` this folder. Make sure that each snakemake run is finished before proceeding to the next one. 

First, we have pre-defined the snakefile directory in the system:
```bash
snake=/n/data1/bch/hemonc/camargo/li/snakeTemplate
```


Get the data

```bash
snakemake -s ${snake}/Snakefile_get_data.py --core 1
```

Run pear, fastqc, multiqc, and CARLIN analysis. Note that only the CARLIN analysis has the option to use the sbatch job mode. This is because we need to make sure that the pear result is done before we proceed with the CARLIN analysis. You can use more cores to parallelize the analysis (make sure that you run `srun -t 0-12:00 --pty -p interactive -n 10 --mem 50G /bin/bash` to get more resources first).

```bash
snakemake -s ${snake}/Snakefile_QC_CARLIN.py --core 10
```

Merge all samples, calculate allele bank, make relevant plots, and transfer data to `o2_data`, where you can git push to share data with your local machine.

```bash
snakemake -s ${snake}/Snakefile_downstream.py --core 1
```

Sometimes, `o2_data` can get too large due to the historical files. You can clean up the `.git` there by running 

```bash
sh ${snake}/clean_o2_data.sh 
```
Note that, to run this, you will need to switch to the login node. 

