# Introduction

This is a pipeline written with snakemake to automatically manage the data preprocessing (e.g., run PEAR to merge R1 and R2), sequence quality control, and run CARLIN pipeline. All you need to provide is the raw fastq file.

This pipeline is especially useful when you have multiple samples from a single sequencing run. All you have to do is to provide the relevant information in a `config.yaml` file.   


This is a brief example:
![image info](https://user-images.githubusercontent.com/4595786/205734971-e4a62308-9d16-4727-9107-36aff168a6d3.png)

This pipeline must be used together with a customized version of the CARLIN pipeline at https://github.com/ShouWenWangLab/Custom_CARLIN

## Clone the code 
First, go to a directory where you want to store the code
```bash
code_directory='your/code/directory'
cd $code_directory
```

Next, run the following bash commands to clone both the customized CARLIN repository as well as this snakemake respository.
```bash
git clone git@github.com:ShouWenWangLab/snakemake_carlin.git
mkdir CARLIN_pipeline
cd CARLIN_pipeline
git clone git@github.com:ShouWenWangLab/Custom_CARLIN.git
```

Setup the environment
```bash
conda install -n base -c conda-forge mamba --yes
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
pip install --user ipykernel
pip install jupyterlab umi_tools seaborn papermill biopython
mamba install -c bioconda fastqc multiqc --yes
python -m ipykernel install --user --name=snakemake
```

Additionally, you will need to install [pear](https://www.h-its.org/downloads/pear-academic/) and `matlab` so that they will be available as commands in terminal. We only consider running this pipeline on a remote server using a SLURM system. The matlab will be loaded with the command
```bash
module load matlab/2019a
```

## Running the pipeline

### file structure
As indicated in the above example, the `config.yaml` file should be at the root folder of a project, and the fastq data should be at the folder `raw_fastq`. 
We assume that the data is generated with Miseq machine from Illumina. Specifically, we assume that the file name starts with a sample_ID, and has both R1 and R2: 
```python
fq_R1=f"{sample}_L001_R1_001.fastq.gz"
fq_R2=f"{sample}_L001_R2_001.fastq.gz"
```
Please rename the files if they are not in this format. An example `config.ymal` file is
```yaml
project_name : 'Li_112219'
project_ID : '144505366'
script_dir: '{code_directory}/snakemake_carlin/source'
CARLIN_dir : '{code_directory}/CARLIN_pipeline'
SampleList : ['HSC','MPP','MyP'] #Remove 1_S*, it will have few reads, affect the output
cfg_type : 'sc10xV3' # available protocol: BulkRNA_Tigre_14UMI, BulkRNA_Rosa_14UMI, BulkRNA_12UMI, scLimeCat,sc10xV3
template : 'cCARLIN' # short_primer_set: {Tigre_2022_v2, Rosa_v2, cCARLIN}, long_primer_set: {Tigre_2022,Rosa,cCARLIN}
read_cutoff_UMI_override : [3,10] # assume to be a list, UMI cutoff is the same as CB cutoff for single-cell protocol
CARLIN_memory_factor : 300 # request memory at X times the size of the pear fastq file.
sbatch : 1 # 1, run sbatch job;  0, run in the interactive mode. 
CARLIN_max_run_time : 12 # hour
```
`code_directory` should be the same directory where you clone the code. 

`SampleList` should be the list of samples that you want to analyze. 

`cfg_type` should match the protocol of the experiment. Some of the provided protocols include
 
 * `BulkRNA_Tigre_14UMI`: Bulk CARLIN library with Tigre locus, with a UMI of 14bp
 * `BulkRNA_Rosa_14UMI`:  Bulk CARLIN library with Rosa locus, with a UMI of 14bp
 * `BulkRNA_12UMI`: Bulk CARLIN library with Col1a1 locus, with a UMI of 12bp
 * `scLimeCat`: Single-cell CARLIN library using the Camellia-seq protocol
 * `sc10xV3`: Single-cell CARLIN library using the 10X v3 protocol
 
 `template` should match the primer set used. We have template corresponding to shorter primers in TC and RC: {`Tigre_2022_v2`, `Rosa_v2`}, and longer primers: {`Tigre_2022`, `Rosa`}. For Col1a1 locus, we only have a single primer set, corresponding to tempalte `cCARLIN`. 

`read_cutoff_UMI_override`: minimum number of reads needed to support a UMI (bulk library) or a cell barcode (single cell library). It should be a list of read cutoff like [3,10].

`CARLIN_memory_factor`: When running on o2, the requested memory should be `CARLIN_memory_factor` times the fastq file size.

`sbatch`: when running on o2, whether to run with sbatch jobs (1) or in interactive mode (0)

`CARLIN_max_run_time`: When running on o2, the maximum run time to request, in the unit of hours

### Getting data from base space
When the fastq files are not downloaded yet in the `raw_fastq` folder, and the data sits at base space of Illumina, you can provide `project_name` and `project_ID` in `config.yaml` to automaically download the data. 

First, check the available fastq data with the terminal command 
```bash
bs auth # this needs to be done only once for authentification
bs list project
```
![image info](https://user-images.githubusercontent.com/4595786/205739547-7439adf6-90a3-45bc-ac36-71758cef4e6c.png)

Next, select the desired project name and ID. In the above `config.yaml` file, we selected the data from the first entry. 

Then, run the snakemake script at the same directory as the `config.yaml` file
```bash
snakemake -s $code_directory/snakemake_carlin/snakefile_get_data.py --configfile config.yaml --core 1
```

### CARLIN analsysis
This command will generate the QC report and process each sample with the CARLIN pipeline
```bash
snakemake -s $code_directory/packages/snakemake_carlin/snakefile_integrate_CARLIN.py  --configfile config.yaml --core 10
```

Finally, you may run this command to get an html report across all samples
```bash
snakemake -s $code_directory/packages/snakemake_carlin/snakefile_downstream_fast.py --configfile config.yaml --core 5 --ri -R generate_report -R plots 
```
The result will show up at the `merge_all` folder as shown in the above image. 


### A single-cell pipeline for libraries with higher amplification heterogneity
We also developed our own CARLIN pipeline that works well for single-cell libraries with higher amplification heterogeneity, e.g., one cell gets 10K reads, while another cell only has 10 reads. This pipeline is written in jupyter notebook (`source/single_cell_CARLIN.ipynb`), and it requires to first install a companion repository `carlin_hf`. 
```bash
cd $code_directory
git clone git@github.com:ShouWenWangLab/carlin_hf.git
cd carlin_hf
python setup.py develop
```

Then, you need to add more parameters in the `config.yaml` file
```yaml
project_name : 'Li_112219'
project_ID : '144505366'
script_dir: '{code_directory}/snakemake_carlin/source'
CARLIN_dir : '{code_directory}/CARLIN_pipeline'
SampleList : ['HSC','MPP','MyP'] #Remove 1_S*, it will have few reads, affect the output
cfg_type : 'sc10xV3' # available protocol: BulkRNA_Tigre_14UMI, BulkRNA_Rosa_14UMI, BulkRNA_12UMI, scLimeCat,sc10xV3
template : 'cCARLIN' # short_primer_set: {Tigre_2022_v2, Rosa_v2, cCARLIN}, long_primer_set: {Tigre_2022,Rosa,cCARLIN}
read_cutoff_UMI_override : [3,10] # assume to be a list, UMI cutoff is the same as CB cutoff for single-cell protocol
CARLIN_memory_factor : 300 # request memory at X times the size of the pear fastq file.
sbatch : 1 # 1, run sbatch job;  0, run in the interactive mode. 
CARLIN_max_run_time : 12 # hour
single_cell_pipeline: # This is an extension, needed only if you run snakefile_single_cell_CARLIN.py
    coarse_grained_readcutoff_floor: 5 # the lower bound of the later read count filtering, after denoising, and re-group reads. 
    distance_relative_threshold: 0.03 # 5% error rate, will be multipled with the sequence length
    read_ratio_threshold: 0.6
    seq_3prime_upper_N: 15
    output_folder: 'custom_method'
```

Finally, run 
```bash
snakemake -s $code_directory/packages/snakemake_carlin/snakefile_single_cell_CARLIN.py  --configfile config.yaml --core 10
```

The result will show up as a jupyter notebook and a corresponding html report
![image](https://user-images.githubusercontent.com/4595786/205761409-2f5678c2-51aa-409b-93f1-ab32509a2c74.png)
