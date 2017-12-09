# riboduct
RNA-Seq pipeline.

## Installation
### Install conda
`riboduct` requires `conda` package manager. To install `conda` for Linux,
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### Set up pipeline
Download `riboduct` and set up the `conda` environment.
```
git clone https://github.com/msfuji/riboduct.git
cd riboduct
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda env create --name riboduct --file environment.yaml
```
Download RNA-SeQC.
```
mkdir bin
wget -P bin/ http://www.broadinstitute.org/cancer/cga/tools/rnaseqc/RNA-SeQC_v1.1.8.jar
```
Download reference human genome (hs37d5) and gene model (GENCODE v19).
```
# edit db_dir
vi config.yaml

source activate riboduct

# Download genome and gtf to db_dir, and index them. This may take a long time.
snakemake --configfile config.yaml --config env_dir=$CONDA_PREFIX setup_db

## Speed up indexing on SGE clusters.
# snakemake --configfile config.yaml --config env_dir=$CONDA_PREFIX \
# --cluster "qsub -pe def_slot {threads} -o {log} -e {log}" --jobs 2 setup_db

source deactivate
```

## Usage
```
# edit FASTQ path
vi config.yaml

# activate the conda environment
source activate riboduct

# run
snakemake --configfile config.yaml --config env_dir=$CONDA_PREFIX

## run on SGE clusters.
# snakemake --configfile config.yaml --config env_dir=$CONDA_PREFIX \
# --cluster "qsub -pe def_slot {threads} -cwd -o {log} -e {log}" --jobs 100

# finish
source deactivate
```
