# riboduct
RNA-Seq pipeline.

## Installation
### Install dependencies
`riboduct` requires `conda` package manager and `Snakemake` workflow manageer.
To install `conda` for Linux,
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
For MacOSX,
```
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```
To install `Snakemake`,
```
conda install -c bioconda snakemake
```

### Set up pipeline
Download `riboduct` and set up the `conda` environment.
```
git clone https://github.com/msfuji/riboduct.git
cd riboduct
conda env create --name riboduct --file environment.yaml
```

Download reference human genome (hs37d5) and gene model (GENCODE v19).
```
# Activate the conda environment.
source activate riboduct

# Download genome, gtf and index them. This may take a long time.
snakemake setup_db
# Speed up indexing when you are working on clusters.
# snakemake --cluster "qsub -pe def_slot {threads} setup_db"

source deactivate
```

## Usage

### Starting riboduct
Activate the conda environment by
```
source activate riboduct
```

### Finishing riboduct
```
source deactivate
```
