# riboduct
RNA-Seq pipeline.

## Installation
### Install dependencies
`riboduct` requires `conda` package manager and `Snakemake` workflow managemer.
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

### Set up riboduct
Download riboduct and set up the conda environment.
```
git clone https://github.com/msfuji/riboduct.git
cd riboduct
conda env create --name riboduct --file environment.yaml
```

Download reference human genome (hs37d5) and gene model (GENCODE v19).
```
source activate riboduct
snakemake setup_db
source deactivate
```
If you are using `pyenv`, `activate` may conflict. In that case, use `conda info -e` to find the full path of `activate`.

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
