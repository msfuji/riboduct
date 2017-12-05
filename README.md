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
Download by
```
git https://github.com/msfuji/riboduct.git
```
Set up by
```
cd riboduct
conda env create --name riboduct --file environment.yaml
```
### Starting riboduct
Activate the conda environment by
```
source activate riboduct
```
If you are using `pyenv`, the above command may conflict. In that case, use `conda info -e` to find the full path of `activate` command.

### Finishing riboduct
```
source deactivate
```
