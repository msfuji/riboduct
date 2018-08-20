#!bin/bash

command=$1
config_yaml=$2


# install
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda env create --name riboduct --file environment.yaml

# index
source activate riboduct
snakemake --configfile $config_yaml --config env_dir=$CONDA_PREFIX setup_db
source deactivate

# index_sge
source activate riboduct
snakemake --configfile $config_yaml --config env_dir=$CONDA_PREFIX \
--cluster "qsub -pe def_slot {threads} -o {log} -e {log}" --jobs 2 setup_db
source deactivate

# run
source activate riboduct
snakemake --configfile $config_yaml --config env_dir=$CONDA_PREFIX
source deactivate

# run_sge
source activate riboduct
snakemake --configfile $config_yaml --config env_dir=$CONDA_PREFIX \
--cluster "qsub -pe def_slot {threads} -cwd -o {log} -e {log}" \
--jobs $max_sge_jobs
source deactivate
