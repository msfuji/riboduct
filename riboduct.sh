#!/bin/bash

command=$1
param=$2

#
# install riboduct
#
if [ $command = "install" ]; then
  conda_env=$param
  echo "Installing riboduct..."
  conda create -n $conda_env -y python=3.10 conda-forge::mamba
  eval "$(conda shell.bash hook)"  # allow script to activate conda
  conda activate $conda_env
  mamba install -y -c bioconda -c conda-forge -c r snakemake-minimal star subread samtools picard rna-seqc r r-dplyr r-data.table r-readr pandas bx-python
  conda deactivate
  echo "DONE."
  exit 0
elif [ -z $param ]; then
  echo "[ERROR] config.yaml not specified."
  exit 1
fi

config_yaml=$param

#
# check config.yaml whether to use SGE
#
sge_jobs=`awk '$1~/^sge_jobs:/ {print $2}' $config_yaml`
if [ -z $sge_jobs ]; then
  use_sge=false
  echo "Use SGE cluster... [No]"
else
  use_sge=true
  echo "Use SGE cluster... [Yes]"
  echo "Max SGE jobs... ["$sge_jobs"]"
fi

#
# add SGE options for snakemake
#
snake_command="snakemake --config env_dir=$CONDA_PREFIX --configfile $config_yaml --cores 4"
if $use_sge; then
  snake_command=`echo $snake_command --cluster \"qsub -terse -cwd -pe def_slot {threads} -l s_vmem={params.memory},mem_req={params.memory} -o {log} -e {log}\" --jobs $sge_jobs`
fi

#
# run snakemake
#
if [ $command = "index" ]; then
  echo "Indexing database..."
  bash -c "$snake_command index"
  echo "DONE."
elif [ $command = "run" ]; then
  echo "Running riboduct..."
  bash -c "$snake_command"
  echo "DONE."
else
  echo "[ERROR] Unknown command ", $command
fi
