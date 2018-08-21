#!bin/bash

command=$1
config_yaml=$2

#
# install riboduct
#
if [ $command = "install" ]; then
  echo "Installing riboduct..."
  conda config --add channels r
  conda config --add channels defaults
  conda config --add channels conda-forge
  conda config --add channels bioconda
  conda env create --name riboduct --file environment.yaml
  echo "DONE."
  exit 0
fi

#
# check config.yaml whether to use SGE
#
sge_jobs=`awk '$1~/^sge_jobs:/ {print $2}' config.yaml`
if [ -z $sge_jobs ]; then
  use_sge=false
  echo "Use SGE cluster... [No]"
else
  use_sge=true
  echo "Use SGE cluster... [Yes]"
  echo "Max SGE jobs... ["$sge_jobs"]"
fi

#
# activate conda
#
source activate riboduct

#
# add SGE options for snakemake
#
snake_command="snakemake --config env_dir=$CONDA_PREFIX --configfile $config_yaml"
if $use_sge; then
  snake_command=`echo $snake_command --cluster \"qsub -pe def_slot {threads} -o {log} -e {log}\" --jobs $sge_jobs`
fi

#
# run snakemake
#
if [ $command = "index" ]; then
  echo "Indexing database..."
  $snake_command index
  echo "DONE."
elif [ $command = "run" ]; then
  echo "Running riboduct..."
  $snake_command
  echo "DONE."
else
  echo "[ERROR] Unknown command ", $command
fi

#
# deactivate conda
#
source deactivate
