#!/bin/bash

set -eu

sample_name=$1
r1_s3=$2
r2_s3=$3

# download STAR index
aws s3 cp --recursive s3://riboduct-db ~/riboduct-db

# download FASTQ
mkdir fastq
aws s3 cp $r1_s3 fastq/r1.fastq.bz2
aws s3 cp $r2_s3 fastq/r2.fastq.bz2

# create config
sed "s/SAMPLE_NAME/$sample_name/" config.yaml.template > config.yaml

eval "$(conda shell.bash hook)"  # allow script to activate conda
conda activate riboduct_env
bash riboduct.sh run config.yaml
conda deactivate

# upload result
aws s3 cp --recursive expression s3://riboduct-output/$sample_name

