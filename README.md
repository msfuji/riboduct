# riboduct
RNA-Seq pipeline. `riboduct` maps sequence reads onto reference genome, counts
reads on annotated transcripts, computes FPKM and FPKM-UQ. Also computes
several metrics for quality control.

<img src="misc/riboduct_flowchart.png" width="300px">

## Installation
### 1. Install conda
`riboduct` requires `conda` package manager. To install `conda` for Linux,
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### 2. Install pipeline
Download and install `riboduct`. Here you can replace `RIBODUCT_ENV` with your favorite name.
```
git clone https://github.com/msfuji/riboduct.git
cd riboduct
bash riboduct.sh install RIBODUCT_ENV
```

### 3. Set up database
First, download reference genome and gene annotation to anywhere you want.
```
# I primarily use the human genome (GRCh37).
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

# Alternatively, you may use the human genome GRCh38.
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.p12.genome.fa.gz
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz

# The mouse genome GRCm38 can be also used.
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/GRCm38.p6.genome.fa.gz
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/gencode.vM18.annotation.gtf.gz
```
Make a database directory, and modify `config.yaml` so that the `db_dir`
parameter points to the database directory. Also modify two other parameters
(`genome_fa_gz`, `annotation_gtf_gz`) for the downloaded files of genome and
gene annotation. Start indexing of database.
```
conda activate RIBODUCT_ENV
bash riboduct.sh index config.yaml
conda deactivate
```

## Usage
Make a local copy of pipeline for each project.
```
git clone https://github.com/msfuji/riboduct.git
cd riboduct
```
Modify `config.yaml`. First, add path of FASTQ files you want to analyze.
Also modify `db_dir` to point the directory where you installed database.
When finished, start running the pipeline.
```
conda activate RIBODUCT_ENV
bash riboduct.sh run config.yaml
conda deactivate 
```
