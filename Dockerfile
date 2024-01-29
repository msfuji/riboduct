FROM ubuntu:22.04

# avoid user interaction during apt update
ARG DEBIAN_FRONTEND=noninteractive

# Install AWS CLI
RUN apt update && \
 apt install -y awscli

# Install miniconda
RUN apt install -y wget
RUN mkdir -p ~/miniconda3
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
RUN bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
RUN rm -rf ~/miniconda3/miniconda.sh
RUN ~/miniconda3/bin/conda init bash
ENV PATH ~/miniconda3/bin:$PATH

# Copy pipeline files
RUN mkdir riboduct
WORKDIR riboduct
COPY . .

# Install riboduct
RUN bash riboduct.sh install riboduct_env

CMD bash run.sh sample_A s3://riboduct-input/sample_A.R1.fastq.bz2 s3://riboduct-input/sample_A.R2.fastq.bz2
