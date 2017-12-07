bin_dir=config["env_dir"]+"/bin/"

rule all:
    input:
        "expression/raw_counts.gmt"
#         "expression/fpkm.gmt"
#         "expression/fpkm_uq.gmt"

################################################################################

rule download_genome:
    output:
        config["db_dir"]+"/genome/hs37d5.fa.gz"
    log:
        config["db_dir"]+"/genome/download_genome.log/"
    shell:
        "curl -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz && "
        "mv hs37d5.fa.gz {output}"

rule decompress_genome:
    input:
        config["db_dir"]+"/genome/hs37d5.fa.gz"
    output:
        config["db_dir"]+"/genome/hs37d5.fa"
    log:
        config["db_dir"]+"/genome/decompress_genome.log/"
    shell:
        "gunzip -c {input} > {output}"

rule download_gtf:
    output:
        config["db_dir"]+"/gene_model/gencode.v19.annotation.gtf.gz"
    log:
        config["db_dir"]+"/gene_model/download_gtf.log/"
    shell:
        "curl -O ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz && "
        "mv gencode.v19.annotation.gtf.gz {output}"

rule format_gtf:
    input:
        config["db_dir"]+"/gene_model/gencode.v19.annotation.gtf.gz"
    output:
        config["db_dir"]+"/gene_model/gencode.v19.annotation.hs37d5_chr.gtf"
    log:
        config["db_dir"]+"/gene_model/format_gtf.log/"
    shell:
        "gunzip -c {input} | tail -n +6 | sed -e \"s/^chrM/MT/g;s/^chr//g\" > {output}"

rule setup_db:
    input:
        genome=config["db_dir"]+"/genome/hs37d5.fa",
        gtf=config["db_dir"]+"/gene_model/gencode.v19.annotation.hs37d5_chr.gtf"
    output:
        config["db_dir"]+"/star_index/SAindex",
        dir=config["db_dir"]+"/star_index/"
    log:
        config["db_dir"]+"/star_index/setup_db.log/"
    threads: 8
    shell:
        bin_dir+"STAR "
        "--runMode genomeGenerate "
        "--genomeDir {output.dir} "
        "--genomeFastaFiles {input.genome} "
        "--sjdbOverhang 100 "
        "--sjdbGTFfile {input.gtf} "
        "--runThreadN {threads} "
        "--outFileNamePrefix {output.dir}"

################################################################################

rule star_1_pass:
    input:
        lambda wildcards: config["fastq"][wildcards.sample],
        index=config["db_dir"]+"/star_index/"
    output:
        "star_1_pass/{sample_id}/SJ.out.tab"
    shell:
        "touch {output}"

rule star_1_pass_index:
    input:
        "star_1_pass/{sample_id}/SJ.out.tab"
    output:
        "star_1_pass/index/SAindex",
    shell:
        "touch {output}"

rule star_2_pass:
    input:
        "star_1_pass/index"
    output:
        "star_2_pass/{sample_id}/Aligned.sortedByCoord.out.bam"
    shell:
        "touch {output}"

################################################################################
#
# TODO
#
rule raw_counts:
    input:
        "star_2_pass/{sample_id}/Aligned.sortedByCoord.out.bam"
    output:
        "expression/raw_counts/{sample_id}.txt"
    shell:
        "touch {output}"

rule merge_raw_counts:
    input:
        expand("expression/raw_counts/{sample_id}.txt", sample_id=config["fastq"])
    output:
        "expression/raw_counts.gmt"
    shell:
        "touch ${output}"
