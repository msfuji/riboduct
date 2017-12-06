rule download_genome:
    output:
        config["db_dir"]+"/genome/hs37d5.fa.gz"
    log:
        config["db_dir"]+"/genome/download_genome.log"
    shell:
        "curl -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz && "
        "mv hs37d5.fa.gz {output}"

rule decompress_genome:
    input:
        config["db_dir"]+"/genome/hs37d5.fa.gz"
    output:
        config["db_dir"]+"/genome/hs37d5.fa"
    log:
        config["db_dir"]+"/genome/decompress_genome.log"
    shell:
        "gunzip -c {input} > {output}"

rule download_gtf:
    output:
        config["db_dir"]+"/gene_model/gencode.v19.annotation.gtf.gz"
    log:
        config["db_dir"]+"/gene_model/download_gtf.log"
    shell:
        "curl -O ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz && "
        "mv gencode.v19.annotation.gtf.gz {output}"

rule format_gtf:
    input:
        config["db_dir"]+"/gene_model/gencode.v19.annotation.gtf.gz"
    output:
        config["db_dir"]+"/gene_model/gencode.v19.annotation.hs37d5_chr.gtf"
    log:
        config["db_dir"]+"/gene_model/format_gtf.log"
    shell:
        "gunzip -c {input} | tail -n +6 | sed -e \"s/^chrM/MT/g;s/^chr//g\" > {output}"

rule setup_db:
    input:
        genome=config["db_dir"]+"/genome/hs37d5.fa",
        gtf=config["db_dir"]+"/gene_model/gencode.v19.annotation.hs37d5_chr.gtf"
    output:
        config["db_dir"]+"/star_index/SAindex",
        dir=config["db_dir"]+"/star_index"
    threads: 8
    log:
        config["db_dir"]+"/star_index/setup_db.log"
    shell:
        "rm -f {output.dir}/* && "
        "star "
        "--runMode genomeGenerate "
        "--genomeDir {output.dir} "
        "--genomeFastaFiles {input.genome} "
        "--sjdbOverhang 125 "
        "--sjdbGTFfile {input.gtf} "
        "--runThreadN {threads}"
