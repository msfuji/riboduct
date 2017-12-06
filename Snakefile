rule download_genome:
    output:
        "reference/genome/hs37d5.fa.gz"
    log:
        "reference/genome/download_genome.log"
    shell:
        "curl -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz 2> {log} && "
        "mv hs37d5.fa.gz {output}"

rule decompress_genome:
    input:
        "reference/genome/hs37d5.fa.gz"
    output:
        "reference/genome/hs37d5.fa"
    log:
        "reference/genome/decompress_genome.log"
    shell:
        "gunzip -c {input} > {output} 2> {log}"

rule download_gtf:
    output:
        "reference/gene_model/gencode.v19.annotation.gtf.gz"
    log:
        "reference/gene_model/download_gtf.log"
    shell:
        "curl -O ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz 2> {log} && "
        "mv gencode.v19.annotation.gtf.gz {output}"

rule format_gtf:
    input:
        "reference/gene_model/gencode.v19.annotation.gtf.gz"
    output:
        "reference/gene_model/gencode.v19.annotation.hs37d5_chr.gtf"
    log:
        "reference/gene_model/format_gtf.log"
    shell:
        "gunzip -c {input} | tail -n +6 | sed -e \"s/^chrM/MT/g;s/^chr//g\" > {output} 2> {log}"

rule setup_db:
    input:
        genome="reference/genome/hs37d5.fa",
        gtf="reference/gene_model/gencode.v19.annotation.hs37d5_chr.gtf"
    output:
        "reference/star_index/SAindex",
        dir="reference/star_index"
    threads: 8
    log:
        "reference/star_index/setup_db.log"
    shell:
        "rm -f {output.dir}/* && "
        "star "
        "--runMode genomeGenerate "
        "--genomeDir {output.dir} "
        "--genomeFastaFiles {input.genome} "
        "--sjdbOverhang 125 "
        "--sjdbGTFfile {input.gtf} "
        "--runThreadN {threads} 2> {log}"
