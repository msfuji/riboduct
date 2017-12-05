rule download_genome:
    output:
        "reference/genome/hs37d5.fa.gz"
    shell:
        "wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"

rule download_gtf:
    output:
        "reference/gene_model/gencode.v19.annotation.gtf.gz"
    shell:
        "wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"

rule format_gtf:
    input:
        "reference/gene_model/gencode.v19.annotation.gtf.gz"
    output:
        "reference/gene_model/gencode.v19.annotation.hs37d5_chr.gtf"
    shell:
        "gunzip -c {input} | tail -n +6 | sed -e \"s/^chrM/MT/g;s/^chr//g\" > {output}"

rule setup_db:
    input:
        genome="reference/genome/hs37d5.fa.gz",
        gtf="reference/gene_model/gencode.v19.annotation.hs37d5_chr.gtf"
    output:
        "reference/star_index/SAindex"
    threads: 4
    shell:
        "dir=`dirname {output}`; rm -f $dir/* &&"
        "star"
        "--runMode genomeGenerate"
        "--genomeDir {output}"
        "--genomeFastaFiles {genome}"
        "--sjdbOverhang 125"
        "--sjdbGTFfile {gtf}"
        "--runThreadN {threads}"
