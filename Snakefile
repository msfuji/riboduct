rule download_genome:
    output:
        "reference/genome/hs37d5.fa.gz"
    shell:
        "wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"

rule download_gene_model:
    output:
        "reference/gene_model/gencode.v19.annotation.gtf.gz"
    shell:
        "wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"

rule setup_db:
