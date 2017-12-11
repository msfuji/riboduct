bin_dir=config["env_dir"]+"/bin/"

rule all:
    input:
        "expression/raw_counts.tsv",
        "expression/fpkm.tsv",
        "expression/fpkm_uq.tsv",
        "qc/rna_seqc/index.html"

################################################################################
#
# database construction
#

rule setup_db:
    input:
        config["db_dir"]+"/genome/hs37d5.fa.fai",
        config["db_dir"]+"/genome/hs37d5.dict",
        config["db_dir"]+"/star_index/SAindex"

rule download_genome:
    output:
        config["db_dir"]+"/genome/hs37d5.fa.gz"
    log:
        config["db_dir"]+"/log/download_genome/"
    shell:
        "wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz && "
        "mv hs37d5.fa.gz {output}"

rule decompress_genome:
    input:
        config["db_dir"]+"/genome/hs37d5.fa.gz"
    output:
        config["db_dir"]+"/genome/hs37d5.fa"
    log:
        config["db_dir"]+"/log/decompress_genome/"
    shell:
        "gunzip -c {input} > {output}"

rule faidx_genome:
    input:
        config["db_dir"]+"/genome/hs37d5.fa"
    output:
        config["db_dir"]+"/genome/hs37d5.fa.fai"
    log:
        config["db_dir"]+"/log/faidx_genome/"
    shell:
        bin_dir+"samtools faidx {input}"

rule dict_genome:
    input:
        config["db_dir"]+"/genome/hs37d5.fa"
    output:
        config["db_dir"]+"/genome/hs37d5.dict"
    log:
        config["db_dir"]+"/log/dict_genome/"
    shell:
        bin_dir+"picard CreateSequenceDictionary "
        "R={input} "
        "O={output}"

rule download_gtf:
    output:
        config["db_dir"]+"/gene_model/gencode.v19.annotation.gtf.gz"
    log:
        config["db_dir"]+"/log/download_gtf/"
    shell:
        "wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz && "
        "mv gencode.v19.annotation.gtf.gz {output}"

rule format_gtf:
    input:
        config["db_dir"]+"/gene_model/gencode.v19.annotation.gtf.gz"
    output:
        config["db_dir"]+"/gene_model/gencode.v19.annotation.hs37d5_chr.gtf"
    log:
        config["db_dir"]+"/log/format_gtf/"
    shell:
        "gunzip -c {input} | tail -n +6 | sed -e \"s/^chrM/MT/g;s/^chr//g\" > {output}"

rule star_index:
    input:
        genome=config["db_dir"]+"/genome/hs37d5.fa",
        gtf=config["db_dir"]+"/gene_model/gencode.v19.annotation.hs37d5_chr.gtf"
    output:
        config["db_dir"]+"/star_index/SAindex",
        dir=config["db_dir"]+"/star_index/"
    log:
        config["db_dir"]+"/log/star_index/"
    threads: 8
    shell:
        bin_dir+"STAR "
        "--runMode genomeGenerate "
        "--genomeDir {output.dir} "
        "--genomeFastaFiles {input.genome} "
        "--sjdbOverhang 100 "
        "--sjdbGTFfile {input.gtf} "
        "--runThreadN {threads} "
        "--outFileNamePrefix {output.dir} "

################################################################################
#
# read mapping
#

def comma_join(files):
    if(isinstance(files,list)):
        return ",".join(files)
    else:
        return files

rule star_1_pass:
    input:
        read1=lambda wildcards: config["fastq"][wildcards.sample_id][0],
        read2=lambda wildcards: config["fastq"][wildcards.sample_id][1],
    output:
        "star_1_pass/{sample_id}/SJ.out.tab",
        dir="star_1_pass/{sample_id}/"
    params:
        read1=lambda wildcards: comma_join(config["fastq"][wildcards.sample_id][0]),
        read2=lambda wildcards: comma_join(config["fastq"][wildcards.sample_id][1]),
        index=config["db_dir"]+"/star_index/",
        readFilesCommand=config["readFilesCommand"]
    log:
        "log/star_1_pass/{sample_id}/"
    threads: 8
    shell:
        bin_dir+"STAR "
        "--genomeDir {params.index} "
        "--readFilesIn {params.read1} {params.read2} "
        " --runThreadN {threads} "
        "--outFilterMultimapScoreRange 1 "
        "--outFilterMultimapNmax 20 "
        "--outFilterMismatchNmax 10 "
        "--alignIntronMax 500000 "
        "--alignMatesGapMax 1000000 "
        "--sjdbScore 2 "
        "--alignSJDBoverhangMin 1 "
        "--genomeLoad NoSharedMemory "
        "--readFilesCommand {params.readFilesCommand} "
        "--outFilterMatchNminOverLread 0.33 "
        "--outFilterScoreMinOverLread 0.33 "
        "--sjdbOverhang 100 "
        "--outSAMstrandField intronMotif "
        "--outSAMtype None "
        "--outSAMmode None "
        "--outFileNamePrefix {output.dir} "

rule star_2_pass:
    input:
        read1=lambda wildcards: config["fastq"][wildcards.sample_id][0],
        read2=lambda wildcards: config["fastq"][wildcards.sample_id][1],
        sj=expand("star_1_pass/{sample_id}/SJ.out.tab", sample_id=config["fastq"]),
    output:
        "star_2_pass/{sample_id}/Aligned.sortedByCoord.out.bam",
        dir="star_2_pass/{sample_id}/"
    params:
        read1=lambda wildcards: comma_join(config["fastq"][wildcards.sample_id][0]),
        read2=lambda wildcards: comma_join(config["fastq"][wildcards.sample_id][1]),
        index=config["db_dir"]+"/star_index/",
        readFilesCommand=config["readFilesCommand"],
        rg_line=lambda wildcards: "ID:1 LB:Library PL:Illumina SM:%s PU:Platform" % (wildcards.sample_id,)
    log:
        "log/star_2_pass/{sample_id}/"
    threads: 8
    shell:
        bin_dir+"STAR "
        "--genomeDir {params.index} "
        "--readFilesIn {params.read1} {params.read2} "
        "--runThreadN {threads} "
        "--outFilterMultimapScoreRange 1 "
        "--outFilterMultimapNmax 20 "
        "--outFilterMismatchNmax 10 "
        "--alignIntronMax 500000 "
        "--alignMatesGapMax 1000000 "
        "--sjdbScore 2 "
        "--alignSJDBoverhangMin 1 "
        "--genomeLoad NoSharedMemory "
        "--limitBAMsortRAM 0 "
        "--readFilesCommand {params.readFilesCommand} "
        "--outFilterMatchNminOverLread 0.33 "
        "--outFilterScoreMinOverLread 0.33 "
        "--sjdbOverhang 100 "
        "--outSAMstrandField intronMotif "
        "--outSAMattributes NH HI NM MD AS XS "
        "--outSAMunmapped Within "
        "--outSAMtype BAM SortedByCoordinate "
        "--outSAMheaderHD @HD VN:1.4 "
        "--outSAMattrRGline {params.rg_line} "
        "--outFileNamePrefix {output.dir} "
        "--sjdbFileChrStartEnd {input.sj} "
        "--limitSjdbInsertNsj 2000000 "

rule mark_duplicates:
    input:
        "star_2_pass/{sample_id}/Aligned.sortedByCoord.out.bam"
    output:
        bam="star_2_pass/{sample_id}/Aligned.sortedByCoord.out.markdup.bam",
        metrics="star_2_pass/{sample_id}/markdup_metrics.txt"
    log:
        "log/mark_duplicates/"
    shell:
        bin_dir+"picard -Xmx5G MarkDuplicates I={input} O={output.bam} M={output.metrics}"

rule index_markdup_bam:
    input:
        bam="star_2_pass/{sample_id}/Aligned.sortedByCoord.out.markdup.bam"
    output:
        bam="star_2_pass/{sample_id}/Aligned.sortedByCoord.out.markdup.bam.bai",
    log:
        "log/index_markdup_bam/"
    shell:
        bin_dir+"samtools index {input}"

################################################################################
#
# QC
#
rule rna_seqc_sample_file:
    input:
        bams=expand("star_2_pass/{sample_id}/Aligned.sortedByCoord.out.markdup.bam", sample_id=config["fastq"]),
        bais=expand("star_2_pass/{sample_id}/Aligned.sortedByCoord.out.markdup.bam.bai", sample_id=config["fastq"])
    output:
        "qc/rna_seqc/sample_file.txt"
    log:
        "log/rna_seqc_sample_file/"
    script:
        "scripts/make_rna_seqc_sample_file.py"

rule rna_seqc:
    input:
        "qc/rna_seqc/sample_file.txt"
    output:
        "qc/rna_seqc/index.html",
        dir="qc/rna_seqc/"
    params:
        java7=config["env_dir"]+"/../../pkgs/java-jdk-7.0.91-1/bin/java",
        rna_seqc_jar=config["env_dir"]+"/share/rna-seqc-1.1.8-0/RNA-SeQC_v1.1.8.jar",
        genome=config["db_dir"]+"/genome/hs37d5.fa",
        gtf=config["db_dir"]+"/gene_model/gencode.v19.annotation.hs37d5_chr.gtf",
    log:
        "log/rna_seqc/"
    shell:
        "{params.java7} -Xmx2G -jar {params.rna_seqc_jar} "
        "-s {input} "
        "-t {params.gtf} "
        "-r {params.genome} "
        "-o {output.dir} "

################################################################################
#
# expression quantification
#

rule feature_counts:
    input:
        expand("star_2_pass/{sample_id}/Aligned.sortedByCoord.out.markdup.bam", sample_id=config["fastq"])
    output:
        "expression/feature_counts/counts.txt"
    params:
        gtf=config["db_dir"]+"/gene_model/gencode.v19.annotation.hs37d5_chr.gtf",
        strandness=2 # 2 for Illumina TruSeq
    log:
        "log/feature_counts/"
    threads: 8
    shell:
        bin_dir+"featureCounts "
        "-p "
        "-T {threads} "
        "-s {params.strandness} "
        "-a {params.gtf} "
        "-o {output} "
        "{input} "

rule calc_fpkm:
    input:
        "expression/feature_counts/counts.txt"
    output:
        raw_counts="expression/raw_counts.tsv",
        fpkm="expression/fpkm.tsv",
        fpkm_uq="expression/fpkm_uq.tsv",
        dir="expression/"
    log:
        "log/calc_fpkm/"
    shell:
        bin_dir+"Rscript scripts/calc_fpkm.R {input} {output.dir}"
