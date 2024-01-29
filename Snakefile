bin_dir=config["env_dir"]+"/bin/"

rule all:
    input:
        "expression/raw_counts.gct",
        "expression/fpkm.gct",
        # "expression/fpkm_uq.gct",
        "qc/rna_seqc/metrics.txt"

################################################################################
#
# database indexing
#

rule index:
    input:
        config["db_dir"]+"/genome/genome.fa.fai",
        config["db_dir"]+"/genome/genome.dict",
        config["db_dir"]+"/star_index/SAindex",
        config["db_dir"]+"/gene_model/gene_name.txt",
        config["db_dir"]+"/gene_model/annotation.exon.gtf",
        config["db_dir"]+"/gene_model/collapsed.gtf"

rule link_genome:
    input:
        config["genome_fa_gz"]
    output:
        config["db_dir"]+"/genome/genome.fa.gz"
    log:
        config["db_dir"]+"/log/link_genome/"
    params:
        memory="5.3G"
    run:
        import os
        os.symlink(input[0], output[0])

rule decompress_genome:
    input:
        config["db_dir"]+"/genome/genome.fa.gz"
    output:
        config["db_dir"]+"/genome/genome.fa"
    log:
        config["db_dir"]+"/log/decompress_genome/"
    params:
        memory="5.3G"
    shell:
        "gunzip -c {input} > {output}"

rule faidx_genome:
    input:
        config["db_dir"]+"/genome/genome.fa"
    output:
        config["db_dir"]+"/genome/genome.fa.fai"
    log:
        config["db_dir"]+"/log/faidx_genome/"
    params:
        memory="5.3G"
    shell:
        bin_dir+"samtools faidx {input}"

rule dict_genome:
    input:
        config["db_dir"]+"/genome/genome.fa"
    output:
        config["db_dir"]+"/genome/genome.dict"
    log:
        config["db_dir"]+"/log/dict_genome/"
    params:
        memory="20G"
    shell:
        bin_dir+"picard CreateSequenceDictionary "
        "R={input} "
        "O={output}"

rule link_gtf:
    input:
        config["annotation_gtf_gz"]
    output:
        config["db_dir"]+"/gene_model/annotation.gtf.gz"
    log:
        config["db_dir"]+"/log/link_gtf/"
    params:
        memory="5.3G"
    run:
        import os
        os.symlink(input[0], output[0])

rule decompress_gtf:
    input:
        config["db_dir"]+"/gene_model/annotation.gtf.gz"
    output:
        config["db_dir"]+"/gene_model/annotation.gtf"
    log:
        config["db_dir"]+"/log/decompress_gtf/"
    params:
        memory="5.3G"
    shell:
        "gunzip -c {input} > {output}"

rule exon_gtf:
    input:
        config["db_dir"]+"/gene_model/annotation.gtf"
    output:
        config["db_dir"]+"/gene_model/annotation.exon.gtf"
    log:
        config["db_dir"]+"/log/exon_gtf/"
    params:
        memory="5.3G"
    shell:
        "awk '$3==\"exon\"' {input} > {output}"

rule collapse_gtf:
    input:
        config["db_dir"]+"/gene_model/annotation.gtf"
    output:
        config["db_dir"]+"/gene_model/collapsed.gtf"
    log:
        config["db_dir"]+"/log/collapse_gtf/"
    params:
        memory="5.3G"
    shell:
        bin_dir+"python3 "
        "scripts/collapse_annotation.py {input} {output}"

rule star_index:
    input:
        genome=config["db_dir"]+"/genome/genome.fa",
        gtf=config["db_dir"]+"/gene_model/annotation.gtf"
    output:
        config["db_dir"]+"/star_index/SAindex"
    log:
        config["db_dir"]+"/log/star_index/"
    threads: 4
    params:
        memory="5.3G",
        outdir=config["db_dir"]+"/star_index/"
    shell:
        bin_dir+"STAR "
        "--runMode genomeGenerate "
        "--genomeDir {params.outdir} "
        "--genomeFastaFiles {input.genome} "
        "--sjdbOverhang 100 "
        "--sjdbGTFfile {input.gtf} "
        "--runThreadN {threads} "
        "--outFileNamePrefix {params.outdir} "

rule extract_gene_name:
    input:
        config["db_dir"]+"/gene_model/annotation.gtf"
    output:
        config["db_dir"]+"/gene_model/gene_name.txt"
    log:
        config["db_dir"]+"/log/extract_gene_name/"
    params:
        memory="5.3G"
    script:
        "scripts/extract_gene_name_from_gtf.py"


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
        "star_1_pass/{sample_id}/SJ.out.tab"
    params:
        read1=lambda wildcards: comma_join(config["fastq"][wildcards.sample_id][0]),
        read2=lambda wildcards: comma_join(config["fastq"][wildcards.sample_id][1]),
        index=config["db_dir"]+"/star_index/",
        readFilesCommand=config["readFilesCommand"],
        memory="5.3G",
        outdir="star_1_pass/{sample_id}/"
    log:
        "log/star_1_pass/{sample_id}/"
    threads: 4
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
        "--outFileNamePrefix {params.outdir} "

rule star_2_pass:
    input:
        read1=lambda wildcards: config["fastq"][wildcards.sample_id][0],
        read2=lambda wildcards: config["fastq"][wildcards.sample_id][1],
        sj=expand("star_1_pass/{sample_id}/SJ.out.tab", sample_id=config["fastq"]),
    output:
        "star_2_pass/{sample_id}/Aligned.sortedByCoord.out.bam"
    params:
        read1=lambda wildcards: comma_join(config["fastq"][wildcards.sample_id][0]),
        read2=lambda wildcards: comma_join(config["fastq"][wildcards.sample_id][1]),
        index=config["db_dir"]+"/star_index/",
        readFilesCommand=config["readFilesCommand"],
        rg_line=lambda wildcards: "ID:1 LB:Library PL:Illumina SM:%s PU:Platform" % (wildcards.sample_id,),
        memory="10.6G",
        outdir="star_2_pass/{sample_id}/"
    log:
        "log/star_2_pass/{sample_id}/"
    threads: 4
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
        "--outFileNamePrefix {params.outdir} "
        "--sjdbFileChrStartEnd {input.sj} "
        "--limitSjdbInsertNsj 2000000 "

rule mark_duplicates:
    input:
        "star_2_pass/{sample_id}/Aligned.sortedByCoord.out.bam"
    output:
        bam="star_2_pass/{sample_id}/Aligned.sortedByCoord.out.markdup.bam",
        metrics="star_2_pass/{sample_id}/markdup_metrics.txt"
    log:
        "log/mark_duplicates/{sample_id}/"
    params:
        memory="16G"
    shell:
        bin_dir+"picard -Xmx8G MarkDuplicates I={input} O={output.bam} M={output.metrics}"

rule index_markdup_bam:
    input:
        bam="star_2_pass/{sample_id}/Aligned.sortedByCoord.out.markdup.bam"
    output:
        bam="star_2_pass/{sample_id}/Aligned.sortedByCoord.out.markdup.bam.bai",
    log:
        "log/index_markdup_bam/{sample_id}/"
    params:
        memory="5.3G"
    shell:
        bin_dir+"samtools index {input}"

################################################################################
#
# QC
#
rule rna_seqc:
    input:
        bam="star_2_pass/{sample_id}/Aligned.sortedByCoord.out.markdup.bam",
        bai="star_2_pass/{sample_id}/Aligned.sortedByCoord.out.markdup.bam.bai"
    output:
        "qc/rna_seqc/{sample_id}/{sample_id}.metrics.tsv"
    params:
        gtf=config["db_dir"]+"/gene_model/collapsed.gtf",
        memory="5.3G",
        outdir="qc/rna_seqc/{sample_id}/"
    log:
        "log/rna_seqc/{sample_id}/"
    shell:
        bin_dir+"rnaseqc "
        "-s \"{wildcards.sample_id}\" "
        "{params.gtf} {input.bam} {params.outdir}"


rule merge_rna_seqc:
    input:
        expand("qc/rna_seqc/{sample_id}/{sample_id}.metrics.tsv", sample_id=config["fastq"])
    output:
        "qc/rna_seqc/metrics.txt"
    log:
        "log/merge_rna_seqc/"
    params:
        memory="5.3G"
    script:
        "scripts/merge_rna_seqc_metrics.py"

################################################################################
#
# expression quantification
#
rule feature_counts:
    input:
        "star_2_pass/{sample_id}/Aligned.sortedByCoord.out.markdup.bam"
    output:
        "expression/feature_counts/{sample_id}/counts.txt"
    params:
        gtf=config["db_dir"]+"/gene_model/annotation.gtf",
        strandness=2,  # 2 for Illumina TruSeq
        memory="5.3G"
    log:
        "log/feature_counts/{sample_id}/"
    threads: 1
    shell:
        bin_dir+"featureCounts "
        "-p "
        "-T {threads} "
        "-s {params.strandness} "
        "-a {params.gtf} "
        "-o {output} "
        "{input}"

rule merge_feature_counts:
    input:
        expand("expression/feature_counts/{sample_id}/counts.txt", sample_id=config["fastq"])
    output:
        "expression/feature_counts/counts.txt"
    params:
        memory="5.3G"
    log:
        "log/merge_feature_counts/"
    shell:
        bin_dir+"Rscript scripts/merge_feature_counts.R {output} {input}"

rule calc_fpkm:
    input:
        fcount="expression/feature_counts/counts.txt",
        name=config["db_dir"]+"/gene_model/gene_name.txt"
    output:
        raw_counts="expression/raw_counts.gct",
        fpkm="expression/fpkm.gct"
        # fpkm_uq="expression/fpkm_uq.gct"
    log:
        "log/calc_fpkm/"
    params:
        memory="5.3G",
        outdir="expression/"
    shell:
        bin_dir+"Rscript scripts/calc_fpkm.R {input.fcount} {input.name} {params.outdir}"
