bin_dir=config["env_dir"]+"/bin/"

rule all:
    input:
        "expression/raw_counts.tsv"
        "expression/fpkm.tsv"
        "expression/fpkm_uq.tsv"

################################################################################
#
# database construction
#

rule download_genome:
    output:
        config["db_dir"]+"/genome/hs37d5.fa.gz"
    log:
        config["db_dir"]+"/log/download_genome/"
    shell:
        "curl -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz && "
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
        "curl -O ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz && "
        "mv gencode.v19.annotation.gtf.gz {output}"

rule format_gtf:
    input:
        config["db_dir"]+"/gene_model/gencode.v19.annotation.gtf.gz"
    output:
        config["db_dir"]+"/gene_model/gencode.v19.annotation.hs37d5_chr.gtf"
    log:
        config["db_dir"]+"log/format_gtf/"
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
        config["db_dir"]+"/log/setup_db/"
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
        readFilesCommand=config["readFilesCommand"]
        rg_line="@RG\tID:1\tLB:Library\tPL:Illumina\tSM:%s\tPU:Platform" % (wildcards.sample_id,)
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
        "--outSAMattrRGline {params.rg_line}"
        "--outFileNamePrefix {output.dir} "
        "--sjdbFileChrStartEnd {input.sj} "

rule mark_duplicates:
    input:
        "star_2_pass/{sample_id}/Aligned.sortedByCoord.out.bam"
    output:
        bam="star_2_pass/{sample_id}/Aligned.sortedByCoord.out.markdup.bam",
        metrics="star_2_pass/{sample_id}/markdup_metrics.txt"
    log:
        "log/mark_duplicates/"
    shell:
        bin_dir+"picard MarkDuplicates I={input} O={output.bam} M={output.metrics}"

rule index_markdup_bam:
    input:
        bam="star_2_pass/{sample_id}/Aligned.sortedByCoord.out.markdup.bam"
    output
        bam="star_2_pass/{sample_id}/Aligned.sortedByCoord.out.markdup.bam.bai",
    log:
        "log/index_markdup_bam/"
    shell:
        bindir+"samtools index {input}"

################################################################################
#
# QC
#


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
    threads: 16
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
    #    raw_counts="expression/raw_counts.tsv",
    #    fpkm="expression/fpkm.tsv",
        fpkm_uq="expression/fpkm_uq.tsv",
        dir="expression/"
    log:
        "log/calc_fpkm/"
    shell:
        bin_dir+"Rscript scripts/calc_fpkm.R {input} {output.dir}"
