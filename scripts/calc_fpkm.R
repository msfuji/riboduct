# Compute FPKM and FPKM-UQ from raw read counts.
# Masashi Fujita, 1208/2017

library(dplyr)
library(data.table)
library(readr)

# args <- commandArgs(trailingOnly=T)
# count_file <- args[1]
# outdir <- args[2]

count_file <- "../raw_counts.txt"
outdir <- "./"

#
# load output of featureCounts
#
count <- fread(count_file, skip=1)

#
# make a matrix of raw counts
#
end <- ncol(count)
m <- count[,7:end] %>% as.matrix
colnames(m) <- count[,7:end] %>% colnames %>% dirname %>% basename
rownames(m) <- count$Geneid

# gene length
len <- count$Length

#
# FPKM-UQ
#
outfile <- paste0(outdir,"/raw_counts.tsv")
data.frame(gene_id=rownames(m), m)%>% write_tsv(outfile)

#
# FPKM-UQ
#
tot <- colSums(m)
norm_factor <- (10**9 * (1/len) %*% t(1/tot)) # normalization factor
fpkm <- m * norm_factor
outfile <- paste0(outdir,"/fpkm.tsv")
data.frame(gene_id=rownames(fpkm), fpkm)%>% write_tsv(outfile)

#
# FPKM-UQ
#
any_read <- rowSums(m)>0
uq <- m[any_read,] %>% apply(2, function(x){quantile(x,probs=0.75)})
norm_factor <- (10**9 * (1/len) %*% t(1/uq)) # normalization factor
norm_factor <- norm_factor / nrow(m) # further normalize by the number of genes
fpkm_uq <- m * norm_factor
outfile <- paste0(outdir,"/fpkm_uq.tsv")
data.frame(gene_id=rownames(fpkm_uq), fpkm_uq)%>% write_tsv(outfile)

