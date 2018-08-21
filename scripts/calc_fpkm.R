# Compute FPKM and FPKM-UQ from raw read counts.
# Masashi Fujita, 1208/2017

library(dplyr)
library(data.table)
library(readr)

args <- commandArgs(trailingOnly=T)
count_file <- args[1]
name_file <- args[2]
outdir <- args[3]

# function for saving the expression levels in the GCT format
write_gct <- function(df, name, outfile) {
  df <- name %>%
    select(gene_name, gene_id) %>%
    inner_join(df) %>%
    rename(Name=gene_name, Description=gene_id)

    "#1.2\n" %>% cat(file=outfile)
    sprintf("%d\t%d\n", nrow(df), ncol(df)-2) %>% cat(file=outfile, append=T)
    df %>% write_tsv(outfile, append=T, col_names=T)
}

#
# load output of featureCounts
#
count <- fread(count_file, skip=1)

#
# load the mapping table between gene_id and gene_name
#
name <- fread(name_file)

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
# raw counts
#
outfile <- paste0(outdir,"/raw_counts.tsv")
data.frame(gene_id=rownames(m), m)%>% write_gct(name, outfile)

#
# FPKM
#
tot <- colSums(m)
norm_factor <- (10**9 * (1/len) %*% t(1/tot)) # normalization factor
fpkm <- m * norm_factor
outfile <- paste0(outdir,"/fpkm.tsv")
data.frame(gene_id=rownames(fpkm), fpkm)%>% write_gct(name, outfile)

#
# FPKM-UQ
#
any_read <- rowSums(m)>0
uq <- m[any_read,] %>% apply(2, function(x){quantile(x,probs=0.75)})
norm_factor <- (10**9 * (1/len) %*% t(1/uq)) # normalization factor
norm_factor <- norm_factor / nrow(m) # further normalize by the number of genes
fpkm_uq <- m * norm_factor
outfile <- paste0(outdir,"/fpkm_uq.tsv")
data.frame(gene_id=rownames(fpkm_uq), fpkm_uq)%>% write_gct(name, outfile)
