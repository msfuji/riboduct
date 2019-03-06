library(dplyr)
library(data.table)
library(readr)

args <- commandArgs(trailingOnly=T)
outfile <- args[1]
infiles <- args[-1]

print(infiles[1])
df <- fread(infiles[1], skip=1)
colnames(df)[7] <- colnames(df)[7] %>% dirname %>% basename
df <- df[,c(1,6,7)]

for(infile in infiles[-1]) {
  print(infile)
  t <- fread(infile, skip=1)
  colnames(t)[7] <- colnames(t)[7] %>% dirname %>% basename
  if(any(df$Geneid != t$Geneid)) {
    stop("[ERROR] Geneid mismatch.")
  }
  if(any(df$Length != t$Length)) {
    stop("[ERROR] Length mismatch.")
  }
  df <- bind_cols(df, t[,7,drop=F])
}

df %>% write_tsv(outfile)
