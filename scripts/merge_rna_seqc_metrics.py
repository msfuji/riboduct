import pandas as pd

#
# load and concat single-sample metrics files
#
m_list=[]
for metrics in snakemake.input:
    m=pd.read_csv(metrics, sep="\t")
    m_list.append(m)
df=pd.concat(m_list)

df=df.ix[:,[
"Sample", "Note", "Total Purity Filtered Reads Sequenced",
"Alternative Aligments", "Failed Vendor QC Check", "Read Length",
"Estimated Library Size", "Mapped", "Mapping Rate", "Mapped Unique",
"Mapped Unique Rate of Total", "Unique Rate of Mapped",
"Duplication Rate of Mapped", "Base Mismatch Rate", "rRNA", "rRNA rate",
"Mapped Pairs", "Unpaired Reads", "End 1 Mapping Rate", "End 2 Mapping Rate",
"End 1 Mismatch Rate", "End 2 Mismatch Rate", "Fragment Length Mean",
"Fragment Length StdDev", "Chimeric Pairs", "Intragenic Rate", "Exonic Rate",
"Intronic Rate", "Intergenic Rate", "Split Reads",
"Expression Profiling Efficiency", "Transcripts Detected", "Genes Detected",
"End 1 Sense", "End 1 Antisense", "End 2 Sense", "End 2 Antisense",
"End 1 \% Sense", "End 2 \% Sense"
]]

df.to_csv(snakemake.output[0], sep="\t")
