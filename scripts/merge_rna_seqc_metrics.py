import pandas as pd

#
# load and concat single-sample metrics files
#
m_list=[]
for metrics in snakemake.input:
    m=pd.read_csv(metrics, sep="\t", index_col=0)
    m_list.append(m)
df=pd.concat(m_list, axis=1)

df.to_csv(snakemake.output[0], sep="\t")
