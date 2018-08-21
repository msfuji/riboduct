import csv
import re

gene_id_pat=re.compile("gene_id \"(\S+)\"")
gene_name_pat=re.compile("gene_name \"(\S+)\"")

with open(snakemake.output, "w") as fo:
    writer = csv.writer(fo, delimiter="\t")
    writer.writerow(["gene_id", "gene_name"])
    with open(snakemake.input) as fi:
        reader=csv.reader(fi, delimiter="\t")
        for row in reader:
            if row[0][0]=="#" or row[2] != "gene":
                continue
            m = gene_id_pat.search(row[8])
            gene_id = m.group(1)
            m = gene_name_pat.search(row[8])
            gene_name = m.group(1)
            writer.writerow([gene_id, gene_name])
