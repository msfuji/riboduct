import csv
import os.path

bams=snakemake.input[0]

with open(snakemake.output[0], "w") as f:
    writer=csv.writer(f, delimiter="\t")
    writer.writerow(["Sample ID", "Bam File", "Notes"])
    for bam in bams:
        sample_id=os.path.basename(os.path.dirname(bam))
        writer.writerow([sample_id, bam, ""])
