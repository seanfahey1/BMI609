#!/usr/bin/env python3
import pandas as pd
import re
import sys


regex = r'gene_id "(.*)"; transcript_id ".*"; exon_number "\d+"; gene_name "(.*)"; gene_biotype'


def extract_gene_name(row):
    match = re.search(regex, row[8])
    if match is not None:
        return match.group(2)
    else:
        return None


def extract_gene_id(row):
    match = re.search(regex, row[8])
    if match is not None:
        return match.group(1)
    else:
        return None


def main():
    with open("DESeq2_results.csv", "r") as f:
        df = pd.read_csv(f, header=0)

    with open("danio_genome/Danio_rerio.Zv9.66.gtf", "r") as gtf:
        gtf = pd.read_csv(gtf, header=None, sep="\t", low_memory=False)

    gtf["gene_name"] = gtf.apply(extract_gene_name, axis=1)
    gtf["gene_id"] = gtf.apply(extract_gene_id, axis=1)
    gtf.set_index("gene_id", inplace=True)
    gtf = gtf["gene_name"]

    df = df.join(gtf)
    with open("DESeq2_results_gene_name.csv", "w") as out:
        out.write(df.to_csv(sep=",", header=True, index=True))


if __name__ == "__main__":
    sys.exit(main())
