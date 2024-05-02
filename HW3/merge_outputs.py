#!/usr/bin/env python3
import re
import sys

import pandas as pd

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
    # open files as data frames
    with open("DESeq2_results.csv", "r") as f:
        df = pd.read_csv(f, header=0)
    df.index.name = "GeneID"

    with open("danio_genome/Danio_rerio.Zv9.66.gtf", "r") as gtf:
        gtf = pd.read_csv(gtf, header=None, sep="\t", low_memory=False)

    # extract desired data into clean columns
    gtf["gene_name"] = gtf.apply(extract_gene_name, axis=1)
    gtf["gene_id"] = gtf.apply(extract_gene_id, axis=1)

    # hold on to only needed columns
    gtf.set_index("gene_id", inplace=True)
    gtf = gtf["gene_name"]

    # convert to dict
    gtf_dict = gtf.to_dict()

    # map dict into new col
    df["gene_name"] = df.index.map(gtf_dict)

    # sort and store
    df = df.sort_values(by="padj", ascending=True)

    with open("DESeq2_results_gene_name.csv", "w") as out:
        out.write(df.to_csv(sep=",", header=True, index=True))


if __name__ == "__main__":
    sys.exit(main())
