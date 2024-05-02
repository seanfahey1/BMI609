import re
import sys
from collections import defaultdict
from pathlib import Path
import numpy as np

import pandas as pd
from tqdm import tqdm
from Bio import Entrez
import re
from statsmodels.stats.multitest import multipletests
from scipy.stats import rankdata, ttest_ind


def collect_initial_df(file_path):
    with open(Path(file_path), "r") as file:
        init_df = pd.read_csv(file, header=0, sep="\t", index_col=0)

    # small cleanup to deal with whitespace in output file
    init_df.index = [re.sub("\s", "", x) for x in init_df.index]
    init_df.index.name = "gene"

    # quick fix for named tuples step later (names can't start with number)
    init_df.columns = [re.sub("1buff", "buffer", x) for x in init_df.columns]

    return init_df


def calc_direction(row, target):
    buffer_col_names = [x for x in row.index if x.startswith("buffer")]
    target_col_names = [x for x in row.index if x.startswith(target)]

    buffer_average = np.mean([getattr(row, x) for x in buffer_col_names])
    target_average = np.mean([getattr(row, x) for x in target_col_names])

    target_diff = np.mean(target_average) - np.mean(buffer_average)
    target_dir = "+" if target_diff > 0 else "-" if target_diff < 0 else "0"

    return target_dir


def unadjusted_pval(group1, group2):
    _, p = ttest_ind(group1, group2, nan_policy="raise")

    # return 1 if the p val can't be calculated (eg: comparing 0,0,0 to 0,0,0)
    if np.isnan(p):
        return 1
    else:
        return p


def correction(p_vals, fdr=False):
    # keeping while deciding which to use
    if fdr:
        ranked_p_values = rankdata(p_vals, method="average")
        corr_p = p_vals * len(p_vals) / ranked_p_values
        corr_p[corr_p > 1] = 1

        return corr_p

    else:
        corr_p = multipletests(p_vals, method="bonferroni", is_sorted=False)

        return corr_p[1]


def get_gene_info(gene):
    Entrez.email = "sfahey1745@sdsu.edu"
    with open("Final/entrez.key", "r") as key:
        Entrez.api_key = key.read()

    gene_clean = gene.strip("\n")
    query = f"{gene_clean}[gene] AND homo sapiens[All Fields]"

    with Entrez.esearch(db="gene", term=query, idtype="acc", usehistory="y") as handle:
        esearch_handler = Entrez.read(handle)
    try:
        with Entrez.efetch(
            db="gene",
            ret_mode="text",
            retmax=1,
            webenv=esearch_handler["WebEnv"],
            query_key=esearch_handler["QueryKey"],
            idtype="acc",
        ) as fetch_handle:
            data = fetch_handle.read()
            data_clean = re.sub("\n", " ", data)

            name_dirty = re.search("prot {\s+name (.*?)desc", data_clean)
            summary_dirty = re.search("summary (.*)location {", data_clean)

            if name_dirty is not None and summary_dirty is not None:
                name = re.sub("\s+", " ", name_dirty.group(1))
                name = re.sub(r"[\{\}\"\']", "", name).strip()
                name = re.sub("[ ,:;]$", "", name)
                name = re.sub(",", "", name)
                name = name.strip(" ")

                summary = re.sub("\s+", " ", summary_dirty.group(1))
                summary = re.sub(r"[\{\}\"\']", "", summary).strip()
                summary = re.sub("[ ,:;]$", "", summary)
                summary = re.sub(",", "", summary)
                summary = summary.strip(" ")

                # just kick out so we can end the session before the next query
                return name, summary

            else:
                return "", ""

    except:
        return "", ""


def calculate_diff_expression(df):
    all_columns = df.columns
    new_columns = ["buffer", "wnt3a", "wnt5a"]
    normalized_data = defaultdict(lambda: defaultdict(list))
    for row in df.itertuples():
        for group in new_columns:
            control = row.noTX
            columns = [x for x in all_columns if x.startswith(group)]
            for column in columns:
                normalized_data[row.Index][group].append(getattr(row, column) - control)

    return normalized_data


def main():
    init_df = collect_initial_df("Final/ALL-NPC_RPKM_VALUES.txt")

    # determine (average) direction of change in expression
    wnt3a_dir = init_df.apply(lambda row: calc_direction(row, "wnt3a"), axis=1)
    wnt5a_dir = init_df.apply(lambda row: calc_direction(row, "wnt5a"), axis=1)

    wnt3a_dir.name = "wnt3a_direction"
    wnt5a_dir.name = "wnt5a_direction"

    # p-value analysis
    p_val_data = {"gene": [], "wnt3a": [], "wnt5a": []}
    normalized_data = calculate_diff_expression(init_df)

    for gene in tqdm(normalized_data.keys()):
        control_group = normalized_data[gene]["buffer"]
        wnt3a_p = unadjusted_pval(normalized_data[gene]["wnt3a"], control_group)
        wnt5a_p = unadjusted_pval(normalized_data[gene]["wnt5a"], control_group)

        p_val_data["gene"].append(gene)
        p_val_data["wnt3a"].append(wnt3a_p)
        p_val_data["wnt5a"].append(wnt5a_p)

    p_val_df = pd.DataFrame.from_dict(p_val_data).set_index("gene")
    adj_p_val_df = p_val_df.apply(correction, axis=0)

    # add a temp sum column for sorting
    adj_p_val_df["sum"] = adj_p_val_df.apply(sum, axis=1)

    # join directions into df
    adj_p_val_df = adj_p_val_df.join(wnt3a_dir, how="outer", validate="one_to_one")
    adj_p_val_df = adj_p_val_df.join(wnt5a_dir, how="outer", validate="one_to_one")

    # add in unadjusted p-values
    adj_p_val_df = adj_p_val_df.join(
        p_val_df, how="outer", validate="one_to_one", rsuffix="_unadjusted"
    )

    # sort and drop temp column
    adj_p_val_df.sort_values(by="sum", inplace=True, ascending=True)
    adj_p_val_df.drop("sum", axis=1, inplace=True)

    # drop genes where the values never changed
    no_change_genes = adj_p_val_df[
        (adj_p_val_df["wnt3a_direction"] == "0")
        & (adj_p_val_df["wnt5a_direction"] == "0")
    ].index
    adj_p_val_df.drop(no_change_genes, axis=0, inplace=True)

    # get top genes info from NCBI
    sample_num = 50
    names = []
    descriptions = []
    genes = list(adj_p_val_df.index)
    for gene in tqdm(genes[:sample_num]):
        name, description = get_gene_info(gene)
        names.append(name)
        descriptions.append(description)

    for _ in genes[sample_num:]:
        names.append("")
        descriptions.append("")

    adj_p_val_df["full_names"] = names
    adj_p_val_df["descriptions"] = descriptions

    # write csv file
    with open("Final/ALL-adj-p-values.csv", "w") as file:
        file.write(adj_p_val_df.to_csv(header=True, index=True))


if __name__ == "__main__":
    sys.exit(main())
