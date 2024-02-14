import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import plotly.express as px

"""
12.5 points - You will write a program in Python/R/C (just pick a language you like – I’m listing ones here that I
prefer) that will take as input a FASTQ file and print the distribution of quality scores across all reads. You can
summarize the distribution of Q scores at each base with a statistic of your choice (e.g. mean, mode, median, quantile
distribution). If you’d like, you can also plot the distribution of Q scores as a box plot much like what’s generated
by FASTQC. You will then run your program on the provided FASTQ file, and obtain the output from it. You'll find this
FASTQ file under Files->Assignment Files->assignment1.fastq.
"""


def read_file():
    positions_dict = defaultdict(list)

    with open("assignment1.fastq", "r") as fastq:
        fastq = fastq.readlines()

    quality_scores = fastq[3::4]  # start @ line 3, grab every 4th line

    for line in quality_scores:
        for i, score in enumerate(line.strip("\n")):
            adjusted_score = ord(score) - 33
            positions_dict[i].append(adjusted_score)
    return positions_dict


def plot(scores_dict):
    scores_df = pd.DataFrame(columns=["position", "score", "method"])

    for position in sorted(scores_dict.keys()):
        median_score = np.median(scores_dict[position])
        scores_df.loc[len(scores_df)] = [position + 1, median_score, "median"]

        mean_score = np.mean(scores_dict[position])
        scores_df.loc[len(scores_df)] = [position + 1, mean_score, "mean"]

    px.line(
        scores_df,
        x="position",
        y="score",
        color="method",
        width=600,
        height=600,
    ).update_layout(
        title="Median Quality Scores",
        xaxis_title="position in read",
        yaxis_title="quality score",
    ).show()


def main():
    scores_dict = read_file()
    plot(scores_dict)


if __name__ == "__main__":
    sys.exit(main())
