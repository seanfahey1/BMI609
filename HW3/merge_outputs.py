#!/usr/bin/env python3
import pandas as pd


with open("DESeq2_results.csv", "r") as f:
    df = pd.read_csv(f, header=0)
df
