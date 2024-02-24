#!/usr/bin/env bash

# Pick an influenza genome from SRA. Analyze the FASTQ files using your own implementation of FASTQC (Assignment 1),
# and provide the output from your tool. Discuss the Q-score distributions (or any other quality metrics) that you glean
# information about - note you will not receive any points if you don't discuss your results (20 points).

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR23143762/SRR23143762
fastq-dump -I --split-files SRR23143762
python ../HW1/Question1.py --i SRR23143762_1.fastq --o SRR23143762_1_plot.html
