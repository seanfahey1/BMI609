#!/usr/bin/bash

# For this exam, you will use the corresponding paired end reads to run through the entire MOTHUR/QIIME2/DADA2/KRAKEN2 pipeline (just pick one). IMPORTANT: No Galaxy use! You will do all analyses in the command line interface of the associated tool, collate all code you used and submit the code with your exam.

# Tutorials, installation instructions for each tool are here (you only need to run your data through one of those pipelines):
# DADA2 - https://benjjneb.github.io/dada2
# MOTHUR - https://mothur.org/wiki/miseq_sop
# QIIME2 - https://docs.qiime2.org/2024.2/interfaces/q2cli
# KRAKEN2 - https://github.com/DerrickWood/kraken2/wiki/Manual

# Thereon, each one of you will use the FASTQ files I've assigned you (see below), and analyze it. Then, submit all plots, result files, and your code used in a Zipped folder. You'll find all FASTQ files under Files->Assignment Files->Midterm 2
# zr13074_4V3V4_R1.fastq.gz and zr13074_4V3V4_R1.fastq.gz

sudo apt-get install mothur
wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v132.tgz

gunzip -k reads/zr13074_4V3V4_R1.fastq.gz
gunzip -k reads/zr13074_4V3V4_R2.fastq.gz

mothur > make.file(inputdir=reads, type=fastq, prefix=zr13074)
