#!/bin/bash

# fetch data
for file in 2cells_1.fastq 2cells_2.fastq 6h_1.fastq 6h_2.fastq
do
    wget ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/$file
done

# install deps. for deseq2
sudo apt-get install -y libxml2-dev libcurl4-openssl-dev

# install python deps
pip install pandas HTSeq

# install and source STAR
cd /home/exouser/tools
wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz && tar -xzf 2.7.11b.tar.gz
cd STAR-2.7.11b/source
make STAR
echo "export PATH=$PATH:/home/exouser/tools/STAR-2.7.11b/source" >> ~/.bashrc && source ~/.bashrc

# build genome index using STAR
cd ~/workspace/BMI609/HW3
gunzip danio_genome/Danio_rerio.Zv9.66.gtf.gz
mkdir danio_genome/STAR
STAR --runMode genomeGenerate --genomeDir danio_genome/STAR --genomeFastaFiles danio_genome/Danio_rerio.Zv9.66.dna.fa --sjdbGTFfile danio_genome/Danio_rerio.Zv9.66.gtf

# align
STAR --genomeDir danio_genome/STAR --readFilesIn 2cells_1.fastq 2cells_2.fastq --runThreadN 12 --outFileNamePrefix 2cell_aligned
STAR --genomeDir danio_genome/STAR --readFilesIn 6h_1.fastq 6h_2.fastq --runThreadN 12 --outFileNamePrefix 6hour_aligned

# htseq counts table
python -m HTSeq.scripts.count 2cell_alignedAligned.out.sam danio_genome/Danio_rerio.Zv9.66.gtf > htseq_count_2_cell.txt
python -m HTSeq.scripts.count 6hour_alignedAligned.out.sam danio_genome/Danio_rerio.Zv9.66.gtf > htseq_count_6_hour.txt

# call deseq2 R file
R ./deseq2.R
