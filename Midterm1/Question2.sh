#!/usr/bin/env bash

# Assemble the above genome using the same k-mer size using
#    (a) velvet (https://github.com/dzerbino/velvetLinks to an external site.)
#    (b) megahit (https://github.com/voutcn/megahitLinks to an external site.)
#    (c) unicycler (https://github.com/rrwick/UnicyclerLinks to an external site.)
#    (d) masurca (https://github.com/alekseyzimin/masurcaLinks to an external site.)

# Then compare the final assemblies with respect to the H3N2 reference genome using QUAST. Discuss the assemblies with
# respect to the 3C criteria (completeness, contiguity, and correctness) - once again, note that you will not receive
# any points if you don't discuss your results (50 points).

# Please also submit your QUAST output (as .html file) with your answers.

# VELVET
velveth velvet_25 25 -shortPaired -separate -fastq SRR23143762_1.fastq SRR23143762_2.fastq
velvetg velvet_25/

# MEGAHIT
megahit --k-list 25 --no-mercy -1 SRR23143762_1.fastq -2 SRR23143762_2.fastq -o megahit_25/

# MASURCA
masurca -g masurca.cfg
# edit mascura.cfg to set .fastq files correctly and GRAPH_KMER_SIZE to 25
masurca masurca.cfg
./assemble.sh

# compare with QUAST
mkdir assemblies
cp velvet_25/contigs.fa assemblies/velvet.fa
cp megahit_25/final.contigs.fa assemblies/megahit.fa
cp masurca_25/CA/primary.genome.scf.fasta assemblies/masurca.fasta
python3 /home/exouser/tools/quast/quast.py --min-contig 200 -f -b -r h3n2.fna.gz assemblies/velvet.fa assemblies/megahit.fa assemblies/masurca.fasta
