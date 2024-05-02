## Assignment 3
### Assignment Instructions
For this assignment, you will analyze the zebrafish RNAseq data that we had used before (with TopHat, Cufflinks + Cuffdiff) and analyze them using the HTseq + DEseq2 pipelines instead. You're welcome to do this using Galaxy, or via [Bioconductor](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). Please submit your history (if you did this using Galaxy), or your code (if you did this with Bioconductor) with your assignment, along with output files of differentially expressed genes.

Describe your results - what outliers do you find? What are the functions of these outlier loci that are differentially expressed between the 2-cell and 6-hour conditions?

### How to use
`bash ./Question1.sh` will output 2 files:

1. DESeq2_results.csv
2. DESeq2_results_gene_name.csv
