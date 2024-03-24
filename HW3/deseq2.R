if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

library("DESeq2")
directory <- "/home/exouser/workspace/BMI609/HW3/"

sampleFiles <- grep(
    "htseq_count_.*txt",
    list.files(directory),
    value=TRUE
    )

sampleCondition <- sub("htseq_count_(.*).txt","\\1",sampleFiles)

sampleTable <- data.frame(
    sampleName = sampleFiles,
    fileName = sampleFiles,
    condition = sampleCondition
    )

sampleTable$condition <- factor(sampleTable$condition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(
    sampleTable = sampleTable,
    directory = directory,
    design= ~ condition
    )

# note: I'm hitting an error here...
dds <- DESeq(ddsHTSeq)
res <- results(dds, contrast=c("condition","2cell","6hour"))
