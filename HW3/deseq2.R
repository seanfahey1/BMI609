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
    design= ~ 1
    )

dds <- DESeq(ddsHTSeq)
res <- results(dds)
res_clean = res[!is.na(res$padj),]
res_clean = subset(res_clean, padj < 0.05)

write.table(
    res_clean,
    file = "DESeq2_results.csv",
    append = FALSE,
    quote = TRUE,
    sep = ",",
    eol = "\n",
    na = "NA",
    row.names = TRUE,
    col.names = TRUE,
    )
