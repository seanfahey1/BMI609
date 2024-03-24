library("DESeq2")
directory <- "/path/to/your/files/"

sampleFiles <- grep(
    "htseq_count_",
    list.files(directory),
    value=TRUE
    )

sampleCondition <- sub("(.*htseq_count_).*","\\1",sampleFiles)

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

dds <- DESeq(ddsHTSeq)
res <- results(dds, contrast=c("condition","2cell","6hour"))
