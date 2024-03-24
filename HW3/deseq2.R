if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library("DESeq2")
directory <- "/home/exouser/workspace/BMI609/HW3/"

sample_files <- grep(
  "htseq_count_.*txt",
  list.files(directory),
  value = TRUE
)

sample_condition <- sub("htseq_count_(.*).txt", "\\1", sample_files)

sample_table <- data.frame(
  sampleName = sample_files,
  fileName = sample_files,
  condition = sample_condition
)

sample_table$condition <- factor(sample_table$condition)

dds_ht_seq <- DESeqDataSetFromHTSeqCount(
  sampleTable = sample_table,
  directory = directory,
  design = ~ 1
)

dds <- DESeq(dds_ht_seq)
res <- results(dds, alpha = 0.05)
res_clean <- res[!is.na(res$padj), ]
res_clean <- subset(res_clean, padj < 0.05)

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
