library(devtools)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

use_package("edgeR")

gene_table <- function(count_file, sample_file ){
  count_table <- read.table(count_file, header = TRUE, as.is = TRUE, row.names = 1, sep = "\t")
  sample_table <- read.table(sample_file, header = TRUE, as.is = TRUE, row.names = 1, sep = "\t")

  sample_table$individual <- factor(sample_table$individual)
  sample_table$sex <- factor(sample_table$sex)
  sample_table$diseas <- factor(sample_table$diseas)
  a <- colnames(count_table) == rownames(sample_table)
  gene_data <- list(count_table,sample_table))
  return(gene_data)}



low_gene_filtering <- function(cutoff) {
  gene_data <- gene_table("E-MTAB-2523.counts.txt", "E-MTAB-2523_sample table.txt")
  count_table <- gene_data[[1]]
  sample_table <- gene_data[[2]]
  filt <- rowMeans(log2(edgeR::cpm(count_table)+1))
  sum(filt <= 1)
  count_table_filt <- count_table[filt > cutoff,]

  msg <- paste("filtering works! data befoe filtering ",
               paste(dim(count_table), collapse = "x"),
               "after filtering",
                paste(dim(count_table_filt), collapse = "x"))

  message(msg)
  return(list(count_table_filt, sample_table))}

low_gene_filtering(1)

