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

  return(count_table)}



low_gene_filtering <- function(cutoff) {
  countdata <- gene_table("E-MTAB-2523.counts.txt", "E-MTAB-2523_sample table.txt")
  filt <- rowMeans(log2(edgeR::cpm(countdata)+1))
  sum(filt <= 1)
  count_table_filt <- countdata[filt > cutoff,]
  return(list(dim(countdata), dim(count_table_filt)))}

low_gene_filtering(1)
#comentar
