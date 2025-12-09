library(devtools)
use_package("edgeR")

gene_table <- function(count_file, sample_file ){
  count_table <- read.table(count_file, header = TRUE, as.is = TRUE, row.names = 1, sep = "\t")
  sample_table <- read.table(sample_file, header = TRUE, as.is = TRUE, row.names = 1, sep = "\t")

  sample_table$individual <- factor(sample_table$individual)
  sample_table$sex <- factor(sample_table$sex)
  sample_table$diseas <- factor(sample_table$diseas)

  table_maching <- colnames(count_table) == rownames(sample_table)

  if (table_maching == FALSE)
    {stop("Count_table and sample table are not matching !")}

  return(list(head(count_table), sample_table, a))}

gene_table("E-MTAB-2523.counts.txt", "E-MTAB-2523_sample table.txt")

low_gene_filtering <- function(cutoff) {
filt <- rowMeans(log2(cpm(count_table)+1))
sum(filt <= 1)
count_table_filt <- count_table[filt > cutoff,]
dim(count_table)
dim(count_table_filt)
}

low_gene_filtering(1)
