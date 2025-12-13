library(devtools)
use_package("edgeR")


#' import count and sample file
#'
#' @param count_file is the count file for RNA-seq experiment
#' @param sample_file is the sample file where you define which group each sample belongs to. assumes the file has header and is separated with tab
#'
#' @returns count and sample file  as data tables assumes the file has header and is separated with tab
#' @export
#'
#' @examples gene_table("E-MTAB-2523.counts.txt", "E-MTAB-2523_sample table.txt")

gene_table <- function(count_file, sample_file ){
  count_table <- read.table(count_file, header = TRUE, as.is = TRUE, row.names = 1, sep = "\t")
  sample_table <- read.table(sample_file, header = TRUE, as.is = TRUE, row.names = 1, sep = "\t")
  sample_table$individual <- factor(sample_table$individual)
  sample_table$sex <- factor(sample_table$sex)
  sample_table$diseas <- factor(sample_table$diseas)
  a <- colnames(count_table) == rownames(sample_table)
  gene_data <- list(count_table,sample_table)
  return(gene_data)}



#' filter out the low count of the count table
#'
#' @param cutoff tells the function to filter out all below the cutoff
#' @param count_table  is the count file for RNA-seq experiment. assumes the file has header and is separated with tab
#' @param sample_table is the sample file where you define which group each sample belongs to.assumes the file has header and is separated with tab
#'
#' @returns filtered counts
#' @export
#'
#' @examples gene_table(cutoff = 1,"E-MTAB-2523.counts.txt", "E-MTAB-2523_sample table.txt")
low_gene_filtering <- function(cutoff = 1, count_table, sample_table) {
  gene_data <- gene_table(count_table, sample_table)
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




