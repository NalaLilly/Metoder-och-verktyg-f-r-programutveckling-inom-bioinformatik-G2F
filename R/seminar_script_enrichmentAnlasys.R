library(devtools)

use_package("org.Hs.eg.db")
use_package("AnnotationDbi")
use_package("clusterProfiler")
use_package("ggplot2")



#' Gene Ontology GO enrichment analysis
#'
#'
#' @param gene_type tells the function what kind Gene ID it is, take gene symbol or ENTREZ ID
#' @param count_file is the count file for RNA-seq experiment, assumes the file has header and is separated with tab
#' @param sample_file is the sample file where you define which group each sample belongs to. assumes the file has header and is separated with tab
#'
#' @returns Export the results of the Gene Ontology enrichment analysis analysis to a CVS file
#' @export
#'
#' @examples GOpathway(gene_type = "SYMBOL")
#' @examples GOpathway(gene_type = "ENTREZID")

GO_pathway <- function(count_table, sample_table, gene_type = "SYMBOL"){

  gene_filterd <-low_gene_filtering(cutoff = 1, count_table, sample_table)
  count_table <- gene_filterd[[1]]
  sample_table <- gene_filterd[[2]]

  # basic argument check so the user knows what to pass
  if (missing(count_table) || missing(sample_table)) {
    stop("Please provide 'count_table' and 'sample_table' (e.g. from gene_table()).")
  }
  statistical_analysis <- run_dge_edger(count_table, sample_table,
                     group_col     = "disease",
                     case_label    = "carcinoma",
                     control_label = "normal",
                     fdr_cutoff    = 0.05,
                     lfc_cutoff    = 1)

  significant_genes <- statistical_analysis[[2]]
  genes_names <- rownames(significant_genes)

  GO <- clusterProfiler::enrichGO(gene = genes_names,
                            OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                            keyType = gene_type,
                            ont = "BP",
                            pvalueCutoff = 0.05)

  write.table(GO, file = "Result_GO_analysisi.csv", sep = "\t",row.names = FALSE)

}



#' KEGG pathway enrichment analysis
#'
#' @param gene_type tells the function what kind Gene ID it is, take gene symbol or ENTREZ ID
#'
#' @returns data table that shows enrich pathway in a genne set
#' @export
#'
#' @examples KEGG_pathway(gene_type = "SYMBOL")
#' @examples KEGG_pathway(gene_type = "ENTREZID")

KEGG_pathway <- function(count_table, sample_table, gene_type = "SYMBOL"){

  gene_filterd <-low_gene_filtering(cutoff = 1, count_table, sample_table)
  count_table <- gene_filterd[[1]]
  sample_table <- gene_filterd[[2]]

  # basic argument check so the user knows what to pass
  if (missing(count_table) || missing(sample_table)) {
    stop("Please provide 'count_table' and 'sample_table' (e.g. from gene_table()).")
  }
  statical_analysis <- run_dge_edger(count_table, sample_table,
                group_col     = "disease",
                case_label    = "carcinoma",
                control_label = "normal",
                fdr_cutoff    = 0.05,
                lfc_cutoff    = 1)

  signficant_genes <- statical_analysis[[2]]
  genes_names <- rownames(signficant_genes)
  genes_Entrez <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                            keys=genes_names,
                                            keytype="SYMBOL",
                                            column="ENTREZID")

  genes_Entrez<-genes_Entrez[!duplicated(genes_Entrez) &
                                         !duplicated(genes_Entrez, fromLast=TRUE)]
  paste0("Number of genes in after removing elements with duplicated names and or NAs = ",
         length(genes_Entrez))
  KEGG <- clusterProfiler::enrichKEGG(gene = genes_Entrez,
                                      organism = "hsa",
                                      keyType = "ncbi-geneid",
                                      pvalueCutoff = 0.05)

  write.table(KEGG, file = "Result_KEGG_analysisi.csv", sep = "\t",row.names = FALSE)

}
