library(devtools)

use_package("org.Hs.eg.db")
use_package("AnnotationDbi")
use_package("clusterProfiler")
use_package("ggplot2")


#' Gene Ontology GO enrichment analysis
#'
#'
#' @param gene_type tells the function what kind Gene ID it is, take gene symbol or ENTREZ ID
#'
#' @returns a table with over represented gene
#' @export
#'
#' @examples GOpathway(gene_type = "SYMBOL")
#' @examples GOpathway(gene_type = "ENTREZID")

GO_pathway <- function( gene_type = "SYMBOL"){
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
  head(GO@result[,c(2,3,6)])}




#' KEGG pathway enrichment analysis
#'
#' @param gene_type tells the function what kind Gene ID it is, take gene symbol or ENTREZ ID
#'
#' @returns data table that shows enrich pathway in a genne set
#' @export
#'
#' @examples KEGG_pathway(gene_type = "SYMBOL")
#' @examples KEGG_pathway(gene_type = "ENTREZID")

KEGG_pathway <- function(gene_type = "SYMBOL"){
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
  head(KEGG@result[,c(2,3,6)])}
