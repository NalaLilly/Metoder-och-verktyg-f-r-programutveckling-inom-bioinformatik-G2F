library(devtools)

BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")
BiocManager::install("clusterProfiler")
BiocManager::install("ggplot2")

use_package("org.Hs.eg.db")
use_package("AnnotationDbi")
use_package("clusterProfiler")
use_package("ggplot2")




GO_pathway <- function(organism = org.Hs.eg.db, gene_type = "SYMBOL"){
  statical_analysis <- run_dge_edger(count_table, sample_table,
                     group_col     = "disease",
                     case_label    = "carcinoma",
                     control_label = "normal",
                     fdr_cutoff    = 0.05,
                     lfc_cutoff    = 1)

  significant_genes <- statical_analysis[[2]]
  genes_names <- rownames(significant_genes)
  library(org.Hs.eg.db)

  GO <- clusterProfiler::enrichGO(gene = genes_names,
                            OrgDb = organism,
                            keyType = gene_type,
                            ont = "BP",
                            pvalueCutoff = 0.05)
  head(GO@result[,c(2,3,6)])}





KEGG_pathway <- function(organism = hsa, organism_data = org.Hs.eg.db, gene_type = "SYMBOL"){
  statical_analysis <- run_dge_edger(count_table, sample_table,
                group_col     = "disease",
                case_label    = "carcinoma",
                control_label = "normal",
                fdr_cutoff    = 0.05,
                lfc_cutoff    = 1)

  signficant_genes <- statical_analysis[[2]]
  genes_names <- rownames(signficant_genes)
  library("org.Hs.eg.db")
  genes_Entrez <- AnnotationDbi::mapIds(organism_data,
                                            keys=genes_names,
                                            keytype="SYMBOL",
                                            column="ENTREZID")

  genes_Entrez<-genes_Entrez[!duplicated(genes_Entrez) &
                                         !duplicated(genes_Entrez, fromLast=TRUE)]
  paste0("Number of genes in after removing elements with duplicated names and or NAs = ",
         length(genes_Entrez))
  KEGG <- clusterProfiler::enrichKEGG(gene = genes_Entrez,
                                      organism = organism,
                                      keyType = "ncbi-geneid",
                                      pvalueCutoff = 0.05)
  head(KEGG@result[,c(2,3,6)])}

