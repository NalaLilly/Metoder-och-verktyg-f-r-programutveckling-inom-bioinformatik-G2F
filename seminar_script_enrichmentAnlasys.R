library(devtools)

BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")
BiocManager::install("clusterProfiler")
BiocManager::install("ggplot2")

use_package("org.Hs.eg.db")
use_package("AnnotationDbi")
use_package("clusterProfiler")
use_package("ggplot2")

count_table <- read.table("E-MTAB-2523.counts.txt", header = TRUE, as.is = TRUE, row.names = 1, sep = "\t")
sample_table <- read.table("E-MTAB-2523_sample table.txt", header = TRUE, as.is = TRUE, row.names = 1, sep = "\t")

GO_pathway <- function(organism = org.Hs.eg.db, gene_type = "SYMBOL"){
  statical_analysis<- run_dge_edger(count_table, sample_table,
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
  head(GO@result[,c(2,3,5)])}


GO_pathway(organism = org.Hs.eg.db, gene_type = "SYMBOL")


KEGG_pathway <- function(organism = org.Hs.eg.db){
  a <- run_dge_edger(count_table, sample_table,
                group_col     = "disease",
                case_label    = "carcinoma",
                control_label = "normal",
                fdr_cutoff    = 0.05,
                lfc_cutoff    = 1)

  signficant_genes <- a[[2]]
  genes_names <- rownames(signficant_genes)
  genes_Entrez <- org.Hs.eg.db::mapIds(org.Hs.eg.db,
                                            keys=genes_names,
                                            keytype="SYMBOL",
                                            column="ENTREZID")
  genes_cl10_Entrez<-genes_cl10_Entrez[!duplicated(genes_cl10_Entrez) &
                                         !duplicated(genes_cl10_Entrez,fromLast=TRUE)]
  paste0("Number of genes in after removing elements with duplicated names and or NAs = ",
         length(genes_cl10_Entrez))
  }
