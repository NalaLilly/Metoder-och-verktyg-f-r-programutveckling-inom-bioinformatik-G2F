library(devtools)

BiocManager::install("org.Hs.eg.db")
use_package("org.Hs.eg.db")

ORA <- function(){
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
