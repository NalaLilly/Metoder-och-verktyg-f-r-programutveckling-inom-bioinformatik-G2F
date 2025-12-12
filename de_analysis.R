#' Differential expression analysis with edgeR
#'
#' Performs differential expression analysis between two groups
#' (e.g. carcinoma vs normal) using edgeR on a count matrix.
#'
#' @param counts Integer matrix of raw counts (genes x samples).
#' @param sample_table Data frame with sample metadata. Row names
#'   must match the column names of \code{counts}.
#' @param group_col Name of the column in \code{sample_table}
#'   that defines the groups (e.g. "disease").
#' @param case_label Level in \code{group_col} you want as "case"
#'   (e.g. "carcinoma").
#' @param control_label Level in \code{group_col} you want as "control"
#'   (e.g. "normal").
#' @param fdr_cutoff Numeric, FDR threshold for significance.
#' @param lfc_cutoff Numeric, absolute log2 fold change cutoff.
#'
#' @return A list with:
#'   \item{results}{Data frame with all genes, including logFC, FDR and DEG flag.}
#'   \item{deg}{Data frame with only significant DEGs.}
#'
#' @export
run_dge_edger <- function(counts,
                          sample_table,
                          group_col     = "disease",
                          case_label    = "carcinoma",
                          control_label = "normal",
                          fdr_cutoff    = 0.05,
                          lfc_cutoff    = 1) {

  if (!all(colnames(counts) %in% rownames(sample_table))) {
    stop("Column names of 'counts' must match row names of 'sample_table'.")
  }

  sample_table <- sample_table[colnames(counts), , drop = FALSE]

  group <- factor(sample_table[[group_col]])

  if (!all(c(case_label, control_label) %in% levels(group))) {
    stop("case_label and control_label must be valid group levels.")
  }

  group <- relevel(group, ref = control_label)

  y <- edgeR::DGEList(counts = counts, group = group)
  y <- edgeR::calcNormFactors(y)

  design <- model.matrix(~ group)
  y <- edgeR::estimateDisp(y, design)

  fit <- edgeR::glmQLFit(y, design)
  qlf <- edgeR::glmQLFTest(fit, coef = 2)

  tt <- edgeR::topTags(qlf, n = Inf)$table
  tt$gene <- rownames(tt)

  tt$DEG <- with(tt, FDR < fdr_cutoff & abs(logFC) >= lfc_cutoff)

  deg_only <- tt[tt$DEG, , drop = FALSE]

  return(list(
    results = tt,
    deg     = deg_only
  ))
}
