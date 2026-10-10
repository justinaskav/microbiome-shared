#' Plot Metadata Correlation Heatmap
#'
#' Calculates and plots correlations between metadata variables using ComplexHeatmap.
#' Handles mixed data types via mia::getCrossAssociation.
#'
#' @param tse A TreeSummarizedExperiment object.
#' @param variables Character vector of variable names from colData(tse).
#' @param pretty_labels Optional named character vector to rename variables for plotting.
#' @param method Correlation method (passed to getCrossAssociation, default "spearman").
#' @param p_adjust_method P-value adjustment method across the unique variable pairs (default "BH").
#' @param p_threshold Significance threshold for marking cells (default 0.05).
#' @param heatmap_title Title for the heatmap (default: "Metadata Cross-Association").
#' @param base_fontsize Base font size for heatmap elements (default: 9).
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap.
#'
#' @return A Heatmap object. Needs to be drawn using ComplexHeatmap::draw().
#'
#' @examples
#' \dontrun{
#' ht <- plot_metadata_correlation_heatmap(
#'   tse = tse_object,
#'   variables = c("Age", "BMI", "DiseaseScore"),
#'   method = "spearman",
#'   p_threshold = 0.05
#' )
#' ComplexHeatmap::draw(ht)
#' }
#'
#' @export
plot_metadata_correlation_heatmap <- function(tse,
                                              variables,
                                              pretty_labels = NULL,
                                              method = "spearman",
                                              p_adjust_method = "BH",
                                              p_threshold = 0.05,
                                              heatmap_title = NULL,
                                              base_fontsize = 9,
                                              ...) {

  # --- Calculate Cross Association ---
  # Raw p-values from mia, adjusted here over the unique variable pairs only: the matrix holds
  # every pair twice plus the diagonal, which would change the size of the BH family.
  res <- mia::getCrossAssociation(
    tse, tse,
    col.var1 = variables, col.var2 = variables,
    method = method,
    p.adjust.method = "none",
    test.signif = TRUE,
    mode = "matrix"
  )

  if (is.null(res) || is.null(res$cor) || is.null(res$pval)) {
    stop("Failed to calculate cross-associations. Check input variables and data.")
  }
  res$pval <- res$pval[rownames(res$cor), colnames(res$cor)]
  ut <- upper.tri(res$pval)
  res$p_adj <- matrix(NA_real_, nrow(res$pval), ncol(res$pval), dimnames = dimnames(res$pval))
  res$p_adj[ut] <- stats::p.adjust(res$pval[ut], method = p_adjust_method)
  res$p_adj[lower.tri(res$p_adj)] <- t(res$p_adj)[lower.tri(res$p_adj)]

  final_heatmap_title <- if (is.null(heatmap_title)) "Metadata Cross-Association" else heatmap_title

  # Apply pretty labels if provided
  if (!is.null(pretty_labels)) {
    current_vars <- rownames(res$cor)
    missing_labels <- setdiff(current_vars, names(pretty_labels))
    if (length(missing_labels) > 0) {
      warning("Missing pretty labels for: ", paste(missing_labels, collapse = ", "))
      names(missing_labels) <- missing_labels
      pretty_labels <- c(pretty_labels, missing_labels)
    }
    rownames(res$cor) <- pretty_labels[rownames(res$cor)]
    colnames(res$cor) <- pretty_labels[colnames(res$cor)]
    rownames(res$p_adj) <- pretty_labels[rownames(res$p_adj)]
    colnames(res$p_adj) <- pretty_labels[colnames(res$p_adj)]
  }

  # --- Create Heatmap ---
  cor_mat <- res$cor
  p_adj_mat <- res$p_adj
  max_abs_cor <- max(abs(cor_mat), na.rm = TRUE)
  if (is.na(max_abs_cor) || max_abs_cor == 0) max_abs_cor <- 1

  col_fun <- circlize::colorRamp2(
    c(-max_abs_cor, 0, max_abs_cor),
    c("steelblue", "white", "firebrick")
  )

  # Layer function for significance markers
  layer_fun_sig <- function(j, i, x, y, width, height, fill) {
    p_adj_val <- ComplexHeatmap::pindex(p_adj_mat, i, j)
    is_significant <- !is.na(p_adj_val) & p_adj_val < p_threshold
    grid::grid.text("*", x[is_significant], y[is_significant],
      gp = grid::gpar(fontsize = base_fontsize + 3, col = "black"))
  }

  ht <- ComplexHeatmap::Heatmap(
    cor_mat,
    name = "Correlation",
    col = col_fun,
    heatmap_legend_param = list(
      title = "Correlation",
      title_gp = grid::gpar(fontsize = base_fontsize, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = base_fontsize),
      legend_height = grid::unit(5, "cm"),
      at = round(seq(-max_abs_cor, max_abs_cor, length.out = 5), 2),
      direction = "vertical"
    ),
    column_names_rot = 45,
    row_names_gp = grid::gpar(fontsize = base_fontsize),
    column_names_gp = grid::gpar(fontsize = base_fontsize),
    column_title = final_heatmap_title,
    column_title_gp = grid::gpar(fontsize = base_fontsize + 2, fontface = "bold"),
    rect_gp = grid::gpar(col = "grey", lty = 1, lwd = 0.5),
    layer_fun = layer_fun_sig,
    ...
  )

  return(ht)
}
