#' Create RDA Analysis with Plot and ANOVA Table
#'
#' Wraps `plot_rda_summary()` and `run_dbrda_analysis()` (Type I table) into a single call.
#'
#' @param tse A TreeSummarizedExperiment object with relabundance assay.
#' @param rda_formula_str RDA formula string (e.g., "~ Group + Type").
#' @param color_by Variable for point color (default: "Group").
#' @param shape_by Variable for point shape (default: "Type").
#' @param strata_var Stratification variable for ANOVA (default: "Patient").
#' @param arrow_scaling Scaling factor for biplot arrows (default: 0.7).
#' @param permutations Number of permutations for the figure's term labels and the ANOVA
#'   (default: 999). Match this to the `run_dbrda_analysis()` call for the same block, or the
#'   figure and the table report the same fit at different permutation counts.
#' @param post_process Optional function to apply to the RDA plot
#'   (e.g., `function(p) p + theme_bw(base_size = 14) + background_grid()`).
#'
#' @return A list with `plot` (ggplot object) and `anova_df` (`run_dbrda_analysis()$type1`,
#'   NULL on failure). Caller does print/save/kable.
create_rda_analysis <- function(
  tse,
  rda_formula_str,
  color_by = "Group",
  shape_by = NULL,
  strata_var = NULL,
  arrow_scaling = 0.7,
  permutations = 999,
  post_process = NULL
) {
  p_rda <- tryCatch({
    p <- plot_rda_summary(
      tse = tse, rda_formula_str = rda_formula_str,
      color_by = color_by, shape_by = shape_by,
      color_label = color_by, shape_label = shape_by,
      arrow_scaling = arrow_scaling, strata_var = strata_var,
      permutations = permutations
    )
    if (!is.null(post_process)) post_process(p) else p
  }, error = function(e) {
    ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::labs(title = paste("RDA failed:", e$message))
  })

  list(plot = p_rda,
       anova_df = run_dbrda_analysis(tse, rda_formula_str, strata_var, permutations)$type1)
}
