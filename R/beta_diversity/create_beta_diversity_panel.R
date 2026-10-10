#' Create Beta Diversity Panel
#'
#' Wraps `analyze_beta_diversity()` + `plot_beta_diversity()` for standard
#' distance metrics (Jaccard, Bray-Curtis, unweighted/weighted UniFrac).
#'
#' With `strata_var` set, the PERMANOVA only validly tests terms that vary within a stratum;
#' between-subject terms need `beta_permanova_split()`.
#'
#' @param tse A TreeSummarizedExperiment object with pre-computed ordinations.
#' @param group_var Grouping variable for PERMANOVA (default: "Group").
#' @param formula PERMANOVA formula.
#' @param strata_var Stratification variable for PERMANOVA (default: "Patient").
#' @param color_by Variable for point color (default: "Group").
#' @param shape_by Variable for point shape (default: "Type").
#' @param permutations Number of PERMANOVA / betadisper permutations (default: 999).
#' @param metric_filter Subset of c("jaccard", "bray", "unifrac_unweighted",
#'   "unifrac_weighted") to include (default: all).
#' @param analyze Run the PERMANOVA, betadisper and pairwise tests (default: TRUE). With FALSE
#'   only the axis variance-explained labels are computed, for panels whose inference lives
#'   elsewhere (e.g. `beta_permanova_split()`).
#' @param post_process Optional function to apply to each individual plot
#'   (e.g., `function(p) p + theme_bw()`).
#' @param seed Optional RNG seed forwarded to `analyze_beta_diversity()` so each metric's
#'   permutation p-values are reproducible independently of upstream RNG state. `NULL`
#'   (default) preserves the previous behaviour of drawing from the ambient stream.
#'
#' @return A list with `analysis` (list of beta diversity analysis results)
#'   and `plots` (list of ggplot objects). Caller does layout/save.
#' @export
create_beta_diversity_panel <- function(
  tse,
  group_var = "Group",
  formula,
  strata_var = "Patient",
  color_by = "Group",
  shape_by = "Type",
  permutations = 999,
  post_process = NULL,
  metric_filter = NULL,
  seed = NULL,
  analyze = TRUE
) {
  metrics <- list(
    jaccard = list(dimred = "MDS_jaccard", distance = "jaccard", title = "Jaccard"),
    bray = list(dimred = "MDS_bray", distance = "bray", title = "Bray-Curtis"),
    unifrac_unweighted = list(dimred = "MDS_unifrac", distance = "unifrac", title = "Unweighted UniFrac"),
    unifrac_weighted = list(dimred = "MDS_unifrac_weighted", distance = "unifrac_weighted", title = "Weighted UniFrac")
  )

  if (!is.null(metric_filter)) {
    metrics <- metrics[intersect(metric_filter, names(metrics))]
  }

  beta_analysis <- lapply(metrics, function(m) {
    if (!analyze) return(.variance_explained(tse, m$distance))
    analyze_beta_diversity(
      tse = tse, dimred_name = m$dimred, distance_label = m$distance,
      group_var = group_var, formula = formula, strata_var = strata_var,
      permutations = permutations, seed = seed
    )
  })

  beta_plots <- lapply(names(metrics), function(mn) {
    m <- metrics[[mn]]
    p <- plot_beta_diversity(
      tse = tse, dimred_name = m$dimred, color_by = color_by,
      shape_by = shape_by, color_label = color_by, shape_label = shape_by,
      analysis_result = beta_analysis[[mn]], inset_type = "none",
      color_palette = "Set2", title = m$title
    )

    if (!is.null(post_process)) {
      p <- post_process(p)
    }

    p
  })

  list(analysis = beta_analysis, plots = beta_plots)
}
