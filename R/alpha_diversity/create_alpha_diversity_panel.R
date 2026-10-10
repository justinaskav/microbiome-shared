#' Display label for an alpha-diversity index column.
#'
#' Index columns are named for the code (`Faith_PD`), not for the reader. Anything not listed
#' falls through unchanged.
#'
#' @param idx Character vector of index column names.
#' @return Character vector of display labels.
#' @export
alpha_index_label <- function(idx) {
  labels <- c(Faith_PD = "Faith’s PD", Shannon = "Shannon", Observed = "Observed ASVs")
  out <- unname(labels[idx])
  ifelse(is.na(out), idx, out)
}

#' Create Alpha Diversity Panel
#'
#' Wraps `plot_alpha_diversity_summary_lmm()` for multiple diversity indices
#' into a combined patchwork layout.
#'
#' @param data A data frame containing alpha diversity values and metadata.
#' @param x_var Character string, the grouping variable for the x-axis.
#' @param available_indices Character vector of diversity index names to plot.
#' @param lmm_formula_str LMM formula string with "INDEX" as placeholder for
#'   the response variable (e.g., "INDEX ~ Group + Sex + (1|Patient)").
#' @param ref_group Reference group for post-hoc comparisons.
#' @param x_label Label for the x-axis (default: x_var).
#' @param color_var Variable for color aesthetic (default: x_var).
#' @param shape_var Variable for shape aesthetic (default: "Batch").
#' @param palette RColorBrewer palette name (default: "Set2").
#' @param color_values Named vector of custom colors (overrides palette if provided).
#' @param post_process Optional function to apply to each individual plot
#'   (e.g., `function(p) p + panel_border() + background_grid()`).
#' @param connect_subjects Logical, connect each subject's points (see
#'   `plot_alpha_diversity_summary_lmm()`; default: FALSE).
#' @param facet_var Optional character string; when set, each index plot is faceted
#'   by this variable using `facet_wrap`. Post-hoc significance brackets are
#'   suppressed when faceting: their y positions come from the whole dataset and would be
#'   repeated identically in every facet.
#' @param show_lmm_label Logical, whether each index plot carries the model-term
#'   annotation label at top-right (default: TRUE).
#' @param style Named list of cosmetic parameters from `alpha_plot_style()`, forwarded to
#'   `plot_alpha_diversity_summary_lmm()`. Defaults reproduce the historical hardcoded values.
#'
#' @return A patchwork object combining all alpha diversity plots.
#' @export
create_alpha_diversity_panel <- function(
  data,
  x_var,
  available_indices,
  lmm_formula_str,
  ref_group,
  x_label = NULL,
  color_var = NULL,
  shape_var = "Batch",
  palette = "Set2",
  color_values = NULL,
  post_process = NULL,
  connect_subjects = FALSE,
  facet_var = NULL,
  show_lmm_label = TRUE,
  style = alpha_plot_style()
) {
  if (is.null(x_label)) x_label <- x_var
  if (is.null(color_var)) color_var <- x_var

  alpha_plots <- lapply(available_indices, function(idx) {
    p <- plot_alpha_diversity_summary_lmm(
      data = data,
      x_var = x_var,
      y_var = idx,
      color_var = color_var,
      shape_var = shape_var,
      x_label = x_label,
      y_label = alpha_index_label(idx),
      color_label = color_var,
      shape_label = shape_var,
      palette = palette,
      lmm_formula_str = gsub("INDEX", idx, lmm_formula_str),
      ref_group = ref_group,
      post_hoc_method = "BH",
      posthoc = is.null(facet_var),
      show_posthoc_label = FALSE,
      show_lmm_label = show_lmm_label,
      connect_subjects = connect_subjects,
      style = style
    )

    # Add faceting after the plot is built
    if (!is.null(facet_var)) {
      p <- p + ggplot2::facet_wrap(reformulate(facet_var))
    }

    # Apply custom colors if provided
    if (!is.null(color_values)) {
      p <- p +
        ggplot2::scale_color_manual(values = color_values, na.value = "grey50") +
        ggplot2::scale_fill_manual(values = color_values, na.value = "grey50")
    }

    # Apply post-processing (theming, axis adjustments, etc.)
    if (!is.null(post_process)) {
      p <- post_process(p)
    }

    p
  })

  ncol_val <- if (!is.null(facet_var)) 1 else 2
  patchwork::wrap_plots(alpha_plots, ncol = ncol_val, guides = "collect")
}
