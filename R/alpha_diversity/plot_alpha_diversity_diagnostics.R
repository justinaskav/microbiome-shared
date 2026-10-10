#' Plot Alpha Diversity Model Diagnostics
#'
#' Visual counterpart to the Levene / Shapiro-Wilk tables returned by
#' `calc_alpha_diversity_summary()`. Builds a residual-vs-fitted plot and a normal
#' QQ plot from the `diagnostics$residuals` data frame.
#'
#' These plots carry more weight than the accompanying p-values: at the sample sizes
#' used here Levene's and Shapiro-Wilk tests detect deviations too small to affect the
#' F-tests, so the tests are best read as confirmation of what the plots show.
#'
#' @param residuals_df The `diagnostics$residuals` data frame from
#'   `calc_alpha_diversity_summary()`, with columns `.fitted`, `.resid`,
#'   `.std_resid` and `group`.
#' @param title Plot title, typically the index and dataset (e.g. "Shannon (gut)").
#' @param group_label Legend label for the grouping variable (default: "Group").
#' @param palette RColorBrewer palette name (default: "Set2").
#'
#' @return A patchwork object: residual-vs-fitted, residual spread by group, and a
#'   normal QQ plot.
plot_alpha_diversity_diagnostics <- function(
  residuals_df,
  title,
  group_label = "Group",
  palette = "Set2"
) {
  if (is.null(residuals_df) || nrow(residuals_df) == 0) return(NULL)

  # Title goes on each subplot rather than via plot_annotation(): these rows get
  # nested inside an outer wrap_plots(), which drops an inner patchwork's annotation.
  p_resid <- ggplot2::ggplot(residuals_df, ggplot2::aes(x = .fitted, y = .resid)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_point(ggplot2::aes(color = group), alpha = 0.7, size = 1.6) +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                         color = "grey20", linewidth = 0.6) +
    ggplot2::scale_color_brewer(palette = palette, na.value = "grey50") +
    ggplot2::labs(title = title, subtitle = "Residuals vs fitted",
                  x = "Fitted value", y = "Residual", color = group_label)

  # Direct visual counterpart to Levene's test, which compares the mean absolute
  # deviation from the group median across these groups.
  p_spread <- ggplot2::ggplot(residuals_df, ggplot2::aes(x = group, y = .resid)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_boxplot(ggplot2::aes(fill = group), alpha = 0.6,
                          outlier.shape = NA, show.legend = FALSE) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.4, size = 1.2) +
    ggplot2::scale_fill_brewer(palette = palette, na.value = "grey50") +
    ggplot2::scale_x_discrete(labels = scales::label_wrap(10)) +
    ggplot2::labs(title = title, subtitle = "Residual spread by group",
                  x = group_label, y = "Residual")

  p_qq <- ggplot2::ggplot(residuals_df, ggplot2::aes(sample = .std_resid)) +
    ggplot2::stat_qq(alpha = 0.7, size = 1.6) +
    ggplot2::stat_qq_line(color = "grey20", linewidth = 0.6) +
    ggplot2::labs(title = title, subtitle = "Normal Q-Q",
                  x = "Theoretical quantiles", y = "Standardised residual")

  patchwork::wrap_plots(list(p_resid, p_spread, p_qq), ncol = 3)
}

#' Build a Diagnostics Panel Across Indices
#'
#' Stacks `plot_alpha_diversity_diagnostics()` output for several diversity indices
#' drawn from one `calc_alpha_diversity_summary()` result set.
#'
#' @param summaries Named list of `calc_alpha_diversity_summary()` results, keyed by
#'   index (e.g. `alpha_results$gut` subset to one grouping variable).
#' @param dataset_label Label for the data subset, used in panel titles.
#' @param group_label Legend label for the grouping variable.
#' @param palette RColorBrewer palette name (default: "Set2").
#'
#' @return A patchwork object with one diagnostics row per index, or NULL if no
#'   index yielded diagnostics.
create_alpha_diagnostics_panel <- function(
  summaries,
  dataset_label,
  group_label = "Group",
  palette = "Set2"
) {
  rows <- purrr::compact(lapply(names(summaries), function(idx) {
    diag <- summaries[[idx]]$diagnostics
    if (is.null(diag)) return(NULL)
    plot_alpha_diversity_diagnostics(
      residuals_df = diag$residuals,
      title = paste0(idx, " (", dataset_label, ")"),
      group_label = group_label,
      palette = palette
    )
  }))

  if (length(rows) == 0) return(NULL)
  patchwork::wrap_plots(rows, ncol = 1)
}
