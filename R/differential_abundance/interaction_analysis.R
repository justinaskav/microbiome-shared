#' Interaction Analysis Module
#'
#' Reusable functions for per-taxon two-factor interaction analysis: forest plots, EMM
#' tables and prevalence tables (EMM plots: `plot_emm_interaction()`).
#' Generic over the two interacting variables (e.g. Diagnosis × Site).

#' Plot interaction forest for a single taxon
#'
#' @param maaslin3_result MaAsLin3 result object (interaction model)
#' @param taxon Character, taxon name
#' @param var1 Character, name of the first variable (e.g. "Diagnosis")
#' @param var2 Character, name of the second variable (e.g. "Site")
#' @param var1_levels Character vector of var1 factor levels in display order
#'   (non-reference levels only, as they appear in coefficient names)
#' @param var2_levels Character vector of var2 factor levels in display order
#'   (non-reference levels only, as they appear in coefficient names)
#' @param var2_ref_label Character, label for the reference level of var2
#'   used in the facet strip (e.g. "<reference level> (ref.)")
#' @param ref_label Character, label for the reference combination in the title
#' @param q_threshold Numeric, significance threshold (default 0.05)
#' @return ggplot object
plot_taxon_interaction_forest <- function(
  maaslin3_result,
  taxon,
  var1 = "Diagnosis",
  var2 = "Site",
  var1_levels,
  var2_levels,
  var2_ref_label,
  ref_label,
  q_threshold = 0.05
) {
  df <- prepare_maaslin3_data(
    maaslin3_result,
    model_type = "abundance",
    taxa = taxon
  )

  # Filter to coefficients involving var1 or var2
  df <- df |>
    dplyr::filter(grepl(var1, coefficient_name) | grepl(var2, coefficient_name)) |>
    dplyr::mutate(
      ci_lower = coef - 1.96 * stderr,
      ci_upper = coef + 1.96 * stderr,
      is_significant = qval < q_threshold,
      coef_type = dplyr::case_when(
        grepl(":", coefficient_name) ~ paste0("Interaction (", var1, " \u00d7 ", var2, ")"),
        grepl(paste0("^", var1), coefficient_name) ~ paste0("Main Effect (", var1, ")"),
        grepl(paste0("^", var2), coefficient_name) ~ paste0("Main Effect (", var2, ")"),
        TRUE ~ "Other"
      )
    )

  # Assign var2 group for faceting
  # Main effects of var1 (no var2 in name) go to the reference group
  df$var2_group <- NA_character_
  for (lev in var2_levels) {
    pattern <- paste0("(:", var2, lev, "|^", var2, lev, ")")
    df$var2_group[grepl(pattern, df$coefficient_name)] <- lev
  }
  # var1 main effects (no interaction) belong to the var2 reference group
  df$var2_group[grepl(paste0("^", var1), df$coefficient_name) &
                  !grepl(":", df$coefficient_name)] <- var2_ref_label

  df <- df |> dplyr::filter(!is.na(var2_group))

  # Create clean labels
  df$coef_label <- gsub(":", " \u00d7 ", gsub(var2, "", gsub(var1, "", df$coefficient_name)))

  # Assign var1 rank for ordering within each facet
  df$var1_part <- ""
  for (lev in var1_levels) {
    df$var1_part[grepl(lev, df$coefficient_name, fixed = TRUE)] <- lev
  }

  # Facet levels: reference label first, then var2_levels
  facet_levels <- c(var2_ref_label, var2_levels)
  df$var2_group <- factor(df$var2_group, levels = facet_levels)

  # Sort by facet group, then by var1 order
  var1_order <- c("", var1_levels)
  df$var1_rank <- match(df$var1_part, var1_order)
  df <- df |>
    dplyr::arrange(var2_group, var1_rank) |>
    dplyr::mutate(coef_label = factor(coef_label, levels = rev(unique(coef_label))))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = coef, y = coef_label, color = coef_type)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = ci_lower, xmax = ci_upper, shape = is_significant),
                             size = 0.5, fatten = 3) +
    ggplot2::scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                                name = paste0("q < ", q_threshold),
                                labels = c("TRUE" = "Yes", "FALSE" = "No")) +
    ggplot2::scale_color_brewer(palette = "Dark2", name = "Term Type") +
    ggplot2::facet_grid(var2_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::labs(
      x = "Beta Coefficient (95% CI)",
      y = NULL,
      title = paste0(taxon, ": Interaction (ref. ", ref_label, ")")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      legend.position = "right",
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", size = 10),
      panel.spacing = ggplot2::unit(0.3, "lines")
    )
  p
}


#' EMM table for a single taxon (two-factor interaction model)
#'
#' @param maaslin3_result MaAsLin3 result object (interaction model)
#' @param taxon Character, taxon name
#' @param var1 Character, name of the first variable (rows)
#' @param var2 Character, name of the second variable (columns)
#' @param var1_levels Character vector of var1 factor levels (NULL = the model's order)
#' @param var2_levels Character vector of var2 factor levels (NULL = the model's order)
#' @param title Optional table title. Defaults to "<taxon> EMM (var1 × var2)".
#' @param subtitle Table subtitle.
#' @return gt table object (or NULL if model not available)
table_taxon_emm <- function(
  maaslin3_result,
  taxon,
  var1 = "Diagnosis",
  var2 = "Site",
  var1_levels = NULL,
  var2_levels = NULL,
  title = NULL,
  subtitle = "Estimated marginal mean (95% CI) on log-abundance scale"
) {
  if (is.null(title)) {
    title <- paste0(taxon, " EMM (", var1, " × ", var2, ")")
  }
  model <- maaslin3_result$fit_data_abundance$fits[[taxon]]
  if (is.null(model)) return(NULL)

  emm_formula <- stats::as.formula(paste("~", var1, "*", var2))
  emm <- emmeans::emmeans(model, emm_formula)
  emm_df <- as.data.frame(summary(emm))
  if (is.null(var1_levels)) var1_levels <- levels(factor(emm_df[[var1]]))
  if (is.null(var2_levels)) var2_levels <- levels(factor(emm_df[[var2]]))

  emm_df <- emm_df |>
    dplyr::mutate(
      !!var2 := factor(.data[[var2]], levels = var2_levels),
      !!var1 := factor(.data[[var1]], levels = var1_levels),
      formatted = ifelse(is.na(emmean), "\u2014",
                         sprintf("%.2f (%.2f \u2013 %.2f)", emmean, lower.CL, upper.CL))
    ) |>
    dplyr::filter(!is.na(.data[[var2]])) |>
    dplyr::select(dplyr::all_of(c(var1, var2)), formatted) |>
    tidyr::pivot_wider(names_from = dplyr::all_of(var2), values_from = formatted)

  emm_df |>
    gt::gt() |>
    gt::tab_header(title = title, subtitle = subtitle)
}


#' Prevalence table for a single taxon, var1 (rows) \u00d7 var2 (columns)
#'
#' Percentage of samples with non-zero relative abundance per var1 \u00d7 var2 cell.
#' Argument order mirrors `table_taxon_emm()` so the two can be called side by side.
#'
#' @param tse_genus Genus-level TSE (with `relabundance` assay).
#' @param taxon Character, taxon name (matched against FeatureID).
#' @param var1 Character, row grouping variable name.
#' @param var2 Character, column grouping variable name.
#' @param var1_levels Character vector of var1 factor levels (row order; NULL = data order).
#' @param var2_levels Character vector of var2 factor levels (column order; NULL = data order).
#' @param title,subtitle gt table header text (title defaults to a generic label).
#' @return gt table object
table_taxon_prevalence <- function(
  tse_genus,
  taxon,
  var1 = "Diagnosis",
  var2 = "Site",
  var1_levels = NULL,
  var2_levels = NULL,
  title = NULL,
  subtitle = "n positive / n total (%)."
) {
  if (is.null(title)) {
    title <- paste0(taxon, " prevalence by ", var1, " \u00d7 ", var2)
  }
  cd <- SummarizedExperiment::colData(tse_genus)
  if (is.null(var1_levels)) var1_levels <- levels(factor(cd[[var1]]))
  if (is.null(var2_levels)) var2_levels <- levels(factor(cd[[var2]]))
  mia::meltSE(
    tse_genus, assay.type = "relabundance", add.col = c(var1, var2)
  ) |>
    dplyr::filter(FeatureID == taxon) |>
    dplyr::group_by(.data[[var1]], .data[[var2]]) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_pos = sum(relabundance > 0),
      prev_pct = round(n_pos / n * 100, 1),
      formatted = paste0(n_pos, "/", n, " (", prev_pct, "%)"),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      !!var2 := factor(.data[[var2]], levels = var2_levels),
      !!var1 := factor(.data[[var1]], levels = var1_levels)
    ) |>
    dplyr::select(dplyr::all_of(c(var1, var2)), formatted) |>
    tidyr::pivot_wider(names_from = dplyr::all_of(var2), values_from = formatted,
                       values_fill = "\u2014") |>
    gt::gt() |>
    gt::tab_header(title = title, subtitle = subtitle)
}
