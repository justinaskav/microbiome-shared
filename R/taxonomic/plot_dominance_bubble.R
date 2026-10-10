#' Berger-Parker dominance per sample
#'
#' @param df Long data frame with one row per sample x taxon and a per-sample relative
#'   abundance column.
#' @param group_cols Columns identifying a sample (kept in the output).
#' @param taxon_col,rel_col Taxon and relative-abundance columns.
#' @return One row per sample: `group_cols`, berger_parker (max relative abundance) and
#'   dominant_genus (the taxon holding it).
#' @export
calc_berger_parker <- function(df, group_cols = c("SampleID", "Patient", "Site", "Diagnosis"),
                               taxon_col = "FeatureID", rel_col = "rel_abundance") {
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(berger_parker = max(.data[[rel_col]]),
                     dominant_genus = as.character(.data[[taxon_col]][which.max(.data[[rel_col]])]),
                     .groups = "drop")
}

#' Name a palette after a set of levels, extending it with hue_pal() if too short
#'
#' @param levels Levels to colour, in order.
#' @param palette Base colours.
#' @return Named colour vector.
#' @export
palette_for_levels <- function(levels, palette) {
  if (length(levels) > length(palette)) {
    palette <- c(palette, scales::hue_pal()(length(levels) - length(palette)))
  }
  stats::setNames(palette[seq_along(levels)], levels)
}

#' Subject x site dominance bubble plot
#'
#' Every sample is a point coloured by its most abundant genus. Samples at or below the
#' threshold are small faded dots; dominated samples (> threshold) are ringed circles sized
#' by Berger-Parker. Subjects are ordered within each facet by their peak dominance.
#'
#' @param bp_data Output of `calc_berger_parker()` (usually filtered to subjects with at
#'   least one dominated sample).
#' @param palette Named colour vector for `dom_genus_levels`.
#' @param dom_genus_levels Genus legend order.
#' @param site_levels Levels of `x_var` (x axis).
#' @param facet_levels Levels of `facet_var` (facet rows).
#' @param threshold Berger-Parker dominance threshold.
#' @param x_var,subject_var,facet_var Column names for the x axis, y axis and facets.
#' @return A ggplot object.
#' @export
plot_dominance_bubble <- function(bp_data, palette, dom_genus_levels,
                                  site_levels, facet_levels, threshold = 0.5,
                                  x_var = "Site", subject_var = "Patient", facet_var = "Diagnosis",
                                  size_range = c(3, 8),
                                  subthreshold_size = 2.6, subthreshold_alpha = 0.55,
                                  dominated_stroke = 0.5,
                                  grid_colour = "grey92", grid_linewidth = 0.4,
                                  caption_size = ggplot2::rel(0.8), caption_colour = "grey40") {
  bp_data <- bp_data |>
    dplyr::mutate(
      !!x_var := factor(.data[[x_var]], levels = site_levels),
      !!facet_var := factor(.data[[facet_var]], levels = facet_levels),
      dominant_genus = factor(dominant_genus, levels = dom_genus_levels),
      dominated = berger_parker > threshold
    )

  # Order subjects within each facet by their peak dominance (top = most dominated)
  subject_order <- bp_data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(facet_var, subject_var)))) |>
    dplyr::summarise(max_bp = max(berger_parker), .groups = "drop") |>
    dplyr::arrange(.data[[facet_var]], max_bp) |>
    dplyr::pull(!!rlang::sym(subject_var))
  bp_data[[subject_var]] <- factor(bp_data[[subject_var]], levels = unique(subject_order))

  ggplot2::ggplot(bp_data, ggplot2::aes(x = .data[[x_var]], y = .data[[subject_var]])) +
    # sub-threshold: small fixed-size, low-opacity dot in the genus colour. shape 16 renders
    # via `colour` (alpha applies cleanly); `fill` is mapped too so genera that never dominate
    # still populate the fill scale, otherwise their legend keys come out without swatches.
    ggplot2::geom_point(
      data = dplyr::filter(bp_data, !dominated),
      ggplot2::aes(colour = dominant_genus, fill = dominant_genus),
      shape = 16, size = subthreshold_size, alpha = subthreshold_alpha
    ) +
    # dominated: full genus fill, sized by Berger-Parker, dark ring
    ggplot2::geom_point(
      data = dplyr::filter(bp_data, dominated),
      ggplot2::aes(size = berger_parker, fill = dominant_genus),
      shape = 21, colour = "grey20", stroke = dominated_stroke
    ) +
    ggplot2::scale_size_continuous(range = size_range, limits = c(threshold, 1),
                                   breaks = c(0.6, 0.75, 0.9),
                                   name = "Berger-Parker\n(dominated)") +
    ggplot2::scale_fill_manual(
      values = palette, name = "Most abundant genus",
      limits = dom_genus_levels, drop = FALSE,
      labels = function(x) paste0("*", x, "*")
    ) +
    ggplot2::scale_colour_manual(values = palette, limits = dom_genus_levels, guide = "none") +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::facet_grid(stats::as.formula(paste(facet_var, "~ .")),
                        scales = "free_y", space = "free_y") +
    ggplot2::guides(size = ggplot2::guide_legend(override.aes = list(shape = 21, fill = "grey50")),
                    fill = ggplot2::guide_legend(override.aes = list(size = 4, shape = 21))) +
    ggplot2::labs(x = x_var, y = subject_var,
                  caption = paste0("Small faded dot = sample present, most-abundant genus shown ",
                                   "but not dominant (≤ ", threshold, "); large dot = dominated (> ",
                                   threshold, ")")) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = grid_colour, linewidth = grid_linewidth),
      legend.text = ggtext::element_markdown(),
      plot.caption = ggplot2::element_text(size = caption_size, colour = caption_colour, hjust = 0),
      axis.text.x = ggplot2::element_text(angle = 0)
    )
}
