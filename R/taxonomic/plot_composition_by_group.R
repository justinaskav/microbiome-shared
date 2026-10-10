#' Mean composition per facet x group, as unweighted means of per-sample relative abundance
#'
#' Each sample's taxon abundance is divided by that sample's total from `totals` (not by
#' the sum of `df`), so features missing from `df` (e.g. unclassified reads dropped at genus
#' level) stay in the denominator instead of inflating the named taxa. Taxa absent from a
#' sample count as 0 before averaging.
#'
#' @param df Long data frame (e.g. from `mia::meltSE()`), one row per sample x feature.
#' @param taxon_col Column holding the (already lumped) taxon label.
#' @param totals Data frame of per-sample totals: `sample_col` and `total`.
#' @param facet_var,x_var Grouping columns: one facet per `facet_var` level, one bar per
#'   `x_var` level within it.
#' @param facet_levels,x_levels Factor levels for the two grouping columns.
#' @param x_abbrev Optional named vector abbreviating `x_var` levels in the bar labels.
#' @param sample_meta Data frame with `sample_col`, `facet_var`, `x_var` used to count
#'   samples per bar (default: the distinct rows of `df`).
#' @param sample_col,count_col Sample-ID and count columns.
#' @return Tibble with `facet_var`, `x_var`, taxon, rel, n and x (bar label: abbreviated
#'   group over its sample count, ordered by group then n).
#' @export
calc_group_mean_composition <- function(df, taxon_col, totals, facet_var, x_var,
                                        facet_levels, x_levels, x_abbrev = NULL,
                                        sample_meta = NULL,
                                        sample_col = "SampleID", count_col = "counts") {
  if (is.null(x_abbrev)) x_abbrev <- stats::setNames(x_levels, x_levels)
  if (is.null(sample_meta)) {
    sample_meta <- dplyr::distinct(df, dplyr::across(dplyr::all_of(c(sample_col, facet_var, x_var))))
  }
  group_n <- sample_meta |> dplyr::count(dplyr::across(dplyr::all_of(c(facet_var, x_var))), name = "n")
  keys <- c(sample_col, facet_var, x_var)

  df |>
    dplyr::left_join(totals, by = sample_col) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys)), taxon = .data[[taxon_col]]) |>
    dplyr::summarise(rel = sum(.data[[count_col]]) / dplyr::first(total), .groups = "drop") |>
    tidyr::complete(tidyr::nesting(!!!rlang::syms(keys)), taxon, fill = list(rel = 0)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(facet_var, x_var))), taxon) |>
    dplyr::summarise(rel = mean(rel), .groups = "drop") |>
    dplyr::left_join(group_n, by = c(facet_var, x_var)) |>
    dplyr::mutate(
      !!facet_var := factor(.data[[facet_var]], levels = facet_levels),
      !!x_var := factor(.data[[x_var]], levels = x_levels)
    ) |>
    dplyr::arrange(.data[[x_var]], n) |>
    dplyr::mutate(x = forcats::fct_inorder(paste0(x_abbrev[as.character(.data[[x_var]])], "\n", n)))
}

#' Shade taxa within their parent group (e.g. genera within phyla)
#'
#' Gives each taxon a shade of its parent's hue, in input order (pass taxa sorted by
#' parent, then abundance, to get dark = more abundant).
#'
#' @param taxa Character vector of taxa.
#' @param parents Parent of each taxon (factor; its level order sets the palette order).
#' @param shades Named list: parent -> character vector of shades, dark to light.
#' @return Named colour vector, names = taxa, ordered by parent then input order.
#' @export
shade_by_parent <- function(taxa, parents, shades) {
  if (anyNA(parents)) {
    stop("taxa without a parent: ", paste(taxa[is.na(parents)], collapse = ", "))
  }
  parents <- droplevels(as.factor(parents))
  missing <- setdiff(levels(parents), names(shades))
  if (length(missing) > 0) {
    stop("no shades defined for: ", paste(missing, collapse = ", "))
  }
  by_parent <- split(taxa, parents)
  cols <- lapply(names(by_parent), function(ph) {
    g <- by_parent[[ph]]
    if (length(g) > length(shades[[ph]])) {
      stop(length(g), " taxa in ", ph, " but only ", length(shades[[ph]]), " shades")
    }
    stats::setNames(shades[[ph]][seq_along(g)], g)
  })
  unlist(cols)
}

#' Stacked-bar composition panel, one facet per group
#'
#' @param df Output of `calc_group_mean_composition()` with `taxon` as a factor (stack order).
#' @param cols Named colour vector for taxa.
#' @param levels Taxa in legend order.
#' @param legend_title Legend title.
#' @param legend_rows Legend rows (filled column-wise).
#' @param facet_var Facet column.
#' @param plain_pattern Regex for labels shown upright (lumped categories); all other taxa are
#'   italicised and have "_" replaced by spaces.
#' @param wrap_width Facet-label wrap width.
#' @param base_size Base font size for `theme_classic()`.
#' @return A ggplot object.
#' @export
plot_group_composition <- function(df, cols, levels, legend_title, legend_rows, facet_var,
                                   plain_pattern = "^Other", wrap_width = 12, base_size = 9) {
  taxon_labels <- function(x) {
    ifelse(grepl(plain_pattern, x), x, paste0("*", gsub("_", " ", x), "*"))
  }
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = rel, fill = taxon)) +
    ggplot2::geom_col(width = 0.85, color = "white", linewidth = 0.25) +
    ggplot2::facet_grid(stats::as.formula(paste("~", facet_var)),
                        scales = "free_x", space = "free_x",
                        labeller = ggplot2::label_wrap_gen(wrap_width)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                expand = ggplot2::expansion(mult = c(0, 0.01)),
                                breaks = seq(0, 1, 0.25)) +
    ggplot2::scale_fill_manual(values = cols, breaks = levels, labels = taxon_labels,
                               name = legend_title) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = legend_rows, byrow = FALSE)) +
    ggplot2::labs(x = NULL, y = "Mean relative abundance") +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
      panel.spacing.x = ggplot2::unit(0.8, "lines"),
      axis.text.x = ggplot2::element_text(color = "grey20", lineheight = 0.9),
      axis.text.y = ggplot2::element_text(color = "grey20"),
      axis.line = ggplot2::element_line(color = "grey40", linewidth = 0.3),
      axis.ticks = ggplot2::element_line(color = "grey40", linewidth = 0.3),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 6)),
      legend.text = ggtext::element_markdown(),
      legend.key.size = ggplot2::unit(0.75, "lines"),
      legend.position = "bottom",
      legend.title.position = "left",
      legend.title = ggplot2::element_text(face = "bold", margin = ggplot2::margin(r = 8)),
      plot.tag = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2))
    )
}
