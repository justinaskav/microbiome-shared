#' MaAsLin3 Forest Plot
#'
#' Creates a forest plot showing beta coefficients with 95% confidence intervals
#' (coef +/- 1.96 * SE) for MaAsLin3 differential abundance results.
#'
#' @param maaslin3_input MaAsLin3 input: result object, data frame, or TSV path.
#'   See `prepare_maaslin3_data()` for details.
#' @param coefficient_names Character vector of coefficient names to display
#'   (e.g., `"TypeStent"` or `c("GroupBenign Stricture", "GroupIntra Biliary Malignancy")`).
#' @param model_type One of `"abundance"`, `"prevalence"`, or `"joint"` (default: `"abundance"`).
#' @param taxa Optional character vector of specific taxa to display.
#' @param top_n Integer, number of taxa to show (by q-value). Ignored if `taxa`
#'   is provided (default: 20).
#' @param q_threshold Numeric threshold for q-value significance (default: 0.05).
#' @param sort_by How to sort taxa: `"coef"` (absolute coefficient), `"qval"`,
#'   or `"taxon"` (alphabetical) (default: `"coef"`).
#' @param color_by_coefficient Logical, color points by coefficient name when
#'   multiple coefficients are shown (default: TRUE).
#' @param show_q_legend Logical, whether to show the q-value significance shape
#'   legend (default: TRUE). When FALSE, shapes are still used but the legend is hidden.
#' @param coefficient_labels Named character vector to rename coefficient names
#'   in the legend (e.g., `c("GroupBenign Stricture" = "Benign Stricture")`).
#' @param color_values Named character vector of colors for coefficients. For
#'   single-coefficient plots, a single color value. For multi-coefficient plots,
#'   maps to coefficient names (after label renaming if `coefficient_labels` is used).
#' @param reference_line Where to draw the vertical dashed reference line(s) on the
#'   coefficient axis. One of:
#'   * `"zero"` — x = 0 (default, kept for backward compatibility). For the abundance
#'     model this is NOT the null: MaAsLin3's median comparison
#'     (`median_comparison_abundance = TRUE`) tests each coefficient against the median
#'     coefficient for that metadatum, not against 0.
#'   * `"null"` — the value MaAsLin3 actually tested against, read from the
#'     `null_hypothesis` column (the median-comparison reference for abundance; 0 for
#'     prevalence). This is the statistically faithful reference. Because each
#'     coefficient has its own null, a multi-coefficient plot draws one line per
#'     distinct null value present among the displayed coefficients.
#'   * `"median"` — the median of the coefficients actually plotted. A convenience
#'     approximation to `"null"` for sources without a `null_hypothesis` column; it is
#'     biased when `top_n` subsets to extreme taxa and pools across coefficients.
#'   * a single numeric x-intercept, or `NA` / `NULL` to omit the line.
#'   Significance stars always key off the reported q-values, never the line position.
#' @param reference_line_color Color for the vertical dashed reference line(s)
#'   (default: `"grey50"`).
#' @param taxon_groups Optional named character vector mapping taxon -> group
#'   label (e.g. genus -> phylum). When supplied for a single-coefficient plot,
#'   points are colored by group and a group legend is shown.
#' @param group_colors Optional named character vector mapping group label ->
#'   color, used together with `taxon_groups`.
#' @param legend_levels Optional character vector of group labels to pin the group
#'   legend to, used together with `taxon_groups`. Every level gets a key even when
#'   the panel plots no taxon from it, so that sibling panels of a multi-panel figure
#'   emit identical legends that `patchwork::plot_layout(guides = "collect")` can
#'   merge into one. Levels are shown in the given order.
#' @param title Optional plot title.
#' @param taxon_labeller Function applied to the taxon (y-axis) labels for display
#'   only, e.g. to shorten over-long names (default: `identity`).
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' plot_maaslin3_forest(maaslin3_result, "TypeStent")
#' plot_maaslin3_forest(
#'   "all_results.tsv",
#'   c("GroupBenign Stricture", "GroupIntra Biliary Malignancy"),
#'   model_type = "joint", top_n = 15,
#'   coefficient_labels = c(
#'     "GroupBenign Stricture" = "Benign Stricture",
#'     "GroupIntra Biliary Malignancy" = "Intra Biliary Malignancy"
#'   )
#' )
#' }
#'
#' @export
plot_maaslin3_forest <- function(
  maaslin3_input,
  coefficient_names,
  model_type = "abundance",
  taxa = NULL,
  top_n = 20,
  q_threshold = 0.05,
  sort_by = c("coef", "qval", "taxon"),
  color_by_coefficient = TRUE,
  show_q_legend = TRUE,
  coefficient_labels = NULL,
  color_values = NULL,
  reference_line = "zero",
  reference_line_color = "grey50",
  taxon_groups = NULL,
  group_colors = NULL,
  legend_levels = NULL,
  italicize_taxa = TRUE,
  color_legend_title = "Coefficient",
  title = NULL,
  base_size = 11,
  taxon_labeller = identity
) {

  sort_by <- match.arg(sort_by)

  if (length(coefficient_names) < 1) {
    stop("coefficient_names must contain at least one coefficient.")
  }

  # Prepare data
  df <- prepare_maaslin3_data(
    maaslin3_input,
    coefficient_names = coefficient_names,
    model_type = model_type,
    taxa = taxa
  )

  if (nrow(df) == 0) {
    stop("No data available after filtering.")
  }

  # Select top_n taxa if no explicit taxa provided
  if (is.null(taxa) && !is.null(top_n)) {
    # Rank by minimum q-value across coefficients
    min_q_per_taxon <- stats::aggregate(qval ~ taxon, data = df, FUN = min, na.rm = TRUE)
    min_q_per_taxon <- min_q_per_taxon[order(min_q_per_taxon$qval), ]
    keep_taxa <- utils::head(min_q_per_taxon$taxon, top_n)
    df <- df[df$taxon %in% keep_taxa, , drop = FALSE]
  }

  # Compute 95% CI
  df$ci_lower <- df$coef - 1.96 * df$stderr
  df$ci_upper <- df$coef + 1.96 * df$stderr

  # Significance indicator (drives the star annotation)
  df$is_significant <- df$qval < q_threshold
  df$sig_label <- ifelse(df$is_significant, "*", "")

  # Apply coefficient label renaming, preserving input order for factor levels
  if (!is.null(coefficient_labels) && "coefficient_name" %in% names(df)) {
    matched <- match(df$coefficient_name, names(coefficient_labels))
    has_match <- !is.na(matched)
    df$coefficient_name[has_match] <- coefficient_labels[matched[has_match]]
    # Set factor levels in the order of coefficient_labels (user-specified order)
    coef_level_order <- unname(coefficient_labels)
    df$coefficient_name <- factor(df$coefficient_name, levels = coef_level_order)
  } else {
    # Preserve order of coefficient_names input
    df$coefficient_name <- factor(df$coefficient_name, levels = unique(df$coefficient_name[match(coefficient_names, df$coefficient_name)]))
  }

  # Sort taxa
  if (sort_by == "coef") {
    # Sort by max absolute coefficient across all shown coefficients
    max_abs_coef <- stats::aggregate(coef ~ taxon, data = df, FUN = function(x) max(abs(x)))
    taxon_order <- max_abs_coef$taxon[order(max_abs_coef$coef)]
  } else if (sort_by == "qval") {
    min_q <- stats::aggregate(qval ~ taxon, data = df, FUN = min, na.rm = TRUE)
    taxon_order <- min_q$taxon[order(min_q$qval, decreasing = TRUE)]
  } else {
    taxon_order <- sort(unique(df$taxon), decreasing = TRUE)
  }
  df$taxon <- factor(df$taxon, levels = taxon_order)

  # Default title
  if (is.null(title)) {
    model_label <- if (model_type == "joint") "Joint" else tools::toTitleCase(model_type)
    title <- paste0("MaAsLin3 Forest Plot (", model_label, ")")
  }

  # Build plot
  single_coef <- length(coefficient_names) == 1
  dodge <- ggplot2::position_dodge(width = 0.6)

  if (single_coef && !is.null(taxon_groups)) {
    # Color each taxon's point by its group (e.g. phylum), with a group legend.
    unclassified_label <- "(Unclassified)"
    df$group <- unname(taxon_groups[as.character(df$taxon)])
    has_unclassified <- anyNA(df$group)
    df$group[is.na(df$group)] <- unclassified_label

    # Without this, the factor() below would map the fallback straight back to NA (grey points,
    # no legend key) whenever the caller's legend_levels omits it.
    if (has_unclassified && !is.null(legend_levels) &&
        !unclassified_label %in% legend_levels) {
      legend_levels <- c(legend_levels, unclassified_label)
    }

    # A legend key only gets a glyph for a group the layer actually plots. Panels of a
    # multi-panel figure show different taxa, so their group legends end up with
    # different keys and patchwork's guides = "collect" cannot merge them, however
    # carefully the scales are pinned. `legend_levels` pads the data with one
    # never-drawn row (coef = NA) per missing group, so every panel emits the same
    # complete key set and the legends merge into one.
    if (!is.null(legend_levels)) {
      missing_groups <- setdiff(legend_levels, unique(df$group))
      if (length(missing_groups) > 0) {
        pad <- df[rep(1L, length(missing_groups)), , drop = FALSE]
        pad$group <- missing_groups
        pad$coef <- NA_real_
        pad$ci_lower <- NA_real_
        pad$ci_upper <- NA_real_
        pad$sig_label <- ""
        df <- rbind(df, pad)
      }
      df$group <- factor(df$group, levels = legend_levels)
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = coef, y = taxon, color = group, fill = group
    )) +
      ggplot2::geom_pointrange(
        ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
        shape = 21, size = 0.4, fatten = 2.5, stroke = 0.4, na.rm = TRUE
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = coef, y = taxon, label = sig_label),
        vjust = -0.07, size = 3.5, fontface = "bold", show.legend = FALSE,
        na.rm = TRUE
      ) +
      ggplot2::labs(color = color_legend_title, fill = color_legend_title)

    if (!is.null(group_colors)) {
      # Give the fallback level its own key rather than letting it fall through to na.value,
      # which would drop it from the legend entirely.
      if (has_unclassified && !unclassified_label %in% names(group_colors)) {
        group_colors <- c(group_colors, stats::setNames("grey60", unclassified_label))
      }
      limits <- if (is.null(legend_levels)) ggplot2::waiver() else legend_levels
      p <- p +
        ggplot2::scale_color_manual(values = group_colors, limits = limits,
                                    drop = FALSE, na.value = "grey60") +
        ggplot2::scale_fill_manual(values = group_colors, limits = limits,
                                   drop = FALSE, na.value = "grey60")
    }
  } else if (single_coef) {
    point_color <- if (!is.null(color_values)) color_values[[1]] else "black"

    p <- ggplot2::ggplot(df, ggplot2::aes(x = coef, y = taxon)) +
      ggplot2::geom_pointrange(
        ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
        shape = 21, size = 0.4, fatten = 2.5,
        color = point_color, fill = point_color, stroke = 0.4
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = coef, y = taxon, label = sig_label),
        color = point_color, vjust = -0.07, size = 3.5, fontface = "bold",
        show.legend = FALSE
      )
  } else {
    if (color_by_coefficient) {
      coef_levels <- unique(df$coefficient_name)
      if (is.null(color_values)) {
        coef_colors <- scales::hue_pal()(length(coef_levels))
        names(coef_colors) <- coef_levels
      } else {
        coef_colors <- color_values
      }

      df$fill_color <- coef_colors[as.character(df$coefficient_name)]

      p <- ggplot2::ggplot(df, ggplot2::aes(
        x = coef, y = taxon, color = coefficient_name, group = coefficient_name
      )) +
        ggplot2::geom_pointrange(
          ggplot2::aes(xmin = ci_lower, xmax = ci_upper, fill = fill_color),
          shape = 21, size = 0.4, fatten = 2.5, stroke = 0.4,
          position = dodge
        ) +
        ggplot2::geom_text(
          ggplot2::aes(x = coef, y = taxon, label = sig_label, group = coefficient_name),
          position = dodge, vjust = -0.07, size = 3.5, fontface = "bold",
          show.legend = FALSE
        ) +
        ggplot2::scale_fill_identity() +
        ggplot2::labs(color = color_legend_title) +
        ggplot2::guides(color = ggplot2::guide_legend(
          override.aes = list(fill = unname(coef_colors[levels(df$coefficient_name)]), shape = 21)
        ))

      if (!is.null(color_values)) {
        p <- p + ggplot2::scale_color_manual(values = color_values, name = color_legend_title) +
          ggplot2::guides(color = ggplot2::guide_legend(
            override.aes = list(fill = unname(coef_colors[levels(df$coefficient_name)]), shape = 21)
          ))
      }
    } else {
      df$fill_color <- "black"

      p <- ggplot2::ggplot(df, ggplot2::aes(x = coef, y = taxon, group = coefficient_name)) +
        ggplot2::geom_pointrange(
          ggplot2::aes(xmin = ci_lower, xmax = ci_upper, fill = fill_color),
          shape = 21, size = 0.4, fatten = 2.5, stroke = 0.4,
          position = dodge
        ) +
        ggplot2::geom_text(
          ggplot2::aes(x = coef, y = taxon, label = sig_label, group = coefficient_name),
          position = dodge, vjust = -0.07, size = 3.5, fontface = "bold",
          show.legend = FALSE
        ) +
        ggplot2::scale_fill_identity()
    }
  }

  # Resolve where the reference line(s) are drawn. `df` still carries the NA padding
  # rows added for legend alignment, so medians/nulls are taken with na.rm = TRUE.
  # "null" can yield several x-intercepts (one per coefficient's median-comparison null).
  ref_x <- if (is.null(reference_line) || (length(reference_line) == 1 && is.na(reference_line))) {
    numeric(0)
  } else if (is.numeric(reference_line)) {
    reference_line[[1]]
  } else {
    switch(match.arg(reference_line, c("zero", "null", "median")),
      zero   = 0,
      null   = {
        if (!"null_hypothesis" %in% names(df)) {
          warning("reference_line = 'null' but no null_hypothesis column is available; ",
                  "falling back to 0.")
          0
        } else {
          sort(unique(round(df$null_hypothesis[!is.na(df$coef)], 8)))
        }
      },
      median = stats::median(df$coef, na.rm = TRUE)
    )
  }

  if (length(ref_x) > 0) {
    p <- p +
      ggplot2::geom_vline(xintercept = ref_x, linetype = "dashed",
                          color = reference_line_color, linewidth = 0.5)
  }

  p <- p +
    ggplot2::scale_y_discrete(labels = taxon_labeller) +
    ggplot2::labs(
      x = "Beta Coefficient (95% CI)",
      y = NULL,
      title = title
    )

  if (italicize_taxa) {
    p <- p + ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"))
  }

  p
}
