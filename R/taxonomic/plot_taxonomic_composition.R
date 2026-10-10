#' Plot Taxonomic Composition as Stacked Barplot
#'
#' Creates a stacked 100% barplot showing taxonomic composition per sample.
#' The TSE object should already be aggregated at the desired taxonomic level
#' (via `aggregate_taxa()`).
#'
#' @param tse A TreeSummarizedExperiment object (already aggregated at desired level).
#' @param assay_type Assay type to use (default: "relabundance").
#' @param x_var Column name in colData for x-axis (default: "Sample").
#' @param fill_var Column name for fill aesthetic (default: "FeatureID").
#' @param facet_var Column name(s) in colData for faceting. If length 1, uses
#'   `facet_grid(. ~ var)` with free_x and free space. If length 2, uses
#'   `facet_grid(var1 ~ var2)`. Set to NULL for no faceting.
#' @param palette Named or unnamed character vector of colors. If NULL, a
#'   default palette is generated from Set2 + Set3.
#' @param remainder_labels Character vector of labels treated as remainder/other
#'   categories (default: `c("(Others)", "(Remainder)")`). These are always
#'   ordered last and colored grey.
#' @param remainder_color Color for remainder taxa (default: "#D4C5A9").
#' @param unclassified_labels Character vector of labels treated as unclassified
#'   categories (default: `c("(Unclassified)")`). Ordered after remainder.
#' @param unclassified_color Color for unclassified taxa (default: "#C0C0C0").
#' @param legend_title Title for the fill legend. If NULL, no title is shown.
#' @param y_label Label for y-axis (default: "Relative Abundance").
#' @param x_angle Angle for x-axis text (default: 45).
#' @param add_cols Additional colData columns to include in melted data frame.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_taxonomic_composition(
#'   altExp(tses$rrna, "phylum_aggregate"),
#'   facet_var = "Group",
#'   palette = my_palette,
#'   legend_title = "Phylum"
#' )
#' }
plot_taxonomic_composition <- function(
    tse,
    assay_type = "relabundance",
    x_var = "Sample",
    fill_var = "FeatureID",
    facet_var = NULL,
    palette = NULL,
    remainder_labels = c("(Others)", "(Remainder)"),
    remainder_color = "#D4C5A9",
    unclassified_labels = c("(Unclassified)"),
    unclassified_color = "#C0C0C0",
    legend_title = NULL,
    y_label = "Relative Abundance",
    x_angle = 45,
    add_cols = NULL,
    italicize_taxa = TRUE
) {
  if (!inherits(tse, "TreeSummarizedExperiment") && !inherits(tse, "SummarizedExperiment")) {
    stop("tse must be a TreeSummarizedExperiment or SummarizedExperiment object")
  }

  if (!(assay_type %in% SummarizedExperiment::assayNames(tse))) {
    stop(paste("Assay", assay_type, "not found. Available:",
               paste(SummarizedExperiment::assayNames(tse), collapse = ", ")))
  }

  # Build column selection for meltSE
  melt_cols <- unique(c(x_var, facet_var, add_cols))
  melt_cols <- melt_cols[!is.null(melt_cols) & melt_cols != ""]

  # Melt TSE to long format
  plot_data <- mia::meltSE(
    tse,
    assay.type = assay_type,
    rownames = "taxon",
    add.col = melt_cols
  )

  # Labels to exclude from abundance sorting
  special_labels <- c(remainder_labels, unclassified_labels)

  # Order taxa by mean abundance (descending), special labels always last
  all_taxa <- unique(as.character(plot_data[[fill_var]]))
  remainder_present <- intersect(remainder_labels, all_taxa)
  unclassified_present <- intersect(unclassified_labels, all_taxa)

  mean_abund <- plot_data |>
    dplyr::filter(!(.data[[fill_var]] %in% special_labels)) |>
    dplyr::group_by(.data[[fill_var]]) |>
    dplyr::summarise(mean_abund = mean(.data[[assay_type]], na.rm = TRUE), .groups = "drop")

  # Factor levels: descending abundance, remainder then unclassified last.
  named_by_abund <- as.character(
    mean_abund[[fill_var]][order(mean_abund$mean_abund, decreasing = TRUE)]
  )
  ordered_levels <- c(named_by_abund, remainder_present, unclassified_present)
  plot_data[[fill_var]] <- factor(plot_data[[fill_var]], levels = ordered_levels)

  # Build named color vector
  named_taxa <- named_by_abund

  # Default palette from Set2 + Set3
  if (is.null(palette)) {
    palette <- c(
      RColorBrewer::brewer.pal(8, "Set2"),
      RColorBrewer::brewer.pal(12, "Set3")
    )
  }

  color_map <- stats::setNames(palette[seq_along(named_taxa)], named_taxa)
  for (rl in remainder_present) {
    color_map[rl] <- remainder_color
  }
  for (ul in unclassified_present) {
    color_map[ul] <- unclassified_color
  }

  # Italicize taxa names in legend (except special labels) using markdown
  if (italicize_taxa) {
    italicize_label <- function(x) {
      ifelse(x %in% special_labels, x, paste0("*", x, "*"))
    }
    new_levels <- italicize_label(levels(plot_data[[fill_var]]))
    levels(plot_data[[fill_var]]) <- new_levels
    names(color_map) <- italicize_label(names(color_map))
  }

  # Build plot
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data[[x_var]], y = .data[[assay_type]], fill = .data[[fill_var]])
  ) +
    ggplot2::geom_bar(position = "fill", stat = "identity", linewidth = 0) +
    ggplot2::labs(x = x_var, y = y_label) +
    ggplot2::scale_x_discrete(labels = scales::label_wrap(10)) +
    # ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 100, accuracy = 1)) +
    ggplot2::scale_fill_manual(values = color_map, name = legend_title) +
    ggplot2::theme(
      legend.text = if (italicize_taxa) ggtext::element_markdown() else ggplot2::element_text()
    )

  if (!is.null(x_angle) && x_angle != 0) {
    p <- p + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = x_angle, hjust = 1)
    )
  }

  # Faceting
  if (!is.null(facet_var)) {
    if (length(facet_var) == 1) {
      p <- p + ggplot2::facet_grid(
        stats::as.formula(paste(". ~", facet_var)),
        scales = "free_x", space = "free_x"
      )
    } else if (length(facet_var) == 2) {
      p <- p + ggplot2::facet_grid(
        stats::as.formula(paste(facet_var[1], "~", facet_var[2])),
        scales = "free_x", space = "free_x"
      )
    }
  }

  p
}
