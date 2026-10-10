#' Plot Taxonomic Composition as Box-Jitter Plot
#'
#' Creates a faceted box-jitter plot showing taxonomic composition for specified taxa.
#' The TSE object should already be aggregated at the desired taxonomic level.
#'
#' @param tse A TreeSummarizedExperiment object (already aggregated at desired level).
#' @param taxa Optional character vector of taxa names to filter. If NULL, uses all taxa.
#' @param x_var Column name in colData for x-axis (required).
#' @param color_var Column name in colData for point color (default: NULL).
#' @param fill_var Column name in colData for boxplot fill (default: NULL).
#' @param shape_var Column name in colData for point shape (default: NULL).
#' @param facet_var Column to facet by (default: "FeatureID" for per-taxon facets).
#' @param assay_type Assay type to use (default: "relabundance").
#' @param log_scale Logical, whether to use a logarithmic y-axis (default: TRUE). A pseudo-log
#'   (asinh) scale is used, linear below the smallest non-zero value and log10 above it, so
#'   samples where the taxon is absent stay in the boxplots (plotted at 0) instead of being
#'   dropped as -Inf, which would make the boxes describe detected samples only.
#' @param y_label Label for y-axis (default: "Relative Abundance").
#' @param jitter_width Width for jittered points (default: 0.1).
#' @param jitter_alpha Alpha for jittered points (default: 0.4).
#' @param point_size Size for jittered points (default: 1.2).
#' @param boxplot_alpha Alpha for boxplots (default: 0.6).
#' @param jitter_seed Seed for the jitter offsets (the global RNG state is left untouched).
#' @param label_var Optional colData column used to label non-zero points.
#' @param x_order Optional level order for `x_var`.
#' @param italicize_taxa Italicise facet (taxon) strip labels (default: TRUE).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # Plot specific genera
#' genera_of_interest <- c("Fusobacterium", "Prevotella", "Bacteroides")
#' plot_taxonomic_boxjitter(
#'   tse = altExp(tses$gut, "genus_agglomerate"),
#'   taxa = genera_of_interest,
#'   x_var = "Diagnosis",
#'   color_var = "Site"
#' )
#'
#' # Plot all taxa in TSE (for small number of taxa)
#' plot_taxonomic_boxjitter(
#'   tse = altExp(tses$gut, "phylum_aggregate"),
#'   x_var = "Group",
#'   taxa = NULL  # Use all
#' )
#' }
plot_taxonomic_boxjitter <- function(
    tse,
    taxa = NULL,
    x_var,
    color_var = NULL,
    fill_var = NULL,
    shape_var = NULL,
    facet_var = "FeatureID",
    assay_type = "relabundance",
    log_scale = TRUE,
    y_label = "Relative Abundance (%)",
    jitter_width = 0.1,
    jitter_alpha = 0.4,
    point_size = 1.2,
    boxplot_alpha = 0.6,
    jitter_seed = 42,
    label_var = NULL,
    x_order = NULL,
    italicize_taxa = TRUE
) {
  # Validate inputs
  if (!inherits(tse, "TreeSummarizedExperiment") && !inherits(tse, "SummarizedExperiment")) {
    stop("tse must be a TreeSummarizedExperiment or SummarizedExperiment object")
  }

  if (missing(x_var) || is.null(x_var)) {
    stop("x_var is required")
  }

  if (!(assay_type %in% SummarizedExperiment::assayNames(tse))) {
    stop(paste("Assay", assay_type, "not found. Available:",
               paste(SummarizedExperiment::assayNames(tse), collapse = ", ")))
  }

  # Build column selection for meltSE
  add_cols <- unique(c(x_var, color_var, fill_var, shape_var, label_var))
  add_cols <- add_cols[!is.na(add_cols) & add_cols != ""]

  # Melt the TSE object
  melted_df <- mia::meltSE(
    tse,
    assay.type = assay_type,
    rownames = "taxon",
    add.col = add_cols
  )

  # Filter to specific taxa if provided
  if (!is.null(taxa) && length(taxa) > 0) {
    melted_df <- dplyr::filter(melted_df, FeatureID %in% taxa)

    if (nrow(melted_df) == 0) {
      warning("No matching taxa found. Check that taxa names match FeatureID values.")
      return(ggplot2::ggplot() + ggplot2::theme_void() +
               ggplot2::labs(title = "No matching taxa found"))
    }
  }

  # Apply custom x ordering if provided
  if (!is.null(x_order)) {
    melted_df[[x_var]] <- factor(melted_df[[x_var]], levels = x_order)
  } else {
    melted_df[[x_var]] <- as.factor(melted_df[[x_var]])
  }

  # Pre-compute jittered x positions, reproducibly and without touching the caller's RNG stream
  melted_df$x_jittered <- as.numeric(melted_df[[x_var]]) +
    withr::with_seed(jitter_seed, stats::runif(nrow(melted_df), -jitter_width, jitter_width))

  # Build the plot using jittered positions
  p <- ggplot2::ggplot(melted_df, ggplot2::aes(x = .data[[x_var]], y = .data[[assay_type]]))

  # Add points at jittered positions
  point_aes <- if (!is.null(color_var) && !is.null(shape_var)) {
    ggplot2::aes(x = x_jittered, color = .data[[color_var]], shape = .data[[shape_var]])
  } else if (!is.null(color_var)) {
    ggplot2::aes(x = x_jittered, color = .data[[color_var]])
  } else if (!is.null(shape_var)) {
    ggplot2::aes(x = x_jittered, shape = .data[[shape_var]])
  } else {
    ggplot2::aes(x = x_jittered)
  }
  p <- p + ggplot2::geom_point(point_aes, alpha = jitter_alpha, size = point_size)

  # Add boxplots
  if (!is.null(fill_var)) {
    p <- p + ggplot2::geom_boxplot(
      ggplot2::aes(fill = .data[[fill_var]]),
      alpha = boxplot_alpha,
      outlier.shape = NA
    )
  } else {
    p <- p + ggplot2::geom_boxplot(
      alpha = boxplot_alpha,
      outlier.shape = NA
    )
  }

  # Add labels
  p <- p + ggplot2::labs(
    x = x_var,
    y = y_label
  )

  # Apply y-axis scale with percent formatting
  if (log_scale) {
    values <- melted_df[[assay_type]]
    sigma <- min(values[values > 0], na.rm = TRUE)
    decades <- 10^seq(floor(log10(sigma)), ceiling(log10(max(values, na.rm = TRUE))))
    p <- p + ggplot2::scale_y_continuous(
      transform = scales::pseudo_log_trans(sigma = sigma, base = 10),
      breaks = c(0, decades), labels = scales::percent_format(scale = 100)
    )
  } else {
    p <- p + ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 100))
  }

  # Wrap x-axis labels
  p <- p + ggplot2::scale_x_discrete(labels = scales::label_wrap(10))

  # Add structural theme adjustments only (decorative theming set by caller via theme_set())
  theme_args <- list(panel.spacing = ggplot2::unit(1, "lines"))
  if (italicize_taxa) {
    theme_args$strip.text <- ggplot2::element_text(face = "italic")
  }
  p <- p + do.call(ggplot2::theme, theme_args)

  # Add faceting if requested
  if (!is.null(facet_var) && facet_var != "") {
    # Use free y-axis scales when not using log scale to better show range differences
    scales_param <- if (log_scale) "fixed" else "free_y"
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)), scales = scales_param)
  }

  # Add labels if requested
  if (!is.null(label_var) && label_var != "") {
    label_aes <- if (!is.null(color_var)) {
      ggplot2::aes(x = x_jittered, label = .data[[label_var]], color = .data[[color_var]])
    } else {
      ggplot2::aes(x = x_jittered, label = .data[[label_var]])
    }
    p <- p + ggrepel::geom_label_repel(
      data = ~ dplyr::filter(.x, .data[[assay_type]] > 0),
      label_aes,
      size = 2, max.overlaps = 5, segment.color = NA,
      label.padding = ggplot2::unit(0.12, "lines"), label.size = 0.3,
      fill = scales::alpha("white", 0.85),
      label.r = ggplot2::unit(0, "lines"), show.legend = FALSE
    )
  }

  return(p)
}
