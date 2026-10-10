#' Plot Taxonomic Composition as Heatmap
#'
#' Creates a heatmap visualization of taxonomic composition using ComplexHeatmap.
#' The TSE object should already be aggregated at the desired taxonomic level.
#'
#' @param tse A TreeSummarizedExperiment object (already aggregated at desired level).
#' @param taxa Optional character vector of taxa names to filter. If NULL, uses all taxa.
#' @param taxa_column Column name in rowData containing taxa names (default: "genus").
#'   Used for filtering when taxa parameter is provided.
#' @param aggregate_by Column name in colData to aggregate samples by (default: NULL), via
#'   `aggregate_to_patient()` (mean of per-sample relative abundances). Requires
#'   R/beta_diversity/aggregate_to_patient.R. Set to NULL to skip aggregation.
#' @param annotation_var Column name in colData for column annotation (default: NULL).
#' @param annotation_colors Optional named list of colors for annotation.
#'   If NULL, default colors are used.
#' @param split_by Column name to split heatmap columns by (default: same as annotation_var).
#' @param transform_method Transformation method (default: "rclr"), applied to all taxa before
#'   `taxa` subsetting. Set to NULL to skip.
#' @param standardize Logical, whether to z-score each taxon across columns after transformation
#'   (default: TRUE).
#' @param cluster_cols Logical, whether to cluster columns (default: FALSE).
#' @param cluster_rows Logical, whether to cluster rows (default: FALSE).
#' @param color_palette Either an RColorBrewer palette name (e.g., "YlOrRd") or a
#'   character vector of colors (e.g., `viridisLite::viridis(100)`).
#' @param show_rownames Logical, whether to show row names (default: TRUE).
#' @param show_colnames Logical, whether to show column names (default: FALSE).
#' @param fontsize_row Font size for row names (default: 8.5).
#' @param fontsize_col Font size for column names (default: 6).
#' @param annotation_height Height of column annotation bar in mm (default: NULL,
#'   uses ComplexHeatmap default).
#' @param legend_title Title for heatmap legend (default: "Abundance").
#'
#' @return A ComplexHeatmap Heatmap object.
#'
#' @examples
#' \dontrun{
#' # Plot specific genera
#' genera_of_interest <- c("Fusobacterium", "Prevotella", "Bacteroides")
#' ht <- plot_taxonomic_heatmap(
#'   tse = altExp(tses$gut, "genus_agglomerate"),
#'   taxa = genera_of_interest,
#'   aggregate_by = "Patient",
#'   annotation_var = "Diagnosis"
#' )
#' ComplexHeatmap::draw(ht)
#' }
plot_taxonomic_heatmap <- function(
    tse,
    taxa = NULL,
    taxa_column = "genus",
    aggregate_by = NULL,
    annotation_var = NULL,
    annotation_colors = NULL,
    split_by = NULL,
    transform_method = "rclr",
    standardize = TRUE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    color_palette = "YlOrRd",
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 8.5,
    fontsize_col = 6,
    annotation_height = NULL,
    legend_title = NULL
) {
  # Validate inputs
  if (!inherits(tse, "TreeSummarizedExperiment") && !inherits(tse, "SummarizedExperiment")) {
    stop("tse must be a TreeSummarizedExperiment or SummarizedExperiment object")
  }

  # Set split_by to annotation_var if not specified
  if (is.null(split_by)) {
    split_by <- annotation_var
  }

  # Aggregate by specified variable (e.g., Patient) as the mean of per-sample proportions, so a
  # deeply sequenced sample does not dominate its group (summing counts would weight by depth).
  if (!is.null(aggregate_by) && aggregate_by != "") {
    tse <- aggregate_to_patient(tse, patient_var = aggregate_by, warn = FALSE)
  }

  # Transform on the FULL table, then subset: CLR-type transforms centre each sample on all of
  # its features, so transforming after subsetting would centre on the selected taxa only.
  # rclr is invariant to per-sample scaling, so relabundance (mean proportions after
  # aggregation) and counts give the same result.
  if (!is.null(transform_method) && transform_method != "") {
    input_assay <- if ("relabundance" %in% SummarizedExperiment::assayNames(tse)) "relabundance" else "counts"
    tse <- mia::transformAssay(tse, method = transform_method, assay.type = input_assay)
    assay_name <- transform_method

    # Z-score each taxon across columns. mia's default MARGIN is "samples", which would
    # standardise each column across taxa instead.
    if (standardize) {
      assay_name <- paste0(transform_method, "_z")
      tse <- mia::transformAssay(tse, method = "standardize", assay.type = transform_method,
                                 MARGIN = "features", name = assay_name)
    }
  } else {
    assay_name <- "counts"
  }

  # Filter to specific taxa if provided
  if (!is.null(taxa) && length(taxa) > 0) {
    rd <- SummarizedExperiment::rowData(tse)
    taxa_mask <- if (taxa_column %in% colnames(rd)) rd[[taxa_column]] %in% taxa else rownames(tse) %in% taxa
    tse <- tse[taxa_mask, ]
    if (nrow(tse) == 0) {
      stop("No matching taxa found. Check that taxa names match rowData or rownames.")
    }
  }

  # Determine legend title
  if (is.null(legend_title)) {
    legend_title <- if (standardize) paste0("Standardized ", transform_method) else transform_method
  }

  # Extract assay matrix
  h_assay <- SummarizedExperiment::assay(tse, assay_name)

  # Create column annotation if annotation_var is provided
  if (!is.null(annotation_var) && annotation_var != "") {
    if (!(annotation_var %in% colnames(SummarizedExperiment::colData(tse)))) {
      warning(paste("Annotation variable", annotation_var, "not found in colData. Skipping annotation."))
      column_ha <- NULL
      column_split <- NULL
    } else {
      # Default colors if not provided
      if (is.null(annotation_colors)) {
        annotation_values <- unique(SummarizedExperiment::colData(tse)[[annotation_var]])
        n_colors <- length(annotation_values)
        default_colors <- RColorBrewer::brewer.pal(max(3, min(n_colors, 8)), "Set2")[1:n_colors]
        names(default_colors) <- annotation_values
        annotation_colors <- list()
        annotation_colors[[annotation_var]] <- default_colors
      }

      ha_args <- list(
        df = data.frame(SummarizedExperiment::colData(tse)[, annotation_var, drop = FALSE]),
        gp = grid::gpar(fontsize = 6),
        annotation_name_gp = grid::gpar(fontsize = 10),
        col = annotation_colors
      )
      if (!is.null(annotation_height)) {
        ha_args$simple_anno_size <- grid::unit(annotation_height, "mm")
      }
      column_ha <- do.call(ComplexHeatmap::HeatmapAnnotation, ha_args)

      # Column split
      if (!is.null(split_by) && split_by %in% colnames(SummarizedExperiment::colData(tse))) {
        column_split <- SummarizedExperiment::colData(tse)[[split_by]]
      } else {
        column_split <- NULL
      }
    }
  } else {
    column_ha <- NULL
    column_split <- NULL
  }

  # Resolve heatmap colors
  if (is.character(color_palette) && length(color_palette) == 1) {
    heatmap_colors <- RColorBrewer::brewer.pal(3, color_palette)
  } else {
    heatmap_colors <- color_palette
  }

  # Create heatmap
  ht <- ComplexHeatmap::pheatmap(
    h_assay,
    top_annotation = column_ha,
    column_split = column_split,
    cluster_cols = cluster_cols,
    cluster_rows = cluster_rows,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    color = heatmap_colors,
    border_color = "grey97",
    heatmap_legend_param = list(title = legend_title)
  )

  return(ht)
}
