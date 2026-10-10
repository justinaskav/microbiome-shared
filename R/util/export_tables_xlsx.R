#' Export Abundance Data to Excel
#'
#' For each taxonomic rank, aggregates the TSE (keeping all taxa, merging
#' NAs/Incertae Sedis via `aggregate_taxa()`), transposes the chosen assay
#' so that rows are samples and columns are taxa, and writes each rank
#' as a separate sheet in an xlsx workbook.
#'
#' @param tse A TreeSummarizedExperiment object (counts assay required).
#' @param ranks Character vector of taxonomic ranks to export (e.g.,
#'   `c("family", "genus")`).
#' @param output_path File path for the output `.xlsx` file.
#' @param assay_type Which assay to export: `"relabundance"` (default) or
#'   `"counts"`.
#' @param unclassified Passed to `aggregate_taxa()` (default: `"group"`).
#' @param ... Additional arguments passed to `aggregate_taxa()`.
#'
#' @return The output path (invisibly).
export_abundance_xlsx <- function(tse, ranks, output_path,
                                  assay_type = "relabundance",
                                  unclassified = "group", ...) {
  assay_type <- match.arg(assay_type, c("relabundance", "counts"))
  wb <- openxlsx::createWorkbook()

  for (rank in ranks) {
    tse_agg <- mia::transformAssay(tse, method = "relabundance") |>
      aggregate_taxa(
        rank = rank, top = NULL,
        unclassified = unclassified, ...
      )

    mat <- as.data.frame(
      t(SummarizedExperiment::assay(tse_agg, assay_type))
    )
    mat <- cbind(Sample = rownames(mat), mat)

    openxlsx::addWorksheet(wb, sheetName = rank)
    openxlsx::writeData(wb, sheet = rank, x = mat)
  }

  openxlsx::saveWorkbook(wb, output_path, overwrite = TRUE)
  invisible(output_path)
}
