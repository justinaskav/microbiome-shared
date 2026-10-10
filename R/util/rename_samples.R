#' Rename Samples in a TreeSummarizedExperiment
#'
#' Updates column names and the Sample column in colData based on a named
#' character vector mapping old names to new names. SampleID remains unchanged.
#'
#' @param tse A TreeSummarizedExperiment or SummarizedExperiment object.
#' @param rename_map A named character vector where names are old sample names
#'   and values are new sample names (e.g., `c("old1" = "new1", "old2" = "new2")`).
#'
#' @return The TSE object with updated column names and Sample column.
rename_samples <- function(tse, rename_map) {
  # One lookup per name, so a map like c(a = "b", b = "c") renames a -> b and b -> c rather
  # than chaining a -> b -> c.
  remap <- function(x) {
    hit <- match(x, names(rename_map))
    ifelse(is.na(hit), x, unname(rename_map)[hit])
  }
  colnames(tse) <- remap(colnames(tse))  # also renames assay columns and colData rows
  SummarizedExperiment::colData(tse)$Sample <- remap(as.character(SummarizedExperiment::colData(tse)$Sample))
  tse
}
