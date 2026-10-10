#' Summarize unclassified and ambiguous taxa across taxonomic ranks
#'
#' For each taxonomic rank in a TreeSummarizedExperiment object, counts the
#' number of features that are NA, empty, bare prefixes (e.g. "g__"),
#' uncultured, contain "Incertae Sedis", or contain "unclassified_" composites.
#' Columns where all ranks show zero are dropped automatically.
#'
#' @param tse A TreeSummarizedExperiment object.
#' @param ranks Character vector of taxonomic ranks to summarize. If NULL
#'   (default), auto-detects standard ranks present in rowData.
#' @param drop_empty Logical. Drop pattern columns that are all "-" across
#'   every rank (default: TRUE).
#' @return A data frame with one row per rank and columns for each detected
#'   pattern.
summarize_unclassified_taxa <- function(tse, ranks = NULL, drop_empty = TRUE) {
  rd <- SummarizedExperiment::rowData(tse)

  # Auto-detect ranks if not provided
  if (is.null(ranks)) {
    standard_ranks <- c("kingdom", "domain", "phylum", "class", "order",
                        "family", "genus", "species")
    ranks <- intersect(standard_ranks, colnames(rd))
  }

  n_total <- nrow(rd)

  fmt <- function(n) {
    if (n == 0) return("-")
    pct <- round(n / n_total * 100, 1)
    paste0(format(n, big.mark = ","), " (", pct, "%)")
  }

  results <- lapply(ranks, function(r) {
    vals <- rd[[r]]
    non_na <- vals[!is.na(vals)]

    data.frame(
      Rank = r,
      Total = format(n_total, big.mark = ","),
      NAs = fmt(sum(is.na(vals))),
      Empty = fmt(sum(non_na == "")),
      `Bare prefix` = fmt(sum(grepl("^[a-z]__$", non_na))),
      Uncultured = fmt(sum(grepl("^uncultured$", non_na, ignore.case = TRUE))),
      `Incertae Sedis` = fmt(sum(grepl("Incertae.Sedis", non_na, ignore.case = TRUE))),
      `unclassified_*` = fmt(sum(grepl("^unclassified_", non_na, ignore.case = TRUE))),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, results)

  if (drop_empty) {
    # Drop pattern columns (not Rank/Total) where every row is "-"
    pattern_cols <- setdiff(names(df), c("Rank", "Total"))
    keep <- vapply(pattern_cols, function(col) !all(df[[col]] == "-"), logical(1))
    df <- df[, c("Rank", "Total", pattern_cols[keep]), drop = FALSE]
  }

  df
}
