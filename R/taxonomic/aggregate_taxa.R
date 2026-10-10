#' Aggregate taxa to a specified rank, with Incertae Sedis, unclassified, and top-N handling
#'
#' @param tse_obj TreeSummarizedExperiment object
#' @param rank Taxonomic rank to aggregate (e.g., "genus")
#' @param top Number of top taxa to keep (others grouped as remainder).
#'   If NULL, all taxa are returned without grouping
#' @param incertae_sedis Character. How to handle Incertae Sedis taxa:
#'   \describe{
#'     \item{"remainder"}{(default) Fold into the remainder group — their abundance
#'       is preserved but they are never shown as a named top taxon}
#'     \item{"remove"}{Remove entirely before top-N selection}
#'     \item{"group"}{Merge under a single Incertae Sedis label, kept as its own
#'       category separate from remainder}
#'   }
#' @param unclassified Character. How to handle unclassified taxa (NA, empty
#'   string, bare prefix like "g__", or "uncultured"):
#'   \describe{
#'     \item{"remainder"}{(default) Fold into the remainder group}
#'     \item{"remove"}{Remove entirely before aggregation}
#'     \item{"group"}{Merge under a single unclassified label, kept as its own
#'       category separate from remainder}
#'   }
#' @param remainder_label Character label for grouped non-top taxa (default: "(Remainder)")
#' @param incertae_sedis_label Character label for grouped Incertae Sedis taxa
#'   (default: "(Incertae Sedis Group)"). Only used when \code{incertae_sedis = "group"}.
#' @param unclassified_label Character label for grouped unclassified taxa
#'   (default: "(Unclassified)"). Only used when \code{unclassified = "group"}.
#' @return Aggregated TreeSummarizedExperiment object
aggregate_taxa <- function(tse_obj, rank, top = 10, incertae_sedis = "remainder",
                           unclassified = "remainder",
                           remainder_label = "(Remainder)",
                           incertae_sedis_label = "(Incertae Sedis Group)",
                           unclassified_label = "(Unclassified)") {
  # Validate parameters
  incertae_sedis <- match.arg(incertae_sedis, c("remainder", "remove", "group"))
  unclassified <- match.arg(unclassified, c("remainder", "remove", "group"))

  # Local helpers
  is_incertae_sedis <- function(taxa_names) {
    grepl("Incertae.Sedis", taxa_names, ignore.case = TRUE)
  }

  is_unclassified <- function(taxa_names) {
    is.na(taxa_names) | taxa_names == "" |
      grepl("^[a-z]__$", taxa_names) |
      grepl("^uncultured$", taxa_names, ignore.case = TRUE)
  }

  # Step 0: Handle unclassified values in rowData BEFORE agglomeration
  rank_values <- SummarizedExperiment::rowData(tse_obj)[[rank]]
  unclass_mask <- is_unclassified(rank_values)

  if (unclassified == "remove") {
    tse_obj <- tse_obj[!unclass_mask, ]
  } else if (unclassified %in% c("group", "remainder")) {
    # Relabel unclassified values so they survive agglomeration as a single group
    if (any(unclass_mask)) {
      rank_values[unclass_mask] <- unclassified_label
      SummarizedExperiment::rowData(tse_obj)[[rank]] <- rank_values
    }
  }

  # Step 1: Agglomerate by the specified rank (na.rm = FALSE since we already handled NAs)
  tse_obj <- mia::agglomerateByRank(tse_obj, rank = rank, na.rm = FALSE)

  # Step 1b: agglomerateByRank may append _N suffixes to deduplicate rows that
  # share the same rank value but differ at other ranks. Collapse any suffixed
  # unclassified labels back to the canonical label so they merge properly later.
  rank_values_post <- SummarizedExperiment::rowData(tse_obj)[[rank]]
  suffix_pattern <- paste0("^\\Q", unclassified_label, "\\E(_\\d+)?$")
  suffixed_mask <- grepl(suffix_pattern, rank_values_post, perl = TRUE)
  if (any(suffixed_mask)) {
    rank_values_post[suffixed_mask] <- unclassified_label
    SummarizedExperiment::rowData(tse_obj)[[rank]] <- rank_values_post
    tse_obj <- mia::agglomerateByVariable(tse_obj, by = "rows", f = rank)
  }

  # Step 2: Handle incertae sedis taxa based on incertae_sedis parameter
  row_names <- rownames(SummarizedExperiment::rowData(tse_obj))

  if (incertae_sedis == "remove") {
    # Remove incertae sedis taxa entirely
    incertae_mask <- !is_incertae_sedis(row_names)
    tse_obj <- tse_obj[incertae_mask, ]
  } else if (incertae_sedis == "group") {
    # Group incertae sedis rows BEFORE getTop.
    # Replace any row name containing Incertae Sedis with a unified label.
    renamed_incertae <- vapply(
      row_names,
      function(x) {
        if (is_incertae_sedis(x)) {
          incertae_sedis_label
        } else {
          x
        }
      },
      character(1)
    )
    SummarizedExperiment::rowData(tse_obj)[[rank]] <- renamed_incertae

    # Agglomerate by variable to merge rows with the same (now updated) taxon names
    tse_obj <- mia::agglomerateByVariable(tse_obj, by = "rows", f = rank)
  }
  # For incertae_sedis == "remainder": keep rows as-is, they'll be folded into
  # remainder during the top-N step below

  # Step 3: Handle top = NULL case (return all taxa without grouping)
  if (is.null(top)) {
    # Don't group any taxa as remainder - keep all original names
    # Remove the "Confidence" column from rowData as it doesn't make sense for aggregated taxa.
    rowdata <- SummarizedExperiment::rowData(tse_obj)
    SummarizedExperiment::rowData(tse_obj) <- rowdata[, !(colnames(rowdata) %in% "Confidence")]
    return(tse_obj)
  }

  # Step 4: Determine top taxa from the aggregated object.
  # Exclude incertae sedis and unclassified rows from top-N selection when they
  # should fold into remainder.
  current_names <- SummarizedExperiment::rowData(tse_obj)[[rank]]
  exclude_mask <- rep(FALSE, length(current_names))

  if (incertae_sedis == "remainder") {
    exclude_mask <- exclude_mask | is_incertae_sedis(current_names)
  }
  if (unclassified %in% c("remainder", "group")) {
    exclude_mask <- exclude_mask | (current_names == unclassified_label)
  }

  n_available <- sum(!exclude_mask)
  if (top > n_available) {
    stop(sprintf(
      "aggregate_taxa(rank = '%s'): requested top = %d but only %d taxa are available after excluding %d Incertae Sedis / unclassified row(s) (total %d at this rank). Lower `top` or pass `top = NULL` to keep all taxa.",
      rank, top, n_available, sum(exclude_mask), length(current_names)
    ), call. = FALSE)
  }
  # Top rows by mean relative abundance (the getTop() default), selected by position:
  # agglomeration can leave homonymous rows (same genus, different family) sharing a rowname,
  # which getTop() reports under de-duplicated names that match no row.
  mean_abund <- rowMeans(SummarizedExperiment::assay(tse_obj, "relabundance"))
  mean_abund[exclude_mask] <- -Inf
  is_top <- seq_along(mean_abund) %in% order(mean_abund, decreasing = TRUE)[seq_len(top)]

  # Step 5: Rename taxa (non-top rows become the remainder unless kept as their own group).
  rank_vals <- as.character(current_names)
  keep_label <- is_top |
    (incertae_sedis == "group" & rank_vals == incertae_sedis_label) |
    (unclassified == "group" & rank_vals == unclassified_label)
  renamed <- ifelse(keep_label, rank_vals, remainder_label)

  # Step 6: Convert to factor with desired ordering.
  unique_taxa <- unique(renamed)
  levels_order <- c(
    setdiff(unique_taxa, c(remainder_label, incertae_sedis_label, unclassified_label)),
    if (incertae_sedis == "group" && incertae_sedis_label %in% unique_taxa) incertae_sedis_label else NULL,
    remainder_label,
    if (unclassified == "group" && unclassified_label %in% unique_taxa) unclassified_label else NULL
  )
  SummarizedExperiment::rowData(tse_obj)[[rank]] <- factor(renamed, levels = levels_order)

  # Step 7: Re-aggregate by variable so that the new factor grouping is applied.
  tse_obj <- mia::agglomerateByVariable(tse_obj, by = "rows", f = rank)

  # Remove the "Confidence" column from rowData as it doesn't make sense for aggregated taxa.
  rowdata <- SummarizedExperiment::rowData(tse_obj)
  SummarizedExperiment::rowData(tse_obj) <- rowdata[, !(colnames(rowdata) %in% "Confidence")]

  tse_obj
}
