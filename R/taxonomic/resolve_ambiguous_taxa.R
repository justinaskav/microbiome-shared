#' Resolve ambiguous rank labels against their parent lineage
#'
#' Reference databases (SILVA, UNITE) reuse a handful of non-informative tokens
#' — most notably \code{Incertae_Sedis} — as the rank value for many unrelated
#' lineages. Because \code{\link{aggregate_taxa}} (and \code{mia::agglomerateByRank}
#' underneath it) key on the rank string alone, every such lineage collapses into a
#' single feature. In environmental 16S data that one bin can hold a third of all reads,
#' drawn from Actinomycetota, Chloroflexota, Pseudomonadota, Acidobacteriota and
#' others, which is neither interpretable nor testable.
#'
#' This function rewrites those values in place as
#' \code{"<nearest informative ancestor> <original token>"} (e.g.
#' \code{"Gaiellales Incertae_Sedis"}), so the subsequent agglomeration keeps them
#' apart as distinct, lineage-labelled features. Call it \emph{before}
#' \code{aggregate_taxa()}.
#'
#' Rows whose ancestors are all uninformative are left untouched.
#'
#' @param tse_obj TreeSummarizedExperiment (or SummarizedExperiment) object
#' @param rank Character vector of rank columns to resolve (e.g. \code{"genus"},
#'   or \code{c("family", "genus")}). Each is resolved independently.
#' @param pattern Regular expression identifying ambiguous rank values
#'   (default matches \code{Incertae Sedis} / \code{Incertae_Sedis}, case-insensitive)
#' @param ranks Ordered rank columns from broadest to narrowest. Defaults to the
#'   intersection of the usual taxonomy columns with \code{rowData()}.
#' @return The object with the targeted rowData rank columns rewritten.
resolve_ambiguous_taxa <- function(tse_obj, rank = "genus",
                                   pattern = "Incertae.Sedis",
                                   ranks = NULL) {
  rd <- SummarizedExperiment::rowData(tse_obj)

  if (is.null(ranks)) {
    ranks <- intersect(
      c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
      colnames(rd)
    )
  }

  # An ancestor is only usable as a disambiguating prefix if it is itself informative:
  # not missing, not a bare rank prefix, not "uncultured", and not ambiguous either.
  is_informative <- function(x) {
    !is.na(x) & x != "" &
      !grepl("^[a-z]__$", x) &
      !grepl("^uncultured$", x, ignore.case = TRUE) &
      !grepl(pattern, x, ignore.case = TRUE)
  }

  for (r in rank) {
    if (!r %in% ranks) {
      stop(sprintf("resolve_ambiguous_taxa(): rank '%s' is not among the taxonomy columns (%s).",
                   r, paste(ranks, collapse = ", ")), call. = FALSE)
    }

    values <- as.character(rd[[r]])
    ambiguous <- !is.na(values) & grepl(pattern, values, ignore.case = TRUE)
    if (!any(ambiguous)) next

    # Walk up from the parent rank to the root, taking the first informative ancestor.
    parents <- rev(ranks[seq_len(match(r, ranks) - 1L)])
    prefix <- rep(NA_character_, length(values))
    for (p in parents) {
      needed <- ambiguous & is.na(prefix)
      if (!any(needed)) break
      candidate <- as.character(rd[[p]])
      usable <- needed & is_informative(candidate)
      prefix[usable] <- candidate[usable]
    }

    # Only prefix labels that are actually shared across lineages. SILVA reuses the bare
    # token `Incertae_Sedis` for hundreds of parents and must be split; UNITE already
    # qualifies its own (`Agaricostilbales_gen_Incertae_sedis`), so prefixing would only
    # produce `Agaricostilbales Agaricostilbales_gen_Incertae_sedis`.
    resolved <- ambiguous & !is.na(prefix)
    n_parents <- tapply(prefix[resolved], values[resolved], function(x) length(unique(x)))
    shared <- resolved & values %in% names(n_parents)[n_parents > 1L]
    values[shared] <- paste(prefix[shared], values[shared])
    rd[[r]] <- values
  }

  SummarizedExperiment::rowData(tse_obj) <- rd
  tse_obj
}
