#' Features that exist only to carry abundance into the TSS denominator
#'
#' Placeholders created by `resolve_ambiguous_taxa()` + `aggregate_taxa(unclassified =
#' "group", incertae_sedis = "remainder")`; they are fitted but never reported.
#' @export
denominator_only <- c("(Unclassified)", "(Incertae Sedis Group)", "(Remainder)", "(Others)")

#' Drop denominator-only rows from a MaAsLin3 result's result tables
#'
#' The q-values are deliberately left as MaAsLin3 computed them: `qval_joint` is a BH fit
#' pooled across the abundance and prevalence models, so re-running BH per table would not
#' reproduce it. Leaving them is conservative (one feature out of several hundred in the
#' FDR pool).
#'
#' @param res Result of `maaslin3::maaslin3()`.
#' @param labels Feature names to drop.
#' @return `res` with the rows removed from both result tables (fits are untouched).
#' @export
drop_denominator_features <- function(res, labels = denominator_only) {
  strip <- function(fit) {
    if (is.null(fit) || is.null(fit$results)) return(fit)
    fit$results <- fit$results[!fit$results$feature %in% labels, , drop = FALSE]
    fit
  }
  res$fit_data_abundance <- strip(res$fit_data_abundance)
  res$fit_data_prevalence <- strip(res$fit_data_prevalence)
  res
}
