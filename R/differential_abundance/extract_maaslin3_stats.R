#' Extract and Format MaAsLin3 Statistics for Display
#'
#' Extracts MaAsLin3 results for a specific coefficient and formats them
#' for integration into taxonomic summary tables.
#'
#' Requires `prepare_maaslin3_data()` (prepare_maaslin3_data.R).
#'
#' @param maaslin3_result MaAsLin3 result object, results data frame or all_results.tsv path
#'   (see `prepare_maaslin3_data()`)
#' @param coefficient_name String matching the coefficient column (e.g., "TypeStent", "AOrBB")
#' @param taxa Character vector of taxa names to filter. If NULL, returns all taxa.
#' @param result_type Either "abundance" or "prevalence" (default: "abundance")
#' @param digits_coef Decimal places for coefficient and SE (default: 2)
#' @param digits_pval Decimal places for p-values/q-values (default: 3)
#'
#' @return data.frame with columns: taxon, coef, stderr, pval_individual,
#'   qval_individual, formatted_coef, formatted_pval, formatted_qval
#'
#' @details
#' Built on `prepare_maaslin3_data()` (prepare_maaslin3_data.R), so it accepts the same inputs
#' and drops features whose model reported an error.
#'
#' @examples
#' \dontrun{
#' # Extract Type coefficient statistics
#' type_stats <- extract_maaslin3_stats(
#'   maaslin3_result = maaslin3_genus,
#'   coefficient_name = "TypeStent",
#'   taxa = c("Clostridium", "Dialister")
#' )
#' }
#'
#' @export
extract_maaslin3_stats <- function(
  maaslin3_result,
  coefficient_name,
  taxa = NULL,
  result_type = "abundance",
  digits_coef = 2,
  digits_pval = 3
) {
  result_type <- match.arg(result_type, c("abundance", "prevalence"))
  all_res <- prepare_maaslin3_data(maaslin3_result, model_type = result_type)
  res <- all_res[all_res$coefficient_name == coefficient_name, , drop = FALSE]
  if (nrow(res) == 0) {
    warning("Coefficient '", coefficient_name, "' not found in MaAsLin3 results. ",
            "Available coefficients: ", paste(unique(all_res$coefficient_name), collapse = ", "))
  }

  format_pval <- function(p) {
    threshold <- 10^(-digits_pval)
    ifelse(is.na(p), "",
           ifelse(p < threshold, paste0("< ", format(threshold, scientific = FALSE)),
                  format(round(p, digits_pval), nsmall = digits_pval)))
  }

  out <- data.frame(
    taxon = res$taxon,
    coef = res$coef,
    stderr = res$stderr,
    pval_individual = res$pval,
    qval_individual = res$qval,
    formatted_coef = sprintf(paste0("%.", digits_coef, "f \u00b1 %.", digits_coef, "f"), res$coef, res$stderr),
    formatted_pval = format_pval(res$pval),
    formatted_qval = format_pval(res$qval)
  )

  # With `taxa`, return one row per requested taxon (empty strings where not tested).
  if (!is.null(taxa) && length(taxa) > 0) {
    out <- dplyr::left_join(data.frame(taxon = taxa), out, by = "taxon")
    out$formatted_coef[is.na(out$coef)] <- ""
    out$formatted_pval[is.na(out$pval_individual)] <- ""
    out$formatted_qval[is.na(out$qval_individual)] <- ""
  }
  out
}
