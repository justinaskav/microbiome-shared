#' Run dbRDA on a TSE object and return Type I and Type III ANOVA tables.
#'
#' Requires `tidy_vegan_anova()` (beta_diversity.R).
#'
#' @param tse SummarizedExperiment / TSE object with relabundance assay and colData.
#' @param rda_formula_str Character formula string, e.g. "~ Diagnosis + Site".
#' @param strata_var Optional character; column in colData to use as strata.
#' @param permutations Number of permutations (default 999).
#' @return List with elements: type1, type3 (`tidy_vegan_anova()` data.frames, numeric p-values;
#'   NULL if the test failed), n_samples.
run_dbrda_analysis <- function(tse, rda_formula_str, strata_var = NULL, permutations = 999) {
  dat <- dbrda_data(tse, rda_formula_str, strata_var)
  if (dat$n <= 5) return(list(type1 = NULL, type3 = NULL, n_samples = dat$n))

  bray_dist <- dat$dist
  fit <- vegan::dbrda(stats::as.formula(paste("bray_dist", rda_formula_str)),
                      data = dat$data, na.action = stats::na.fail)
  anova_by <- function(by) tryCatch(
    tidy_vegan_anova(stats::anova(fit, by = by, permutations = permutations, strata = dat$strata)),
    error = function(e) NULL)

  list(type1 = anova_by("terms"), type3 = anova_by("margin"), n_samples = dat$n)
}

#' Complete-case data and Bray-Curtis distance for a dbRDA formula
#'
#' @inheritParams run_dbrda_analysis
#' @return list(dist, data, strata, n): the Bray-Curtis distance of the `relabundance` assay
#'   over complete cases of the formula (and strata) variables, their colData, the strata
#'   vector (NULL if none) and the number of complete cases.
#' @keywords internal
dbrda_data <- function(tse, rda_formula_str, strata_var = NULL) {
  metadata <- as.data.frame(SummarizedExperiment::colData(tse))
  vars <- c(all.vars(stats::as.formula(rda_formula_str)), strata_var)
  complete_idx <- stats::complete.cases(metadata[, vars, drop = FALSE])
  rel <- SummarizedExperiment::assay(tse, "relabundance")[, complete_idx, drop = FALSE]
  data <- metadata[complete_idx, , drop = FALSE]
  list(
    dist = vegan::vegdist(t(rel), method = "bray"),
    data = data,
    strata = if (!is.null(strata_var)) data[[strata_var]] else NULL,
    n = sum(complete_idx)
  )
}

#' Format a dbRDA ANOVA data.frame as a gt table.
#'
#' @param anova_df Data.frame as returned by run_dbrda_analysis()$type1 or $type3.
#' @param title Table title string.
format_rda_table <- function(anova_df, title) {
  anova_df |>
    dplyr::select(term, df, sum_of_sqs, f_statistic, p_value) |>
    gt::gt() |>
    gt::tab_header(title = title) |>
    gt::cols_label(term = "Variable", df = "Df", sum_of_sqs = "Sum of Squares",
                   f_statistic = "F", p_value = "p-value") |>
    gt::fmt_number(columns = sum_of_sqs, decimals = 4) |>
    gt::fmt_number(columns = f_statistic, decimals = 3) |>
    gt::fmt(columns = p_value, fns = function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))) |>
    print_interactive_table()
}
