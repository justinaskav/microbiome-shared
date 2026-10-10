#' Calculate Taxonomic Abundance Summary Statistics
#'
#' Generates summary statistics for relative abundances of taxa by group.
#' Supports multiple summary types: mean+CI, median+IQR, median+CI, mean+SD.
#'
#' @param tse TreeSummarizedExperiment object (already aggregated at desired level)
#' @param taxa Character vector of taxa names to include. If NULL, uses all taxa.
#' @param group_var Column name in colData for grouping (e.g., "Type", "Group")
#' @param statistics Character vector of statistics to calculate. Options:
#'   \itemize{
#'     \item \code{"mean_ci"}: Mean with 95% confidence interval (normal approximation or bootstrap)
#'     \item \code{"median_iqr"}: Median with interquartile range (Q1-Q3)
#'     \item \code{"median_ci"}: Median with 95% confidence interval, calculated using the
#'       percentile bootstrap method — resamples with replacement \code{n_boot} times,
#'       computes the median of each resample, then takes the 2.5th and 97.5th percentiles
#'       of the bootstrap distribution as CI bounds
#'     \item \code{"mean_sd"}: Mean with standard deviation
#'   }
#'   Default: \code{"mean_ci"}. Multiple can be specified, e.g. \code{c("mean_ci", "median_ci")}.
#' @param assay_type Assay type to use (default: "relabundance")
#' @param ci_method Method for mean 95% CI: "normal" (default), "t" or "bootstrap"
#' @param n_boot Number of bootstrap replicates for bootstrap-based CIs (default: 1000).
#'   Used by \code{ci_method = "bootstrap"} for mean CI and always used for median CI.
#' @param as_percentage Logical, express values as percentages (default: TRUE)
#' @param digits Number of decimal places (default: 1)
#'
#' @return Tibble with `taxon`, `group` and the columns described in `describe_groups()`.
#'
#' @export
calc_taxonomic_summary <- function(
  tse,
  taxa = NULL,
  group_var,
  statistics = "mean_ci",
  assay_type = "relabundance",
  ci_method = "normal",
  n_boot = 1000,
  as_percentage = TRUE,
  digits = 1
) {
  # Validate inputs
  if (!inherits(tse, "TreeSummarizedExperiment") && !inherits(tse, "SummarizedExperiment")) {
    stop("tse must be a TreeSummarizedExperiment or SummarizedExperiment object")
  }

  if (!(assay_type %in% SummarizedExperiment::assayNames(tse))) {
    stop(paste("Assay", assay_type, "not found"))
  }

  if (!(group_var %in% colnames(SummarizedExperiment::colData(tse)))) {
    stop(paste("Group variable", group_var, "not found in colData"))
  }

  melted <- mia::meltSE(tse, assay.type = assay_type, add.col = group_var)
  if (!is.null(taxa) && length(taxa) > 0) {
    melted <- dplyr::filter(melted, FeatureID %in% taxa)
  }
  melted$abundance <- melted[[assay_type]] * if (as_percentage) 100 else 1

  describe_groups(
    melted, value_col = "abundance", group_cols = c("FeatureID", group_var),
    statistics = statistics, ci_method = ci_method, n_boot = n_boot, digits = digits,
    unit = if (as_percentage) "%" else "", nonneg = as_percentage
  ) |>
    dplyr::rename(taxon = FeatureID, group = dplyr::all_of(group_var)) |>
    dplyr::relocate(taxon, group)
}

#' Descriptive statistics per group, with formatted display strings
#'
#' Shared by `calc_taxonomic_summary()` and `calc_alpha_diversity_summary()`.
#'
#' @param df Long data frame.
#' @param value_col Column holding the values.
#' @param group_cols Grouping columns (kept in the output).
#' @param statistics Any of "mean_ci", "median_iqr", "median_ci", "mean_sd".
#' @param ci_method Mean CI: "t", "normal" or "bootstrap" (percentile).
#' @param n_boot Bootstrap replicates (median CI, and mean CI when `ci_method = "bootstrap"`).
#' @param digits Decimal places in the formatted strings.
#' @param unit Suffix after point estimates in the formatted strings (e.g. "%").
#' @param nonneg Clamp lower CI bounds at 0 (for non-negative quantities).
#' @return One row per group: `n`, the columns of each requested statistic, its formatted
#'   string (`formatted_mean`, `formatted_median`, `formatted_mean_sd`; with both median
#'   statistics, `formatted_median` is the CI version and `formatted_median_iqr` the IQR one),
#'   and `formatted`, an alias for the first of mean_ci > median_ci > median_iqr > mean_sd.
#' @keywords internal
describe_groups <- function(df, value_col, group_cols,
                            statistics = c("mean_ci", "median_ci"),
                            ci_method = "t", n_boot = 1000, digits = 2,
                            unit = "", nonneg = FALSE) {
  statistics <- match.arg(statistics, c("mean_ci", "median_iqr", "median_ci", "mean_sd"),
                          several.ok = TRUE)
  ci_method <- match.arg(ci_method, c("t", "normal", "bootstrap"))
  grouped <- dplyr::group_by(df, dplyr::across(dplyr::all_of(group_cols)))

  out <- dplyr::summarise(
    grouped,
    n = dplyr::n(),
    mean = mean(.data[[value_col]], na.rm = TRUE),
    sd = stats::sd(.data[[value_col]], na.rm = TRUE),
    median = stats::median(.data[[value_col]], na.rm = TRUE),
    q25 = stats::quantile(.data[[value_col]], 0.25, na.rm = TRUE, names = FALSE),
    q75 = stats::quantile(.data[[value_col]], 0.75, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  ) |>
    dplyr::mutate(se = sd / sqrt(n), iqr = q75 - q25)

  # Bootstrap passes run mean-then-median, group by group, so a given seed reproduces the
  # draws of the per-statistic implementation this replaced.
  if ("mean_ci" %in% statistics) {
    if (ci_method == "bootstrap") {
      ci <- dplyr::summarise(grouped, ci = list(unlist(calc_bootstrap_ci(.data[[value_col]], n_boot))),
                             .groups = "drop")$ci
      out$ci_lower <- vapply(ci, `[`, numeric(1), 1)
      out$ci_upper <- vapply(ci, `[`, numeric(1), 2)
    } else {
      crit <- if (ci_method == "t") stats::qt(0.975, df = out$n - 1) else 1.96
      out$ci_lower <- out$mean - crit * out$se
      out$ci_upper <- out$mean + crit * out$se
    }
    if (nonneg) out$ci_lower <- pmax(out$ci_lower, 0)
  }
  if ("median_ci" %in% statistics) {
    ci <- dplyr::summarise(grouped, ci = list(unlist(calc_bootstrap_median_ci(.data[[value_col]], n_boot))),
                           .groups = "drop")$ci
    out$median_ci_lower <- vapply(ci, `[`, numeric(1), 1)
    out$median_ci_upper <- vapply(ci, `[`, numeric(1), 2)
    if (nonneg) out$median_ci_lower <- pmax(out$median_ci_lower, 0)
  }

  num <- paste0("%.", digits, "f")
  est <- paste0(num, gsub("%", "%%", unit, fixed = TRUE))
  ci_fmt <- paste0(est, " (", num, "-", num, ")")
  keep <- c(group_cols, "n")
  if ("mean_ci" %in% statistics) {
    out$formatted_mean <- sprintf(ci_fmt, out$mean, out$ci_lower, out$ci_upper)
    keep <- c(keep, "mean", "sd", "se", "ci_lower", "ci_upper", "formatted_mean")
  }
  if ("median_iqr" %in% statistics) {
    iqr_col <- if ("median_ci" %in% statistics) "formatted_median_iqr" else "formatted_median"
    out[[iqr_col]] <- sprintf(paste0(est, " [", est, "]"), out$median, out$iqr)
    keep <- c(keep, "median", "q25", "q75", "iqr", iqr_col)
  }
  if ("median_ci" %in% statistics) {
    out$formatted_median <- sprintf(ci_fmt, out$median, out$median_ci_lower, out$median_ci_upper)
    keep <- c(keep, "median", "median_ci_lower", "median_ci_upper", "formatted_median")
  }
  if ("mean_sd" %in% statistics) {
    out$formatted_mean_sd <- sprintf(paste0(est, " \u00b1 ", est), out$mean, out$sd)
    keep <- c(keep, "mean", "sd", "formatted_mean_sd")
  }
  primary <- intersect(c("mean_ci", "median_ci", "median_iqr", "mean_sd"), statistics)[1]
  out$formatted <- out[[c(mean_ci = "formatted_mean", median_ci = "formatted_median",
                          median_iqr = "formatted_median", mean_sd = "formatted_mean_sd")[[primary]]]]
  out[, c(unique(keep), "formatted")]
}

#' Percentile bootstrap CI for the mean
#' @keywords internal
calc_bootstrap_ci <- function(x, n_boot = 1000, conf_level = 0.95) {
  ci <- .bootstrap_ci(x, mean, n_boot, conf_level)
  list(ci_lower = ci[1], ci_upper = ci[2])
}

#' Bootstrap CI helper for the median (percentile method)
#'
#' Computes a 95% confidence interval for the median using the percentile
#' bootstrap method: resamples the data with replacement \code{n_boot} times,
#' computes the median of each bootstrap sample, then takes the 2.5th and
#' 97.5th percentiles of the resulting bootstrap distribution as the CI bounds.
#'
#' @param x Numeric vector of observations
#' @param n_boot Number of bootstrap replicates (default: 1000)
#' @param conf_level Confidence level (default: 0.95)
#' @return Named list with \code{median_ci_lower} and \code{median_ci_upper}
#' @keywords internal
calc_bootstrap_median_ci <- function(x, n_boot = 1000, conf_level = 0.95) {
  ci <- .bootstrap_ci(x, stats::median, n_boot, conf_level)
  list(median_ci_lower = ci[1], median_ci_upper = ci[2])
}

.bootstrap_ci <- function(x, stat, n_boot, conf_level) {
  x <- stats::na.omit(x)
  if (length(x) < 2) return(c(NA_real_, NA_real_))
  alpha <- (1 - conf_level) / 2
  boots <- replicate(n_boot, stat(sample(x, replace = TRUE)))
  unname(stats::quantile(boots, probs = c(alpha, 1 - alpha)))
}
