#' Two-stage variable selection for dbRDA: VIF pruning then forward R² selection.
#'
#' Requires `dbrda_data()` (run_dbrda_analysis.R) and `tidy_vegan_anova()` (beta_diversity.R).
#'
#' Stage A: iteratively drop the highest-VIF variable (vegan::vif.cca) until all
#' remaining non-anchor VIFs are <= vif_threshold.
#' Stage B: forward selection on adjusted R² with Blanchet et al. double stopping
#' criterion (vegan::ordiR2step), anchoring the scientific variables.
#'
#' @param tse SummarizedExperiment / TSE with relabundance assay and colData.
#' @param rda_formula_str Character formula string for the saturated model, e.g.
#'   "~ Diagnosis + Site + Age_z + BMI".
#' @param anchor_vars Character vector of variables protected from both stages.
#' @param strata_var Optional column in colData used as permutation strata.
#' @param vif_threshold Maximum VIF allowed after Stage A (default 10).
#' @param permutations Number of permutations for ANOVA (default 999).
#' @return Named list with: vif_log, dropped_by_vif, kept_after_vif, selected_vars,
#'   step_result, type1_reduced, type3_reduced, n_samples, adj_r2_full,
#'   adj_r2_reduced, reduced_formula.
dbrda_variable_selection <- function(tse, rda_formula_str,
                                     anchor_vars = c("Diagnosis", "Site"),
                                     strata_var = NULL,
                                     vif_threshold = 10,
                                     permutations = 999) {
  dat <- dbrda_data(tse, rda_formula_str, strata_var)
  if (dat$n <= 5) return(NULL)
  formula_vars <- all.vars(stats::as.formula(rda_formula_str))
  bray_dist <- dat$dist
  metadata_subset <- dat$data
  strata_val <- dat$strata
  perm_ctrl <- if (!is.null(strata_val)) {
    permute::how(nperm = permutations, blocks = strata_val)
  } else permutations

  # Max VIF per term. Model-matrix columns are mapped to their term through the `assign`
  # attribute, so a variable whose name prefixes another's (Age / Age_z) cannot absorb its VIFs.
  agg_term_vif <- function(vif_vec, terms) {
    mm <- stats::model.matrix(stats::as.formula(paste("~", paste(terms, collapse = " + "))),
                              data = metadata_subset)
    col_term <- stats::setNames(c("(Intercept)", terms)[attr(mm, "assign") + 1], colnames(mm))
    vapply(terms, function(v) {
      x <- vif_vec[names(vif_vec) %in% names(col_term)[col_term == v]]
      if (length(x) == 0 || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    }, numeric(1))
  }

  # Stage A: VIF pruning
  current_vars <- formula_vars
  vif_log <- list()
  dropped <- character()
  iter <- 0
  repeat {
    iter <- iter + 1
    fmla <- stats::as.formula(paste("bray_dist ~", paste(current_vars, collapse = " + ")))
    dbrda_iter <- vegan::dbrda(fmla, data = metadata_subset, na.action = stats::na.fail)
    vifs <- vegan::vif.cca(dbrda_iter)
    term_vif <- agg_term_vif(vifs, current_vars)
    vif_log[[iter]] <- data.frame(
      iteration = iter,
      variable = names(term_vif),
      vif = round(as.numeric(term_vif), 2),
      stringsAsFactors = FALSE
    )
    droppable <- term_vif[!names(term_vif) %in% anchor_vars]
    if (length(droppable) == 0 || all(is.na(droppable)) ||
        max(droppable, na.rm = TRUE) <= vif_threshold) break
    drop_var <- names(which.max(droppable))
    dropped <- c(dropped, drop_var)
    current_vars <- setdiff(current_vars, drop_var)
    if (length(current_vars) <= length(anchor_vars)) break
  }
  vif_table <- do.call(rbind, vif_log)
  kept_after_vif <- current_vars

  # Stage B: forward R² selection anchored at scientific variables
  effective_anchors <- intersect(anchor_vars, kept_after_vif)
  step_result <- NULL
  selected_vars <- kept_after_vif
  if (length(kept_after_vif) > length(effective_anchors) && length(effective_anchors) > 0) {
    lower_fmla <- stats::as.formula(paste("bray_dist ~", paste(effective_anchors, collapse = " + ")))
    upper_fmla <- stats::as.formula(paste("bray_dist ~", paste(kept_after_vif, collapse = " + ")))
    dbrda_lower <- vegan::dbrda(lower_fmla, data = metadata_subset, na.action = stats::na.fail)
    dbrda_upper <- vegan::dbrda(upper_fmla, data = metadata_subset, na.action = stats::na.fail)
    step_result <- tryCatch(
      vegan::ordiR2step(dbrda_lower, scope = stats::formula(dbrda_upper),
                        direction = "forward", permutations = perm_ctrl,
                        R2scope = TRUE, trace = FALSE),
      error = function(e) NULL
    )
    if (!is.null(step_result)) {
      selected_vars <- all.vars(stats::formula(step_result))
      selected_vars <- selected_vars[selected_vars %in% kept_after_vif]
    }
  }

  reduced_fmla <- stats::as.formula(paste("bray_dist ~", paste(selected_vars, collapse = " + ")))
  dbrda_reduced <- vegan::dbrda(reduced_fmla, data = metadata_subset, na.action = stats::na.fail)
  dbrda_full <- vegan::dbrda(stats::as.formula(paste("bray_dist", rda_formula_str)),
                              data = metadata_subset, na.action = stats::na.fail)

  anova_by <- function(by) tryCatch(
    tidy_vegan_anova(stats::anova(dbrda_reduced, by = by, permutations = permutations,
                                  strata = strata_val)),
    error = function(e) NULL)
  type1 <- anova_by("terms")
  type3 <- anova_by("margin")

  list(
    vif_log = vif_table,
    dropped_by_vif = dropped,
    kept_after_vif = kept_after_vif,
    selected_vars = selected_vars,
    step_result = step_result,
    type1_reduced = type1,
    type3_reduced = type3,
    n_samples = dat$n,
    adj_r2_full = vegan::RsquareAdj(dbrda_full)$adj.r.squared,
    adj_r2_reduced = vegan::RsquareAdj(dbrda_reduced)$adj.r.squared,
    reduced_formula = paste("~", paste(selected_vars, collapse = " + "))
  )
}

#' Render a per-iteration VIF log as a wide gt table (variables × iterations).
#'
#' @param vif_log Data.frame with columns: iteration, variable, vif.
#' @param title Table title string.
format_vif_log_table <- function(vif_log, title) {
  if (is.null(vif_log) || nrow(vif_log) == 0) return(invisible(NULL))
  wide <- tidyr::pivot_wider(vif_log, names_from = iteration,
                              values_from = vif, names_prefix = "iter_")
  wide |>
    gt::gt() |>
    gt::tab_header(title = title) |>
    gt::sub_missing(missing_text = "—") |>
    gt::fmt_number(columns = dplyr::starts_with("iter_"), decimals = 2) |>
    print_interactive_table()
}

#' Print the textual summary of a dbRDA variable-selection result.
#'
#' Emits sample count, saturated/reduced adjusted R², VIF-dropped variables, the
#' selected formula, and the per-iteration VIF log table. No-op if `sel` is NULL.
#'
#' @param sel Result of `dbrda_variable_selection()` (or NULL).
#' @param vif_title Title for the VIF-pruning gt table.
print_dbrda_selection_summary <- function(sel, vif_title) {
  if (is.null(sel)) return(invisible(NULL))
  # "observations", not "samples": callers pass a patient-level object as often as a sample-level one.
  cat(sprintf("Complete observations: n = %d\n", sel$n_samples))
  cat(sprintf("Adjusted R² — saturated: %.4f, reduced: %.4f\n",
              sel$adj_r2_full, sel$adj_r2_reduced))
  # Report the max VIF actually observed, so a Stage A that pruned nothing is legible as a no-op
  # rather than as pruning the reader has to take on trust.
  max_vif <- if (!is.null(sel$vif_log) && nrow(sel$vif_log) > 0) {
    last_iter <- sel$vif_log[sel$vif_log$iteration == max(sel$vif_log$iteration), ]
    max(last_iter$vif, na.rm = TRUE)
  } else NA_real_
  if (length(sel$dropped_by_vif) > 0) {
    cat("Dropped by VIF (in order): ",
        paste(sel$dropped_by_vif, collapse = ", "), "\n", sep = "")
    cat(sprintf("Max VIF after pruning: %.2f\n", max_vif))
  } else {
    cat(sprintf("No variables exceeded the VIF threshold — max VIF %.2f, nothing pruned.\n",
                max_vif))
  }
  cat("Selected formula: ", sel$reduced_formula, "\n", sep = "")
  format_vif_log_table(sel$vif_log, vif_title)
}
