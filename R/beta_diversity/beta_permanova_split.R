#' Patient-aware PERMANOVA (between- vs within-subject terms)
#'
#' Repeated biopsies per patient make a naive sample-level PERMANOVA pseudoreplicated for
#' *between-subject* factors (Diagnosis, Sex, Age, clinical covariates). Restricting
#' permutations within patient — whether written as `vegan::adonis2(strata = Patient)` or as
#' `permute::how(blocks = Patient)`, which are equivalent — does not rescue them: such a factor
#' takes one value per patient, so it is *identical in every permutation*. Its permutation
#' distribution is therefore degenerate and its p-value carries no information (the observed
#' statistic is compared against a null built by reshuffling something that never changed).
#' Leaving them unrestricted instead treats each biopsy as an independent patient, which is
#' anticonservative. Neither option is a valid test, so this helper splits the inference:
#'
#'   * **Between-subject terms** are tested on a **patient-level** object
#'     (`aggregate_to_patient()`: mean of per-sample relative abundances, one observation per
#'     patient; free permutation).
#'   * **Within-subject terms** (Site) are tested at the **sample level** with permutations
#'     restricted *within* patient (`permute::how(blocks = Patient)`), adjusting for the
#'     between-subject terms.
#'
#' Between-subject distances are recomputed on the aggregated patient profiles (averaging stored
#' sample-level distances would not give a valid patient-level distance); the stored sample-level
#' distance in `metadata` is reused for the within-subject test.
#'
#' Requires `aggregate_to_patient()` / `check_between_subject()` (aggregate_to_patient.R) and
#' `tidy_vegan_anova()` (beta_diversity.R).
#'
#' @param tse TreeSummarizedExperiment with a `counts` assay, a `rowTree` (for UniFrac), and
#'   sample-level distances in `metadata` (keys `jaccard`, `bray`, `unifrac`, `unifrac_weighted`).
#' @param between_formula RHS string of between-subject terms, e.g. "Diagnosis + Sex + Age_z".
#' @param within_formula RHS string of within-subject terms tested at sample level, e.g. "Site";
#'   `NULL` to skip.
#' @param metric_filter Subset of c("jaccard", "bray", "unifrac_unweighted", "unifrac_weighted").
#' @param patient_var colData column identifying the subject (default "Patient").
#' @param permutations Number of permutations (default 9999).
#' @param by adonis2 term-testing mode, `"margin"` (default; each term adjusted for all others)
#'   or `"terms"` (sequential). Note the sample-level `analyze_beta_diversity()` uses `"terms"`,
#'   so R² partitioning is not directly comparable between the two outputs unless this matches.
#' @param seed RNG seed for reproducible permutation p-values.
#'
#' @return list with `between` and `within`, each shaped like `analyze_beta_diversity()` output
#'   (`$analysis[[metric_key]]$permanova` data.frame with columns term, df, sum_of_sqs, r2,
#'   f_statistic, p_value, q_value), plus `n_patients` (those in the between test) and
#'   `n_samples` (those in the within test).
#' @export
beta_permanova_split <- function(tse,
                                 between_formula,
                                 within_formula = NULL,
                                 metric_filter = c("bray", "unifrac_unweighted"),
                                 patient_var = "Patient",
                                 permutations = 9999,
                                 by = c("margin", "terms"),
                                 seed = 2026) {
  stopifnot(inherits(tse, "SummarizedExperiment"))
  by <- match.arg(by)
  # metric key in `metadata` (sample-level/within) and the UniFrac weighting flag (between)
  meta_key  <- c(jaccard = "jaccard", bray = "bray",
                 unifrac_unweighted = "unifrac", unifrac_weighted = "unifrac_weighted")
  unknown <- setdiff(metric_filter, names(meta_key))
  if (length(unknown) > 0) {
    stop("unknown metric(s) in metric_filter: ", paste(unknown, collapse = ", "),
         ". Valid: ", paste(names(meta_key), collapse = ", "))
  }
  if (length(metric_filter) == 0) stop("metric_filter is empty; nothing to test.")

  # distance from a patient-level relabundance object, per metric
  .patient_dist <- function(ptse, rel, mk) {
    switch(mk,
      jaccard            = vegan::vegdist(t(rel), method = "jaccard", binary = TRUE),
      bray               = vegan::vegdist(t(rel), method = "bray"),
      unifrac_unweighted = mia::getDissimilarity(ptse, method = "unifrac",
                                                 assay.type = "relabundance", weighted = FALSE),
      unifrac_weighted   = mia::getDissimilarity(ptse, method = "unifrac",
                                                 assay.type = "relabundance", weighted = TRUE))
  }

  cd <- as.data.frame(SummarizedExperiment::colData(tse))
  between_vars <- all.vars(stats::as.formula(paste("~", between_formula)))
  # Validate: warn if any "between" term actually varies within patient.
  check_between_subject(tse, between_vars, patient_var = patient_var, warn = TRUE)

  # ---- patient-level aggregation (between-subject test) ----
  ptse <- aggregate_to_patient(tse, patient_var = patient_var)
  # Patients missing a between-subject covariate are dropped (adonis2 would abort on NA).
  complete <- stats::complete.cases(as.data.frame(SummarizedExperiment::colData(ptse))[, between_vars, drop = FALSE])
  if (any(!complete)) {
    message(sum(!complete), " patient(s) with missing between-subject covariates dropped from the between test.")
    ptse <- ptse[, complete]
  }
  rel  <- as.matrix(SummarizedExperiment::assay(ptse, "relabundance"))
  md   <- as.data.frame(SummarizedExperiment::colData(ptse))[, between_vars, drop = FALSE]
  between_analysis <- list()
  for (mk in metric_filter) {
    dmat <- .patient_dist(ptse, rel, mk)
    set.seed(seed)
    ad <- vegan::adonis2(stats::as.formula(paste("dmat ~", between_formula)),
                         data = md, permutations = permutations, by = by)
    between_analysis[[mk]] <- list(permanova = tidy_vegan_anova(ad))
  }

  # ---- sample-level within-subject test (within-patient permutation) ----
  within_analysis <- NULL
  n_samples <- nrow(cd)
  if (!is.null(within_formula)) {
    within_terms <- all.vars(stats::as.formula(paste("~", within_formula)))
    full_rhs <- paste(within_formula, "+", between_formula)
    model_vars <- unique(c(within_terms, between_vars, patient_var))
    smd <- cd[, model_vars, drop = FALSE]
    ok <- stats::complete.cases(smd)
    smd <- smd[ok, , drop = FALSE]
    n_samples <- nrow(smd)  # the n the within test actually ran on, not the unfiltered column count
    within_analysis <- list()
    md_all <- S4Vectors::metadata(tse)
    for (mk in metric_filter) {
      dist_key <- meta_key[[mk]]
      if (is.null(md_all[[dist_key]])) {
        stop("metadata(tse)[['", dist_key, "']] not found — the sample-level distance for '",
             mk, "' is required for the within-subject test. Available: ",
             paste(names(md_all), collapse = ", "))
      }
      d <- as.matrix(md_all[[dist_key]])
      keep <- rownames(smd)
      missing_rows <- setdiff(keep, rownames(d))
      if (length(missing_rows) > 0) {
        stop("distance '", dist_key, "' is missing ", length(missing_rows),
             " of the samples being tested (e.g. ", missing_rows[1], ").")
      }
      d <- stats::as.dist(d[keep, keep])
      # Permutations restricted within patient: Site varies inside a patient, so shuffling
      # labels within each patient breaks its association with the response while holding the
      # patient effect and every between-subject covariate fixed.
      ctrl <- permute::how(blocks = factor(smd[[patient_var]]), nperm = permutations)
      set.seed(seed)
      ad <- vegan::adonis2(stats::as.formula(paste("d ~", full_rhs)),
                           data = smd, permutations = ctrl, by = by)
      within_analysis[[mk]] <- list(permanova = tidy_vegan_anova(ad, keep_terms = within_terms))
    }
  }

  list(
    between = list(analysis = between_analysis),
    within = if (!is.null(within_analysis)) list(analysis = within_analysis) else NULL,
    n_patients = ncol(ptse),
    n_samples = n_samples
  )
}
