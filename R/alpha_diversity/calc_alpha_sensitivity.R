#' Sensitivity Analyses for an Alpha Diversity LMM Term
#'
#' Refits an alpha-diversity LMM under alternatives that relax the assumptions the
#' residual diagnostics test, so that a reported effect can be shown not to depend on
#' them. Intended for terms whose effects are reported.
#'
#' The headline arm is `permutation`, which makes no distributional assumption at all.
#' It is only valid for a WITHIN-subject term (one that varies across a subject's
#' samples, e.g. Site): labels are shuffled within each subject, which holds the
#' subject random effect and every between-subject covariate fixed while breaking only
#' the term's association with the response. A between-subject term (Diagnosis, Sex,
#' Age) is constant within a subject and cannot be permuted this way -- the same
#' restriction that governs the PERMANOVA split -- so `permute_within` must name the
#' subject variable and `term` must vary within it. The function checks this and errors
#' rather than returning an invalid p-value.
#'
#' @param data Data frame with the alpha index and metadata.
#' @param y_var Character, name of the alpha diversity column.
#' @param term Character, the model term to test (must vary within `permute_within`).
#' @param lmm_formula_str Full LMM formula string, e.g.
#'   "Shannon ~ Diagnosis + Site + Sex + Age_z + (1|Patient)".
#' @param permute_within Character, the subject/grouping variable to permute inside
#'   (default "Patient").
#' @param n_perm Number of permutations (default 999).
#' @param exclude_var Optional name of a logical column in `data`; rows where it is
#'   TRUE are dropped in an extra sensitivity arm (e.g. dominated samples). NULL
#'   skips that arm.
#' @param exclude_label Label describing the `exclude_var` arm.
#' @param seed RNG seed for the permutation arm.
#'
#' Requires `classify_significance()` (R/util/classify_significance.R).
#'
#' @return Data frame with one row per arm: `arm`, `statistic`, `p_value`, `n`, `note`.
calc_alpha_sensitivity <- function(
  data,
  y_var,
  term,
  lmm_formula_str,
  permute_within = "Patient",
  n_perm = 999,
  exclude_var = NULL,
  exclude_label = "subset",
  seed = 2026
) {
  fml <- stats::as.formula(lmm_formula_str)
  prep <- function(d) {
    d |>
      dplyr::select(dplyr::all_of(unique(c(all.vars(fml), permute_within)))) |>
      stats::na.omit() |>
      dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor))
  }
  model_data <- prep(data)

  # A within-subject term has >1 distinct level inside at least some subjects.
  varies_within <- model_data |>
    dplyr::group_by(.data[[permute_within]]) |>
    dplyr::summarise(k = dplyr::n_distinct(.data[[term]]), .groups = "drop") |>
    dplyr::pull(k) |>
    max()
  if (varies_within < 2) {
    stop(sprintf(
      paste("term '%s' is constant within '%s', so it is a between-subject term and",
            "cannot be tested by within-subject permutation. Test it at the subject",
            "level instead."),
      term, permute_within))
  }

  # One fit per data set: returns c(F, p) for `term`.
  get_fp <- function(d) {
    suppressMessages(suppressWarnings(tryCatch(
      unlist(stats::anova(lmerTest::lmer(fml, data = d))[term, c("F value", "Pr(>F)")]),
      error = function(e) c(NA_real_, NA_real_))))
  }

  arms <- list()

  obs <- get_fp(model_data)
  f_obs <- obs[[1]]
  arms[["LMM (as reported)"]] <- list(
    stat = f_obs, p = obs[[2]], n = nrow(model_data),
    note = "Satterthwaite F-test")

  # Assumption-free arm.
  # NB: `.x[sample.int(length(.x))]`, not `sample(.x)` — for a length-1 numeric group
  # `sample(5)` returns a permutation of 1:5 rather than the value itself, which would
  # blow up the mutate. Indexing is unambiguous for every type and group size.
  set.seed(seed)
  f_perm <- replicate(n_perm, {
    dp <- model_data |>
      dplyr::group_by(.data[[permute_within]]) |>
      dplyr::mutate(dplyr::across(dplyr::all_of(term), ~ .x[sample.int(length(.x))])) |>
      dplyr::ungroup()
    get_fp(dp)[[1]]
  })
  n_ok <- sum(!is.na(f_perm))
  arms[[paste0("Within-", tolower(permute_within), " permutation")]] <- list(
    stat = f_obs,
    p = (1 + sum(f_perm >= f_obs, na.rm = TRUE)) / (1 + n_ok),
    n = nrow(model_data),
    note = sprintf("%d permutations, no distributional assumption", n_ok))

  # Rank-based arm: inverse-normal transform removes the tail entirely.
  d_int <- model_data
  d_int[[y_var]] <- qnorm((rank(d_int[[y_var]]) - 0.5) / nrow(d_int))
  fp_int <- get_fp(d_int)
  arms[["LMM on inverse-normal ranks"]] <- list(
    stat = fp_int[[1]], p = fp_int[[2]], n = nrow(d_int),
    note = "response rank-transformed to exact normality")

  if (!is.null(exclude_var)) {
    # Filter on the caller's frame, then prep, so rows stay aligned with the flag.
    d_sub <- prep(data[!data[[exclude_var]] %in% TRUE, , drop = FALSE])
    fp_sub <- get_fp(d_sub)
    arms[[exclude_label]] <- list(
      stat = fp_sub[[1]], p = fp_sub[[2]], n = nrow(d_sub),
      note = "sensitivity to influential samples; lower n reduces power")
  }

  dplyr::bind_rows(lapply(names(arms), function(nm) {
    data.frame(arm = nm, statistic = arms[[nm]]$stat, p_value = arms[[nm]]$p,
               n = arms[[nm]]$n, note = arms[[nm]]$note)
  })) |>
    dplyr::mutate(term = term, significance = classify_significance(p_value)) |>
    dplyr::relocate(term)
}
