#' Requires `classify_significance()` (R/util/classify_significance.R) and `describe_groups()`
#' (R/taxonomic/calc_taxonomic_summary.R).

#' Calculate Alpha Diversity Summary Statistics with LMM
#'
#' Generates descriptive statistics and model-based inference for alpha diversity.
#' Returns descriptive stats, model ANOVA, post-hoc comparisons and assumption
#' diagnostics. `plot_alpha_diversity_summary_lmm()` draws its significance brackets from
#' this function's `posthoc`, so figures and tables share one inference path.
#'
#' @param data Data frame containing alpha diversity values and metadata
#' @param y_var Character, name of the alpha diversity column
#' @param group_var Character, main grouping variable for comparisons
#' @param lmm_formula_str LMM formula string (e.g., "Shannon ~ Group + (1|Patient)")
#' @param ref_group Reference group (optional). With `test = "kruskal"` the post-hoc table is
#'   restricted to comparisons involving it; with `test = "anova"` all pairwise comparisons
#'   are returned (the reference only switches the post-hoc on, see `posthoc`).
#' @param test Inference method: "anova" (default) fits a linear (mixed) model and
#'   reports an ANOVA F-test with emmeans post-hoc comparisons; "kruskal" reports a
#'   non-parametric Kruskal-Wallis omnibus test with Dunn's test post-hoc comparisons.
#'   The "kruskal" path ignores any random-effect term in \code{lmm_formula_str}.
#' @param post_hoc_method P-value adjustment method (default: "BH"), applied across all
#'   pairwise comparisons.
#' @param ci_method Method for the descriptive mean CI: "t" (default), "normal" or "bootstrap"
#' @param statistics Character vector of descriptive statistics (see `describe_groups()`):
#'   "mean_ci", "median_iqr", "median_ci" (percentile bootstrap, `n_boot` replicates),
#'   "mean_sd". Default: \code{c("mean_ci", "median_ci")}.
#' @param n_boot Number of bootstrap replicates for median CI (default: 1000)
#' @param digits Number of decimal places (default: 2)
#' @param posthoc Logical, compute post-hoc comparisons (default: when `ref_group` is given).
#'
#' @return List with four elements:
#'   - descriptive: n, mean, sd, se, ci_lower, ci_upper, median, ..., formatted per group
#'   - model_anova: ANOVA table with F, p-value, effect sizes
#'   - posthoc: Pairwise comparisons with adjusted p-values
#'   - diagnostics: List of three data.frames checking the model's error assumptions --
#'     `levene` (homogeneity of residual variance across `group_var`), `normality`
#'     (Shapiro-Wilk on residuals) and `residuals` (fitted values and residuals for
#'     diagnostic plots). Both tests are run on model residuals rather than raw `y_var`
#'     values, since the assumption concerns the error term, not the marginal spread
#'     within each group. For mixed models these are approximate: `residuals()` returns
#'     conditional residuals (computed against shrunken BLUPs), and within-patient
#'     residuals are correlated, which Levene's internal ANOVA does not account for.
#'     Treat them as diagnostics, not as formal tests.
#'
#' @export
calc_alpha_diversity_summary <- function(
  data,
  y_var,
  group_var,
  lmm_formula_str,
  ref_group = NULL,
  test = "anova",
  post_hoc_method = "BH",
  ci_method = "t",
  statistics = c("mean_ci", "median_ci"),
  n_boot = 1000,
  digits = 2,
  posthoc = !is.null(ref_group)
) {
  # --- Prepare Data ---
  lmm_vars <- all.vars(stats::as.formula(lmm_formula_str))
  model_data <- data |>
    dplyr::select(dplyr::all_of(unique(c(lmm_vars, group_var)))) |>
    stats::na.omit() |>
    dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor))

  if (nrow(model_data) < 10) {
    warning("Insufficient complete cases for analysis")
    return(NULL)
  }

  # --- Descriptive Statistics ---
  descriptive <- describe_groups(
    model_data, value_col = y_var, group_cols = group_var, statistics = statistics,
    ci_method = ci_method, n_boot = n_boot, digits = digits
  ) |>
    dplyr::rename(group = dplyr::all_of(group_var)) |>
    dplyr::relocate(group, .after = dplyr::last_col())

  # --- Non-parametric path: Kruskal-Wallis omnibus + Dunn's post-hoc ---
  if (identical(test, "kruskal")) {
    kw_formula <- stats::as.formula(paste(y_var, "~", group_var))
    kw <- stats::kruskal.test(kw_formula, data = model_data)

    eps_sq <- tryCatch(
      rstatix::kruskal_effsize(model_data, kw_formula)$effsize[1],
      error = function(e) NA_real_
    )

    model_anova <- data.frame(
      term = group_var,
      statistic = unname(kw$statistic),
      df = unname(kw$parameter),
      p_value = kw$p.value,
      epsilon_sq = eps_sq
    ) |>
      dplyr::mutate(significance = classify_significance(p_value))

    # Dunn's test post-hoc (all pairwise), with the requested p-value adjustment
    posthoc_df <- if (posthoc) tryCatch({
      dunn <- rstatix::dunn_test(model_data, kw_formula, p.adjust.method = post_hoc_method)
      res <- data.frame(
        comparison = paste(dunn$group1, "-", dunn$group2),
        statistic = dunn$statistic,
        p_adj = dunn$p.adj
      ) |>
        dplyr::mutate(significance = classify_significance(p_adj))
      if (!is.null(ref_group)) res <- res[contrast_involves(res$comparison, ref_group), , drop = FALSE]
      res
    }, error = function(e) {
      warning(paste("Dunn's test failed:", e$message))
      NULL
    })

    return(list(
      descriptive = descriptive,
      model_anova = model_anova,
      posthoc = posthoc_df
    ))
  }

  # --- Fit Model ---
  is_mixed <- grepl("\\|", lmm_formula_str)

  fitted_model <- tryCatch({
    if (is_mixed) {
      lmerTest::lmer(stats::as.formula(lmm_formula_str), data = model_data)
    } else {
      stats::lm(stats::as.formula(lmm_formula_str), data = model_data)
    }
  }, error = function(e) {
    warning(paste("Model fitting failed:", e$message))
    NULL
  })

  if (is.null(fitted_model)) {
    return(list(descriptive = descriptive, model_anova = NULL, posthoc = NULL,
                diagnostics = NULL))
  }

  # --- Model ANOVA ---
  model_anova <- as.data.frame(stats::anova(fitted_model))
  model_anova$term <- rownames(model_anova)

  # Calculate effect size (partial eta-squared).
  #
  # partial eta^2 = SS_effect / (SS_effect + SS_error), so it needs SS_error. An
  # lmerTest::anova table has NO Residuals row -- its columns are Sum Sq, Mean Sq,
  # NumDF, DenDF, F value, Pr(>F) -- so SS_error must come from the residuals. An
  # earlier version divided by sum(Sum Sq), i.e. the total over FIXED TERMS ONLY,
  # which is each term's share of the fixed-effect SS and not an eta^2 at all; it
  # overstated the effect roughly tenfold (Site: 0.885 vs 0.090 for Shannon).
  #
  # For lm this is the textbook quantity. For a mixed model it is approximate:
  # residuals() gives conditional residuals, so SS_error excludes the between-patient
  # variance absorbed by the random intercept, and the value should be read as the
  # share of WITHIN-patient residual variance explained. There is no consensus
  # definition of partial eta^2 for mixed models; prefer the post-hoc contrast with
  # its 95% CI (the `posthoc` slot) when reporting effect magnitude.
  if ("Sum Sq" %in% colnames(model_anova)) {
    ss_error <- sum(stats::residuals(fitted_model)^2)
    model_anova$partial_eta_sq <-
      model_anova$`Sum Sq` / (model_anova$`Sum Sq` + ss_error)
    # anova.lm carries a Residuals row (anova.lmerTest does not); an effect size for
    # the error term is meaningless and would evaluate to a constant 0.5.
    model_anova$partial_eta_sq[model_anova$term == "Residuals"] <- NA_real_
  }

  # Rename columns for consistency
  col_map <- c(
    "Df" = "df",
    "NumDF" = "df_num",
    "DenDF" = "df_den",
    "Sum Sq" = "sum_sq",
    "Mean Sq" = "mean_sq",
    "F value" = "f_statistic",
    "Pr(>F)" = "p_value"
  )
  hit <- colnames(model_anova) %in% names(col_map)
  colnames(model_anova)[hit] <- col_map[colnames(model_anova)[hit]]

  model_anova <- model_anova |>
    dplyr::select(term, dplyr::everything()) |>
    dplyr::mutate(significance = classify_significance(p_value))

  # --- Model Diagnostics ---
  # Both tests run on residuals, not raw y_var: the assumption is constant variance of
  # the error term, and raw within-group spread is inflated by the covariates and the
  # patient random intercept, which the model already fits.
  diagnostics <- tryCatch({
    resid_vec <- stats::residuals(fitted_model)
    grp_vec <- droplevels(factor(model_data[[group_var]]))

    if (length(resid_vec) != length(grp_vec)) {
      stop("residual/group length mismatch")
    }

    levene_result <- car::leveneTest(y = resid_vec, group = grp_vec, center = stats::median)
    shapiro_result <- stats::shapiro.test(resid_vec)

    list(
      levene = data.frame(
        group_var = group_var,
        df_num = levene_result$Df[1],
        df_den = levene_result$Df[2],
        f_statistic = levene_result$`F value`[1],
        p_value = levene_result$`Pr(>F)`[1]
      ) |>
        dplyr::mutate(significance = classify_significance(p_value)),
      normality = data.frame(
        group_var = group_var,
        statistic = unname(shapiro_result$statistic),
        p_value = shapiro_result$p.value,
        n = length(resid_vec)
      ) |>
        dplyr::mutate(significance = classify_significance(p_value)),
      residuals = data.frame(
        .fitted = stats::fitted(fitted_model),
        .resid = resid_vec,
        .std_resid = resid_vec / stats::sd(resid_vec),
        group = grp_vec
      )
    )
  }, error = function(e) {
    warning(paste("Diagnostics failed:", e$message))
    NULL
  })

  # --- Post-hoc Comparisons ---
  posthoc_df <- if (posthoc) tryCatch({
    emm_obj <- emmeans::emmeans(fitted_model, specs = stats::as.formula(paste("pairwise ~", group_var)),
                                lmer.df = "satterthwaite")
    as.data.frame(summary(emm_obj$contrasts, adjust = post_hoc_method, infer = c(TRUE, TRUE))) |>
      dplyr::mutate(
        significance = classify_significance(p.value),
        formatted = sprintf(
          paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
          estimate, lower.CL, upper.CL
        )
      ) |>
      dplyr::rename(comparison = contrast, p_adj = p.value)
  }, error = function(e) {
    warning(paste("Post-hoc failed:", e$message))
    NULL
  })

  list(
    descriptive = descriptive,
    model_anova = model_anova,
    posthoc = posthoc_df,
    diagnostics = diagnostics
  )
}

#' Split "A - B" contrast labels into their two levels
#'
#' emmeans wraps levels containing special characters in parentheses, e.g. "(Re-stenting)";
#' those are stripped. Levels that themselves contain " - " are not supported.
#'
#' @param comparison Character vector of contrast labels.
#' @return Data frame with `group1` and `group2`.
#' @keywords internal
split_contrast <- function(comparison) {
  parts <- strsplit(as.character(comparison), " - ", fixed = TRUE)
  unwrap <- function(x) gsub("^\\((.*)\\)$", "\\1", trimws(x))
  data.frame(group1 = unwrap(vapply(parts, `[`, character(1), 1)),
             group2 = unwrap(vapply(parts, `[`, character(1), 2)))
}

#' Which contrasts have `level` on either side (exact match, not substring)
#' @keywords internal
contrast_involves <- function(comparison, level) {
  g <- split_contrast(comparison)
  g$group1 == level | g$group2 == level
}
