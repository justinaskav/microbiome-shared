#' Plotting helpers for mikropml-based ML results
#'
#' All helpers consume the outputs of `run_ml_task()` for one or several
#' algorithms on the same binary classification task.

# ---------------------------------------------------------------------------
# 0. Shared AUROC helpers --------------------------------------------------
# ---------------------------------------------------------------------------

#' Class levels in the order pROC needs them, derived from the data.
#'
#' Always taken from `levels(truth)` rather than from a caller-supplied negative class:
#' a `negative_class` argument that disagrees with the factor's level order silently
#' flips the ROC orientation and yields `1 - AUC`.
#'
#' @keywords internal
.roc_levels <- function(truth, positive_class) {
  lvls <- levels(truth)
  if (is.null(lvls)) stop("`truth` must be a factor.")
  if (!positive_class %in% lvls) {
    stop("positive_class '", positive_class, "' is not a level of `truth`. Levels: ",
         paste(lvls, collapse = ", "))
  }
  neg <- setdiff(lvls, positive_class)
  if (length(neg) != 1) {
    stop("expected exactly 2 levels in `truth`, found ", length(lvls), ".")
  }
  c(neg, positive_class)
}

#' AUROC per repeat from one algorithm's pooled test predictions.
#'
#' The single definition of "per-repeat AUROC" used by `plot_ml_roc()`,
#' `plot_ml_per_repeat_auroc()` and `tbl_ml_metrics()`, so the figure legend, the
#' distribution plot and the table cannot drift apart. A repeat containing only one
#' class yields `NA` rather than aborting.
#'
#' @param df A `run_ml_task()$test_preds` tibble (needs `truth`, `prob_pos`, `repeat_id`).
#' @param positive_class Character. Positive factor level.
#' @return Named numeric vector of AUROCs, one per repeat.
#' @keywords internal
.per_repeat_auc <- function(df, positive_class) {
  if (is.null(df$repeat_id)) {
    stop("`test_preds` has no `repeat_id` column; per-repeat AUROC cannot be computed. ",
         "Pass the tibble returned by run_ml_task()$test_preds.")
  }
  lv <- .roc_levels(df$truth, positive_class)
  vapply(split(df, df$repeat_id), function(d) {
    if (length(unique(stats::na.omit(d$truth))) < 2) return(NA_real_)
    tryCatch(
      as.numeric(pROC::auc(d$truth, d$prob_pos, levels = lv,
                           direction = "<", quiet = TRUE)),
      error = function(e) NA_real_
    )
  }, numeric(1))
}

# ---------------------------------------------------------------------------
# 1. ROC overlay -----------------------------------------------------------
# ---------------------------------------------------------------------------

#' Overlay per-algorithm ROC curves from pooled held-out test predictions.
#'
#' Pools `prob_pos` across all repeats per algorithm to draw a single ROC curve
#' per method (for curve shape). The AUROC reported in the legend is the **median
#' across the repeats with the IQR** — not a DeLong CI on the pooled predictions,
#' which would be anticonservative because the same patients recur across repeats'
#' test sets (violating DeLong's independence assumption). This matches the
#' median [IQR] reported by `tbl_ml_metrics()`.
#'
#' @param test_preds_list Named list of tibbles returned by
#'   `run_ml_task()$test_preds` (e.g. `list(RF = ..., XGBoost = ..., ElasticNet = ...)`).
#'   Each tibble must carry a `repeat_id` column (errors otherwise, rather than
#'   rendering "NA [NA, NA]" into the legend).
#' @param positive_class Character. Positive factor level (must match across inputs).
#' @param title Optional plot title.
#' @param colors Optional named character vector of colors keyed by method name.
plot_ml_roc <- function(test_preds_list, positive_class, title = NULL,
                         colors = NULL) {
  stopifnot(is.list(test_preds_list), length(test_preds_list) >= 1)

  rocs <- lapply(names(test_preds_list), function(nm) {
    df <- test_preds_list[[nm]]
    roc_obj <- pROC::roc(
      response = df$truth,
      predictor = df$prob_pos,
      levels = .roc_levels(df$truth, positive_class),
      direction = "<",
      quiet = TRUE
    )
    # Across-repeat AUROC distribution (valid summary of between-repeat variability).
    rep_auc <- .per_repeat_auc(df, positive_class)
    rep_auc <- rep_auc[is.finite(rep_auc)]
    if (length(rep_auc) == 0) {
      stop("no repeat of '", nm, "' yielded a computable AUROC.")
    }
    list(method = nm, roc = roc_obj,
         med = stats::median(rep_auc), lo = stats::quantile(rep_auc, 0.25, names = FALSE),
         hi = stats::quantile(rep_auc, 0.75, names = FALSE))
  })

  roc_df <- do.call(rbind, lapply(rocs, function(r) {
    data.frame(
      method = r$method,
      fpr = 1 - r$roc$specificities,
      tpr = r$roc$sensitivities
    )
  }))

  legend_lbls <- vapply(rocs, function(r) {
    sprintf("%s — AUROC %.3f [%.3f, %.3f]", r$method, r$med, r$lo, r$hi)
  }, character(1))
  names(legend_lbls) <- vapply(rocs, `[[`, character(1), "method")

  roc_df$method <- factor(roc_df$method, levels = names(legend_lbls))

  p <- ggplot2::ggplot(roc_df, ggplot2::aes(x = fpr, y = tpr, color = method)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "grey60") +
    ggplot2::geom_step(linewidth = 0.9, direction = "vh") +
    ggplot2::scale_x_continuous("1 − Specificity (FPR)",
                                limits = c(0, 1), expand = c(0.01, 0.01)) +
    ggplot2::scale_y_continuous("Sensitivity (TPR)",
                                limits = c(0, 1), expand = c(0.01, 0.01)) +
    ggplot2::coord_equal() +
    ggplot2::theme(legend.position = c(0.98, 0.02),
                   legend.justification = c(1, 0),
                   legend.background = ggplot2::element_rect(
                     fill = ggplot2::alpha("white", 0.8),
                     color = "grey80"))
  p <- if (is.null(colors)) {
    p + ggplot2::scale_color_discrete(labels = legend_lbls, name = NULL)
  } else {
    p + ggplot2::scale_color_manual(values = colors, labels = legend_lbls,
                                     name = NULL)
  }

  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p
}

# ---------------------------------------------------------------------------
# 2. Feature importance plot ----------------------------------------------
# ---------------------------------------------------------------------------

#' Top-N permutation-importance taxa per algorithm, with taxonomic labels.
#'
#' @param importance_list Named list of tibbles from `run_ml_task()$importance`.
#'   Each tibble must contain at minimum `feat`, `perf_metric_diff`,
#'   `repeat_id`, `method` (the standard mikropml output).
#' @param taxonomy `prep$taxonomy` from `prepare_ml_data()`.
#' @param top_n Number of top taxa per algorithm. Default 15.
#' @param label_rank Taxonomic rank used for the bar label. Falls back to
#'   feature_id if the requested rank is NA.
plot_ml_importance <- function(importance_list, taxonomy,
                                top_n = 15, label_rank = "genus") {
  stopifnot(is.list(importance_list), length(importance_list) >= 1)

  # Overwrite inner $method with list-name for consistent labelling.
  all_imp <- dplyr::bind_rows(
    lapply(names(importance_list), function(nm) {
      d <- importance_list[[nm]]
      d$method <- nm
      d
    })
  )
  if (!"feat" %in% names(all_imp) || !"perf_metric_diff" %in% names(all_imp)) {
    stop("Importance tibbles must contain 'feat' and 'perf_metric_diff'.")
  }

  # Aggregate per (method, feature) across repeats
  agg <- all_imp |>
    dplyr::group_by(method, feat) |>
    dplyr::summarise(
      mean_diff = mean(perf_metric_diff, na.rm = TRUE),
      sd_diff   = stats::sd(perf_metric_diff, na.rm = TRUE),
      n_reps    = dplyr::n(),
      .groups = "drop"
    )

  top <- agg |>
    dplyr::group_by(method) |>
    dplyr::slice_max(mean_diff, n = top_n, with_ties = FALSE) |>
    dplyr::ungroup()

  # Join taxonomy for nice labels
  tax_lookup <- taxonomy
  tax_lookup$.label <- ifelse(
    !is.na(tax_lookup[[label_rank]]) & tax_lookup[[label_rank]] != "",
    tax_lookup[[label_rank]], tax_lookup$feature_id
  )
  # If multiple ASVs share the same label, suffix with short hash
  short_hash <- substr(tax_lookup$feature_id, 1, 6)
  dup <- tax_lookup$.label %in% tax_lookup$.label[duplicated(tax_lookup$.label)]
  tax_lookup$.label <- ifelse(dup,
                              paste0(tax_lookup$.label, " (", short_hash, ")"),
                              tax_lookup$.label)

  top <- top |>
    dplyr::left_join(tax_lookup[, c("feature_id", ".label")],
                     by = c("feat" = "feature_id"))
  top$.label <- ifelse(is.na(top$.label), top$feat, top$.label)

  # Order bars within each facet by mean_diff
  top <- top |>
    dplyr::mutate(
      facet_rank = paste(method, .label, sep = "::"),
      .label = forcats::fct_reorder(.label, mean_diff)
    )

  ggplot2::ggplot(top, ggplot2::aes(x = mean_diff, y = .label, fill = method)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::facet_wrap(~ method, scales = "free_y") +
    ggplot2::labs(
      x = "Mean Δ AUROC when feature permuted",
      y = NULL,
      caption = paste0("Top ", top_n,
                       " features per algorithm, aggregated over repeats.")
    ) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"))
}

# ---------------------------------------------------------------------------
# 3. Metrics summary table ------------------------------------------------
# ---------------------------------------------------------------------------

#' Build a `gt` table summarising performance per algorithm.
#'
#' Aggregates per-repeat test predictions into one set of metrics per
#' (algorithm × repeat), then reports median [IQR] across repeats. mikropml's
#' performance tibble is used only for the inner-CV training AUC
#' (`cv_metric_AUC`) — test metrics are computed here from the held-out
#' predictions to avoid mikropml's idiosyncratic `multiClassSummary` defaults.
#'
#' Binary metrics are computed at each repeat's `threshold` column (Youden's J on the inner-CV
#' predictions, see `run_ml_task()`); results without that column fall back to 0.5. The cut-off
#' is never chosen on the test predictions themselves, which would bias the metrics upwards.
#'
#' @param performance_list Named list of tibbles (`run_ml_task()$performance`)
#'   used only for `cv_metric_AUC`.
#' @param test_preds_list Named list of tibbles (`run_ml_task()$test_preds`).
#'   Names must match `performance_list`.
#' @param positive_class Character. Positive class level.
#' @param title Caption for the table.
tbl_ml_metrics <- function(performance_list, test_preds_list, positive_class,
                            title = "Model performance") {
  stopifnot(is.list(performance_list), is.list(test_preds_list),
            identical(sort(names(performance_list)), sort(names(test_preds_list))))

  med_iqr <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA_character_)
    sprintf("%.3f [%.3f, %.3f]",
            stats::median(x),
            stats::quantile(x, 0.25, names = FALSE),
            stats::quantile(x, 0.75, names = FALSE))
  }

  # Use the list-name as the canonical method label (the inner $method column
  # from mikropml is the caret method code, e.g. "glmnet", which may differ
  # from the user-facing display name like "ElasticNet").
  cv_summ <- dplyr::bind_rows(
    lapply(names(performance_list), function(nm) {
      d <- performance_list[[nm]]
      d$method <- nm
      d
    })
  ) |>
    dplyr::group_by(method) |>
    dplyr::summarise(`CV AUROC` = med_iqr(cv_metric_AUC), .groups = "drop")

  # Per-repeat metrics, computed outside dplyr group_modify (which has had
  # scoping issues with closures in this code base).
  per_repeat <- lapply(names(test_preds_list), function(nm) {
    df <- test_preds_list[[nm]]
    reps <- unique(df$repeat_id)
    rows <- lapply(reps, function(r) {
      sub <- df[df$repeat_id == r, , drop = FALSE]
      m <- .compute_binary_metrics(sub, positive_class)
      m$method <- nm
      m$repeat_id <- r
      m
    })
    dplyr::bind_rows(rows)
  })
  test_summ <- dplyr::bind_rows(per_repeat) |>
    dplyr::group_by(method) |>
    dplyr::summarise(
      `Test AUROC`   = med_iqr(auc),
      `F1`           = med_iqr(f1),
      `Sensitivity`  = med_iqr(sens),
      `Specificity`  = med_iqr(spec),
      `Balanced acc.`= med_iqr(bal_acc),
      .groups = "drop"
    )

  summ <- dplyr::left_join(cv_summ, test_summ, by = "method") |>
    dplyr::rename(Model = method)

  cutoff_note <- if (all(vapply(test_preds_list, function(d) "threshold" %in% names(d), logical(1)))) {
    "Binary metrics computed at Youden's J cut-off chosen on each repeat's inner-CV predictions."
  } else {
    "Binary metrics computed at a 0.5 probability cut-off."
  }
  gt::gt(summ) |>
    gt::tab_header(title = title) |>
    gt::tab_source_note(
      gt::md(paste("Median [IQR] across repeats.", cutoff_note,
                   "CV AUROC = inner cross-validation; gap vs Test AUROC indicates overfitting."))
    )
}

# ---------------------------------------------------------------------------
# 4. Per-repeat AUROC distribution ----------------------------------------
# ---------------------------------------------------------------------------

#' Boxplot of held-out test AUROC across repeats, one box per algorithm.
#'
#' Computes per-repeat AUROC from pooled test predictions and plots the
#' distribution. Mirrors the per-repeat block previously inlined in each ML
#' section.
#'
#' @param test_preds_list Named list of `run_ml_task()$test_preds` tibbles.
#' @param positive_class Character. Positive factor level. The negative level is derived from
#'   `levels(truth)` rather than taken as an argument, so the orientation cannot disagree with
#'   `plot_ml_roc()` / `tbl_ml_metrics()`.
#' @param ylim Numeric length-2 y-axis limits (default c(0.4, 1.0)).
plot_ml_per_repeat_auroc <- function(test_preds_list, positive_class,
                                     ylim = c(0.4, 1.0)) {
  per_rep <- dplyr::bind_rows(lapply(names(test_preds_list), function(nm) {
    aucs <- .per_repeat_auc(test_preds_list[[nm]], positive_class)
    data.frame(method = nm, repeat_id = names(aucs), auc = unname(aucs))
  }))

  ggplot2::ggplot(per_rep, ggplot2::aes(x = method, y = auc, color = method)) +
    ggplot2::geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.4, na.rm = TRUE) +
    ggplot2::geom_jitter(width = 0.1, height = 0, size = 2, alpha = 0.8, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60") +
    ggplot2::scale_y_continuous("Held-out test AUROC", limits = ylim) +
    ggplot2::labs(x = NULL) +
    ggplot2::theme(legend.position = "none")
}

# ---------------------------------------------------------------------------
# 5. Paired AUROC comparison across repeats ---------------------------------
# ---------------------------------------------------------------------------

#' Paired per-repeat AUROC comparison of two models evaluated on the same splits.
#'
#' Replaces a DeLong test on predictions pooled across repeats: the same patients recur in
#' several repeats' test sets, so pooled predictions are not independent and DeLong's p-value
#' would be anticonservative. Each repeat instead contributes one paired difference
#' (AUROC of `b` minus AUROC of `a` on the identical test set).
#'
#' @param a_preds,b_preds `run_ml_task()$test_preds` tibbles run with the same seeds/splits.
#' @param pos Positive class label; the negative level is derived from `levels(truth)`.
#' @return A list: `per_repeat` (repeat_id, auc_a, auc_b, diff) and `summary` (one-row
#'   tibble: median AUROC of each model, median and IQR of the paired difference, number of
#'   repeats, and the share of repeats in which `b` beat `a`).
compare_auroc_by_repeat <- function(a_preds, b_preds, pos) {
  auc_a <- .per_repeat_auc(a_preds, pos)
  auc_b <- .per_repeat_auc(b_preds, pos)
  ids <- intersect(names(auc_a), names(auc_b))
  per_repeat <- tibble::tibble(repeat_id = ids, auc_a = auc_a[ids], auc_b = auc_b[ids],
                               diff = auc_b[ids] - auc_a[ids])
  d <- per_repeat$diff[is.finite(per_repeat$diff)]
  list(
    per_repeat = per_repeat,
    summary = tibble::tibble(
      auc_a = stats::median(per_repeat$auc_a, na.rm = TRUE),
      auc_b = stats::median(per_repeat$auc_b, na.rm = TRUE),
      diff_median = stats::median(d),
      diff_q25 = stats::quantile(d, 0.25, names = FALSE),
      diff_q75 = stats::quantile(d, 0.75, names = FALSE),
      n_repeats = length(d),
      share_b_better = mean(d > 0)
    )
  )
}

#' @keywords internal
.compute_binary_metrics <- function(df, positive_class) {
  truth <- df$truth
  prob  <- df$prob_pos
  lv <- .roc_levels(truth, positive_class)
  neg_class <- lv[1]
  if (length(unique(stats::na.omit(truth))) < 2) {
    return(tibble::tibble(auc = NA_real_, f1 = NA_real_,
                          sens = NA_real_, spec = NA_real_, bal_acc = NA_real_))
  }
  r <- pROC::roc(truth, prob, levels = lv, direction = "<", quiet = TRUE)
  thr <- if (!is.null(df$threshold)) df$threshold[1] else 0.5
  pred <- factor(ifelse(prob >= thr, positive_class, neg_class),
                 levels = c(neg_class, positive_class))
  tp <- sum(pred == positive_class & truth == positive_class)
  fp <- sum(pred == positive_class & truth == neg_class)
  fn <- sum(pred == neg_class      & truth == positive_class)
  tn <- sum(pred == neg_class      & truth == neg_class)
  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  f1   <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0)
    2 * prec * sens / (prec + sens) else NA_real_
  bal_acc <- mean(c(sens, spec), na.rm = TRUE)
  tibble::tibble(auc = as.numeric(r$auc), f1 = f1, sens = sens,
                 spec = spec, bal_acc = bal_acc)
}
