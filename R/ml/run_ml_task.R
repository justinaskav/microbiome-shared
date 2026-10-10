#' Run a mikropml binary-classification task across multiple seeds
#'
#' Repeatedly trains a single algorithm on patient-blocked train/test splits and
#' aggregates performance, feature importance, and held-out predicted
#' probabilities across repeats. Patient blocking is delegated to
#' `mikropml::run_ml(groups = ...)`, which uses `caret::groupKFold` internally
#' for both the train/test split and the inner CV used for hyperparameter
#' tuning.
#'
#' @param prep A list returned by `prepare_ml_data()`.
#' @param method One of mikropml's supported caret methods. We use
#'   `"glmnet"` (elastic-net logistic), `"rf"` (random forest via ranger), or
#'   `"xgbTree"` (XGBoost).
#' @param repeats Integer. Number of independent train/test splits (different
#'   seeds). Each repeat is one held-out test AUROC.
#' @param kfold Integer. Inner-CV folds for hyperparameter tuning. Default 5.
#' @param cv_times Integer. Inner-CV repetitions. Default 5.
#' @param training_frac Fraction of *patients* allocated to training in each
#'   repeat. Default 0.8.
#' @param hyperparameters Optional named list of tuning grids (passed straight
#'   through to `run_ml`). If NULL, the package defaults are used.
#' @param base_seed Integer. Repeat r uses `seed = base_seed + r`.
#' @param find_feature_importance Logical. Whether to compute permutation
#'   importance (slow but uniform across algorithms). Computed on the held-out
#'   test set per repeat by permuting each feature `nperms` times and measuring
#'   the resulting AUROC drop.
#' @param importance_nperms Integer. Permutations per feature for importance.
#' @param positive_class Optional character. The factor level treated as
#'   "positive" when computing predicted probabilities. Defaults to the second
#'   level of `prep$df[[outcome_colname]]`.
#'
#' @return A list with:
#'   * `performance` — tibble; one row per repeat with all mikropml metrics
#'     plus `repeat_id`, `method`.
#'   * `importance` — tibble; one row per (feature × repeat) from mikropml's
#'     permutation importance, plus `repeat_id`, `method`.
#'   * `test_preds` — tibble; one row per (test sample × repeat) with
#'     `repeat_id`, `method`, `truth` (factor), `prob_pos` (numeric in [0,1]) and
#'     `threshold`: that repeat's Youden's J cut-off, chosen on the inner-CV
#'     (out-of-fold, training-only) predictions at the selected hyperparameters, so
#'     binary test metrics do not use the test labels to pick the cut-off.
#'   * `meta` — settings echoed back for the plotting helpers.
run_ml_task <- function(
  prep,
  method,
  repeats = 10,
  kfold = 5,
  cv_times = 5,
  training_frac = 0.8,
  hyperparameters = NULL,
  base_seed = 2026,
  find_feature_importance = TRUE,
  importance_nperms = 20,
  positive_class = NULL
) {
  stopifnot(is.list(prep), all(c("df", "groups", "meta") %in% names(prep)))
  outcome <- prep$meta$outcome_colname %||% "dx"
  y <- prep$df[[outcome]]
  if (!is.factor(y) || nlevels(y) != 2) {
    stop("This wrapper supports binary classification only.")
  }
  if (is.null(positive_class)) positive_class <- levels(y)[2]
  if (!positive_class %in% levels(y)) {
    stop("positive_class '", positive_class, "' not in outcome levels.")
  }

  perf_rows <- vector("list", repeats)
  imp_rows  <- vector("list", repeats)
  pred_rows <- vector("list", repeats)

  for (r in seq_len(repeats)) {
    seed_r <- base_seed + r
    parts <- .stratified_group_partition(y, prep$groups, training_frac, seed_r)
    # mikropml's own feature_importance has an upstream bug (factor→character
    # coercion in calc_perf_metrics yields a degenerate 4-level outcome). We
    # always disable it and compute permutation importance ourselves below.
    res <- mikropml::run_ml(
      dataset = prep$df,
      method = method,
      outcome_colname = outcome,
      hyperparameters = hyperparameters,
      find_feature_importance = FALSE,
      kfold = kfold,
      cv_times = cv_times,
      training_frac = training_frac,
      groups = prep$groups,
      group_partitions = parts,
      seed = seed_r
    )

    perf_rows[[r]] <- tibble::as_tibble(res$performance) |>
      dplyr::mutate(repeat_id = r, method = method)

    test_data <- res$test_data
    truth <- factor(test_data[[outcome]], levels = levels(y))
    probs <- stats::predict(
      res$trained_model, newdata = test_data, type = "prob"
    )
    pred_rows[[r]] <- tibble::tibble(
      repeat_id = r,
      method = method,
      truth = truth,
      prob_pos = probs[[positive_class]],
      threshold = .cv_youden_threshold(res$trained_model, levels(y), positive_class)
    )

    if (find_feature_importance) {
      imp_rows[[r]] <- .permutation_importance(
        trained_model = res$trained_model,
        test_data = test_data,
        outcome_colname = outcome,
        positive_class = positive_class,
        nperms = importance_nperms,
        seed = seed_r
      ) |> dplyr::mutate(repeat_id = r, method = method)
    }
  }

  list(
    performance = dplyr::bind_rows(perf_rows),
    importance  = if (find_feature_importance) dplyr::bind_rows(imp_rows) else NULL,
    test_preds  = dplyr::bind_rows(pred_rows),
    meta = list(
      method = method,
      repeats = repeats,
      kfold = kfold,
      cv_times = cv_times,
      training_frac = training_frac,
      positive_class = positive_class,
      negative_class = setdiff(levels(y), positive_class),
      base_seed = base_seed
    )
  )
}


# Youden's J cut-off from caret's saved out-of-fold predictions (mikropml trains with
# savePredictions = TRUE), restricted to the selected hyperparameters.
.cv_youden_threshold <- function(trained_model, class_levels, positive_class) {
  cv <- trained_model$pred
  for (p in names(trained_model$bestTune)) cv <- cv[cv[[p]] == trained_model$bestTune[[p]], , drop = FALSE]
  r <- pROC::roc(factor(cv$obs, levels = class_levels), cv[[positive_class]],
                 levels = c(setdiff(class_levels, positive_class), positive_class),
                 direction = "<", quiet = TRUE)
  pROC::coords(r, "best", best.method = "youden", ret = "threshold", transpose = FALSE)$threshold[1]
}


# Permutation-importance on the held-out test set.
# For each feature: permute it `nperms` times, recompute test AUROC, average the
# drop from the baseline. Replaces mikropml::get_feature_importance which has
# an upstream factor-coercion bug. Output columns match the subset mikropml
# returns that the plotting helper needs: `feat`, `perf_metric_diff`.
.permutation_importance <- function(trained_model, test_data, outcome_colname,
                                     positive_class, nperms = 20, seed = 1) {
  set.seed(seed)
  truth <- test_data[[outcome_colname]]
  neg_class <- setdiff(levels(truth), positive_class)
  baseline_probs <- stats::predict(trained_model, test_data, type = "prob")[[positive_class]]
  baseline_auc <- as.numeric(pROC::auc(truth, baseline_probs,
                                        levels = c(neg_class, positive_class),
                                        direction = "<", quiet = TRUE))

  features <- setdiff(names(test_data), outcome_colname)
  rows <- lapply(features, function(f) {
    drops <- vapply(seq_len(nperms), function(p) {
      td <- test_data
      td[[f]] <- sample(td[[f]])
      pp <- stats::predict(trained_model, td, type = "prob")[[positive_class]]
      a <- tryCatch(
        as.numeric(pROC::auc(truth, pp,
                              levels = c(neg_class, positive_class),
                              direction = "<", quiet = TRUE)),
        error = function(e) NA_real_
      )
      baseline_auc - a
    }, numeric(1))
    data.frame(
      feat = f,
      perf_metric_diff = mean(drops, na.rm = TRUE),
      sd_diff = stats::sd(drops, na.rm = TRUE),
      baseline_auc = baseline_auc,
      nperms = nperms,
      stringsAsFactors = FALSE
    )
  })
  tibble::as_tibble(do.call(rbind, rows))
}


# Stratified train/test split at the *patient* level so that each test set
# contains members of both classes. mikropml's default random group split can
# produce single-class test sets when the number of unique patients is small.
.stratified_group_partition <- function(y, groups, training_frac, seed) {
  stopifnot(length(y) == length(groups), training_frac > 0, training_frac < 1)
  set.seed(seed)
  # Map each patient to a class (assume one class per patient, which is the
  # case here — patient is nested within Diagnosis).
  patient_class <- tapply(as.character(y), groups, function(v) {
    u <- unique(v); if (length(u) > 1) NA_character_ else u
  })
  if (anyNA(patient_class)) {
    stop("Some patients have multiple class labels; stratified split assumes one class per patient.")
  }
  train_pat <- character(0)
  for (cls in unique(patient_class)) {
    pats <- names(patient_class)[patient_class == cls]
    n_train <- max(1, round(length(pats) * training_frac))
    n_train <- min(n_train, length(pats) - 1)  # leave ≥1 patient for test
    train_pat <- c(train_pat, sample(pats, n_train))
  }
  list(
    train = train_pat,
    test  = setdiff(names(patient_class), train_pat)
  )
}
