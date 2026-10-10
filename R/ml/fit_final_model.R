#' Train a production model on the full `prep` dataset (no held-out split).
#'
#' Unlike `run_ml_task()` — which carves out a stratified patient-blocked test
#' set per repeat to estimate generalisation — this trains a single model on
#' the entire `prep$df`. Hyperparameters are tuned on that same data by
#' repeated patient-grouped CV (`kfold` x `cv_times`), using mikropml's default
#' grid and AUC metric. The result is ready to be persisted for external scoring.
#'
#' @param prep List from `prepare_ml_data()` (possibly augmented by
#'   `combine_ml_features()`).
#' @param method A caret method supported by mikropml, e.g. "rf" or "glmnet".
#' @param kfold Number of CV folds for tuning. Default 5.
#' @param cv_times Number of CV repeats for tuning. Default 5.
#' @param seed Integer RNG seed. Default 2026.
#' @param hyperparameters Optional named list of tuning grids (mikropml format).
#'   If NULL, `mikropml::get_hyperparams_list()` defaults are used.
#'
#' @return list(model, method, prep_meta, covariate_meta, training_summary).
fit_final_model <- function(prep, method,
                            kfold = 5, cv_times = 5,
                            seed = 2026,
                            hyperparameters = NULL) {
  stopifnot(is.list(prep), all(c("df", "groups", "meta") %in% names(prep)))
  outcome <- prep$meta$outcome_colname %||% "dx"
  y <- prep$df[[outcome]]
  if (!is.factor(y) || nlevels(y) != 2) {
    stop("fit_final_model supports binary classification only.")
  }

  set.seed(seed)
  hp <- hyperparameters %||% mikropml::get_hyperparams_list(prep$df, method)
  cv <- mikropml::define_cv(prep$df, outcome, hp,
                            perf_metric_function = mikropml::get_perf_metric_fn("binary"),
                            class_probs = TRUE, kfold = kfold, cv_times = cv_times,
                            groups = prep$groups)
  final_model <- mikropml::train_model(prep$df, outcome, method, cv,
                                       perf_metric_name = "AUC",
                                       tune_grid = mikropml::get_tuning_grid(hp, method))

  list(
    model = final_model,
    method = method,
    prep_meta = prep$meta,
    covariate_meta = prep$covariate_meta,  # NULL if microbiome-only
    training_summary = list(
      cv_auroc = max(final_model$results$AUC, na.rm = TRUE),
      best_tune = final_model$bestTune,
      training_n = nrow(prep$df),
      date = Sys.Date(),
      package_versions = list(
        mikropml = as.character(utils::packageVersion("mikropml")),
        caret    = as.character(utils::packageVersion("caret")),
        ranger   = tryCatch(as.character(utils::packageVersion("ranger")),
                            error = function(e) NA_character_),
        glmnet   = tryCatch(as.character(utils::packageVersion("glmnet")),
                            error = function(e) NA_character_)
      )
    )
  )
}
