#' Save / load / score portable ML artifacts.
#'
#' The artifact bundles the trained caret model with the preprocessing
#' metadata needed to reproduce the same feature space on a fresh TSE
#' (taxonomic rank, CLR pseudocount, ordered microbiome feature list, and
#' — if applicable — covariate imputation values / one-hot levels).
#'
#' Scoring requires `agglomerate_for_ml()` (prepare_ml_data.R) and, for covariate models,
#' `prepare_ml_covariates()`.

# ---------------------------------------------------------------------------
# save_ml_model -------------------------------------------------------------
# ---------------------------------------------------------------------------

#' Persist a fit from `fit_final_model()` to disk as an RDS artifact.
#'
#' @param fit List returned by `fit_final_model()`.
#' @param path Output path (typically ending in `.rds`).
save_ml_model <- function(fit, path) {
  stopifnot(is.list(fit), !is.null(fit$model), !is.null(fit$prep_meta))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  feature_names <- setdiff(colnames(fit$model$trainingData), c(".outcome", fit$prep_meta$outcome_colname))
  cov_names <- if (!is.null(fit$covariate_meta)) fit$covariate_meta$encoded_names else character(0)
  microbiome_features <- setdiff(feature_names, cov_names)

  class_levels <- levels(fit$model$trainingData$.outcome %||% fit$model$trainingData[[1]])

  artifact <- list(
    model           = fit$model,
    method          = fit$method,
    rank            = fit$prep_meta$rank,
    features        = microbiome_features,
    pseudocount     = fit$prep_meta$pseudocount %||% 1,
    assay_type      = fit$prep_meta$assay_type %||% "counts",
    class_levels    = class_levels,
    positive_class  = class_levels[2],
    covariate_meta  = fit$covariate_meta,
    training_n      = fit$training_summary$training_n,
    training_summary = fit$training_summary,
    saved_at        = Sys.time()
  )
  saveRDS(artifact, path)
  invisible(path)
}

# ---------------------------------------------------------------------------
# load_ml_model -------------------------------------------------------------
# ---------------------------------------------------------------------------

#' Load an artifact previously written by `save_ml_model()`.
#'
#' Warns if any of the package versions present on this machine differ from
#' those recorded at training time (model objects may not be fully
#' deserialisable across major-version jumps).
load_ml_model <- function(path) {
  artifact <- readRDS(path)
  pv <- artifact$training_summary$package_versions
  if (!is.null(pv)) {
    for (pkg in names(pv)) {
      now <- tryCatch(as.character(utils::packageVersion(pkg)),
                      error = function(e) NA_character_)
      if (!is.na(now) && !is.na(pv[[pkg]]) && now != pv[[pkg]]) {
        warning(sprintf(
          "Package '%s' version differs (saved: %s, current: %s).",
          pkg, pv[[pkg]], now
        ), call. = FALSE)
      }
    }
  }
  artifact
}

# ---------------------------------------------------------------------------
# score_ml_model ------------------------------------------------------------
# ---------------------------------------------------------------------------

#' Score a new TSE using a saved artifact.
#'
#' @param artifact List from `load_ml_model()`.
#' @param new_tse A SummarizedExperiment / TreeSummarizedExperiment to score.
#'
#' @return tibble(sample_id, prob_pos, predicted_class).
score_ml_model <- function(artifact, new_tse) {
  stopifnot(is.list(artifact), inherits(new_tse, "SummarizedExperiment"))

  # 1. Agglomerate the incoming TSE to the artifact's rank, naming features as in training.
  rank <- artifact$rank %||% "asv"
  if (rank != "asv") new_tse <- agglomerate_for_ml(new_tse, rank)

  # 2. Per-sample CLR over exactly the training feature set. prepare_ml_data() centred each
  # sample on its prevalence-filtered features, so the centre here must be taken over the same
  # features: centring over all features of the new TSE would shift every value by a
  # sample-specific constant. Features absent from the new TSE count as zero reads.
  counts <- as.matrix(SummarizedExperiment::assay(new_tse, artifact$assay_type %||% "counts"))
  expected <- artifact$features
  counts_f <- matrix(0, nrow = length(expected), ncol = ncol(counts),
                     dimnames = list(expected, colnames(counts)))
  present <- intersect(expected, rownames(counts))
  counts_f[present, ] <- counts[present, , drop = FALSE]
  if (length(present) < length(expected)) {
    warning(length(expected) - length(present), " of ", length(expected),
            " model features are absent from new_tse and were scored as zero counts.",
            call. = FALSE)
  }
  logc <- log(counts_f + (artifact$pseudocount %||% 1))
  X_mb <- t(sweep(logc, 2, colMeans(logc), FUN = "-"))  # samples x features

  # 3. Re-apply covariate transforms if any.
  if (!is.null(artifact$covariate_meta)) {
    cov_res <- prepare_ml_covariates(
      tse = new_tse,
      samples_idx = NULL,
      vars = artifact$covariate_meta$vars,
      fit_meta = artifact$covariate_meta
    )
    # Align covariate columns to the encoded_names from training.
    enc <- artifact$covariate_meta$encoded_names
    X_cov <- matrix(0, nrow = nrow(X_mb), ncol = length(enc),
                    dimnames = list(rownames(X_mb), enc))
    common <- intersect(enc, colnames(cov_res$X))
    if (length(common)) X_cov[, common] <- cov_res$X[, common, drop = FALSE]
    X_full <- cbind(X_mb, X_cov)
  } else {
    X_full <- X_mb
  }

  newdata <- as.data.frame(X_full, check.names = FALSE)
  probs <- stats::predict(artifact$model, newdata = newdata, type = "prob")
  pos <- artifact$positive_class
  preds <- stats::predict(artifact$model, newdata = newdata)
  tibble::tibble(
    sample_id = rownames(newdata),
    prob_pos = probs[[pos]],
    predicted_class = as.character(preds)
  )
}
