#' Encode clinical/lifestyle covariates for ML
#'
#' Numeric variables (including integer): median-impute, z-score.
#' Logical / character / factor variables: mode-impute, one-hot encode with the
#' first level (alphabetical for character, level order for factor, FALSE for
#' logical) dropped as the reference.
#'
#' All imputation values, scaling parameters and one-hot level sets are stored
#' in the returned `meta` so the same transform can be re-applied to an
#' external cohort at scoring time.
#'
#' @param tse A SummarizedExperiment / TreeSummarizedExperiment.
#' @param samples_idx Logical or integer vector indexing into `colData(tse)`
#'   to select the samples used for fitting. If NULL, all samples are used.
#' @param vars Character vector of column names in `colData(tse)`.
#' @param fit_meta Optional. If supplied (from a previous call), reuse stored
#'   imputation values / centres / scales / onehot levels instead of computing
#'   new ones — used at scoring time on an external cohort.
#'
#' @return A list with:
#'   * `X` — numeric matrix (samples × encoded covariates), no NAs.
#'   * `meta` — list(vars, numeric_vars, categorical_vars, impute_vals,
#'              centers, scales, onehot_levels, encoded_names).
prepare_ml_covariates <- function(tse, samples_idx = NULL, vars,
                                   fit_meta = NULL) {
  stopifnot(is.character(vars), length(vars) >= 1)
  cd <- as.data.frame(SummarizedExperiment::colData(tse))
  if (is.null(samples_idx)) samples_idx <- seq_len(nrow(cd))
  cd <- cd[samples_idx, , drop = FALSE]

  missing <- setdiff(vars, names(cd))
  if (length(missing)) {
    stop("Covariate(s) not found in colData: ", paste(missing, collapse = ", "))
  }

  numeric_vars     <- character(0)
  categorical_vars <- character(0)
  for (v in vars) {
    x <- cd[[v]]
    if (is.numeric(x)) numeric_vars <- c(numeric_vars, v)
    else categorical_vars <- c(categorical_vars, v)
  }

  impute_vals   <- list()
  centers       <- list()
  scales        <- list()
  onehot_levels <- list()

  # ---- numeric: median-impute, z-score -----------------------------------
  num_mat <- matrix(NA_real_, nrow = nrow(cd), ncol = length(numeric_vars))
  colnames(num_mat) <- numeric_vars
  for (v in numeric_vars) {
    x <- as.numeric(cd[[v]])
    if (!is.null(fit_meta)) {
      med <- fit_meta$impute_vals[[v]]
      ctr <- fit_meta$centers[[v]]
      scl <- fit_meta$scales[[v]]
    } else {
      med <- stats::median(x, na.rm = TRUE)
      x_imp <- x; x_imp[is.na(x_imp)] <- med
      ctr <- mean(x_imp, na.rm = TRUE)
      scl <- stats::sd(x_imp, na.rm = TRUE)
      if (!is.finite(scl) || scl == 0) scl <- 1
    }
    x[is.na(x)] <- med
    num_mat[, v] <- (x - ctr) / scl
    impute_vals[[v]] <- med; centers[[v]] <- ctr; scales[[v]] <- scl
  }

  # ---- categorical: mode-impute, one-hot ---------------------------------
  cat_blocks <- list()
  for (v in categorical_vars) {
    x <- cd[[v]]
    if (is.logical(x)) x <- ifelse(x, "TRUE", "FALSE")
    x <- as.character(x)
    if (!is.null(fit_meta)) {
      mode_val <- fit_meta$impute_vals[[v]]
      lvls     <- fit_meta$onehot_levels[[v]]
    } else {
      tab <- sort(table(x, useNA = "no"), decreasing = TRUE)
      mode_val <- if (length(tab)) names(tab)[1] else NA_character_
      lvls <- sort(unique(stats::na.omit(x)))
    }
    x[is.na(x) | x == ""] <- mode_val
    impute_vals[[v]] <- mode_val
    onehot_levels[[v]] <- lvls
    # drop first level as reference
    keep_lvls <- if (length(lvls) > 1) lvls[-1] else character(0)
    if (length(keep_lvls) == 0) next
    mat <- vapply(keep_lvls, function(L) as.numeric(x == L), numeric(length(x)))
    if (!is.matrix(mat)) mat <- matrix(mat, ncol = length(keep_lvls))
    colnames(mat) <- paste0(v, "_", keep_lvls)
    cat_blocks[[v]] <- mat
  }
  cat_mat <- if (length(cat_blocks)) do.call(cbind, cat_blocks) else
    matrix(numeric(0), nrow = nrow(cd), ncol = 0)

  X <- cbind(num_mat, cat_mat)
  # Sanitise names: caret-friendly
  colnames(X) <- gsub("[^A-Za-z0-9_]+", "_", colnames(X))
  rownames(X) <- rownames(cd)

  list(
    X = X,
    meta = list(
      vars = vars,
      numeric_vars = numeric_vars,
      categorical_vars = categorical_vars,
      impute_vals = impute_vals,
      centers = centers,
      scales = scales,
      onehot_levels = onehot_levels,
      encoded_names = colnames(X)
    )
  )
}
