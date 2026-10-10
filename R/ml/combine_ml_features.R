#' Merge encoded covariates into a `prep` list from `prepare_ml_data()`.
#'
#' Row order of `covariates$X` must match `prep$df` row order (the caller is
#' responsible for selecting the same `samples_idx` when calling
#' `prepare_ml_covariates`).
#'
#' @param prep List from `prepare_ml_data()`.
#' @param covariates List from `prepare_ml_covariates()` (`$X`, `$meta`).
#' @return A new `prep` with the augmented `df`, an extra `covariate_meta`
#'   element, and `meta$has_covariates = TRUE`.
combine_ml_features <- function(prep, covariates) {
  stopifnot(is.list(prep), is.list(covariates),
            all(c("df", "groups", "meta") %in% names(prep)),
            all(c("X", "meta") %in% names(covariates)))
  if (nrow(covariates$X) != nrow(prep$df)) {
    stop("Row count mismatch: prep$df has ", nrow(prep$df),
         " rows, covariates$X has ", nrow(covariates$X), ".")
  }
  clash <- intersect(colnames(prep$df), colnames(covariates$X))
  if (length(clash)) {
    stop("Covariate columns clash with feature columns: ",
         paste(clash, collapse = ", "))
  }
  df2 <- cbind(prep$df, as.data.frame(covariates$X, check.names = FALSE))
  prep$df <- df2
  prep$covariate_meta <- covariates$meta
  prep$meta$has_covariates <- TRUE
  prep$meta$n_covariates <- ncol(covariates$X)
  prep
}
