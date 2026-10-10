#' Prepare a TreeSummarizedExperiment for mikropml ML training
#'
#' Subsets a TSE to the requested diagnostic classes (optionally relabelling them
#' into a coarser binary outcome), applies a prevalence filter at the feature
#' level, CLR-transforms the filtered counts, and returns a flat data frame
#' ready to pass to `mikropml::run_ml()`.
#'
#' CLR is recomputed on the subset/filtered table to avoid leakage from features
#' that were dropped by the prevalence filter.
#'
#' @param tse A TreeSummarizedExperiment (e.g. `tses$gut`).
#' @param target_levels Character vector of `Diagnosis` levels to retain
#'   (e.g. `c("Control", "Carcinoma")` or all four levels).
#' @param target_relabel Optional named list mapping new class labels to vectors
#'   of original Diagnosis levels. Example:
#'   `list(Control = "Control", Disease = c("Adenoma", "Advanced adenoma", "Carcinoma"))`.
#'   If NULL, original labels are used (possibly with spaces collapsed to "_").
#' @param prevalence Numeric in [0, 1]. Features must be present (count > 0) in
#'   at least this fraction of subset samples to be retained. Default 0.2.
#' @param assay_type Source assay for CLR. Default "counts".
#' @param pseudocount Pseudocount added before CLR. Default 1.
#' @param outcome_colname Name of the outcome column in the returned data frame.
#'   Default "dx".
#' @param aggregate_rank Optional. If NULL (default), use the TSE as-is
#'   (ASV-level). If a character like `"genus"`, agglomerate the TSE to that
#'   rank via `mia::agglomerateByRank()` before subsetting/CLR. Used for
#'   portable production models whose features must align across cohorts.
#'
#' @return A list with:
#'   * `df` — data.frame whose first column (`outcome_colname`) is a factor
#'     outcome and remaining columns are CLR-transformed feature values
#'     (one column per retained feature, named by feature ID).
#'   * `groups` — character vector of patient IDs aligned to `df` rows
#'     (used by `mikropml::run_ml(groups = ...)` for patient-blocked CV).
#'   * `taxonomy` — data.frame of `rowData(tse)` restricted to retained features,
#'     with a `feature_id` column for joining onto importance results.
#'   * `meta` — list with `n_samples`, `n_features_in`, `n_features_out`,
#'     `class_counts`, `prevalence`, `target_levels`, `target_relabel`.
prepare_ml_data <- function(
  tse,
  target_levels,
  target_relabel = NULL,
  prevalence = 0.2,
  assay_type = "counts",
  pseudocount = 1,
  outcome_colname = "dx",
  aggregate_rank = NULL
) {
  stopifnot(
    inherits(tse, "SummarizedExperiment"),
    is.character(target_levels), length(target_levels) >= 2,
    is.numeric(prevalence), prevalence >= 0, prevalence <= 1
  )

  # --- Optional pre-agglomeration to a taxonomic rank ---
  if (!is.null(aggregate_rank)) tse <- agglomerate_for_ml(tse, aggregate_rank)

  # --- Row subset (samples) ---
  cd <- SummarizedExperiment::colData(tse)
  if (!"Diagnosis" %in% names(cd)) {
    stop("colData(tse) must contain a 'Diagnosis' column.")
  }
  if (!"Patient" %in% names(cd)) {
    stop("colData(tse) must contain a 'Patient' column.")
  }
  keep_samples <- as.character(cd$Diagnosis) %in% target_levels
  if (sum(keep_samples) == 0) {
    stop("No samples match target_levels: ", paste(target_levels, collapse = ", "))
  }
  tse_sub <- tse[, keep_samples]

  # --- Relabel outcome ---
  orig <- as.character(SummarizedExperiment::colData(tse_sub)$Diagnosis)
  if (is.null(target_relabel)) {
    y <- gsub("[^A-Za-z0-9]+", "_", orig)
    class_levels <- gsub("[^A-Za-z0-9]+", "_", target_levels)
  } else {
    y <- rep(NA_character_, length(orig))
    for (new_label in names(target_relabel)) {
      y[orig %in% target_relabel[[new_label]]] <- new_label
    }
    if (anyNA(y)) {
      stop(
        "target_relabel does not cover all target_levels. Unlabelled: ",
        paste(unique(orig[is.na(y)]), collapse = ", ")
      )
    }
    class_levels <- names(target_relabel)
  }
  y <- factor(y, levels = class_levels)
  if (nlevels(y) < 2 || any(table(y) < 2)) {
    stop("Outcome must have at least 2 classes each with ≥2 samples. Got: ",
         paste(paste0(names(table(y)), "=", table(y)), collapse = ", "))
  }

  # --- Prevalence filter (on subset) ---
  counts <- SummarizedExperiment::assay(tse_sub, assay_type)
  n_samp <- ncol(counts)
  prev_n <- Matrix::rowSums(counts > 0)
  keep_feat <- prev_n >= prevalence * n_samp
  if (sum(keep_feat) < 2) {
    stop("Prevalence filter retained <2 features. Lower 'prevalence' or check assay.")
  }
  counts_f <- counts[keep_feat, , drop = FALSE]

  # --- CLR transform (per-sample, on filtered counts) ---
  # log(x + pseudocount) - rowMean(log(x + pseudocount))   computed per sample
  logc <- log(as.matrix(counts_f) + pseudocount)
  sample_means <- colMeans(logc)
  clr_mat <- sweep(logc, 2, sample_means, FUN = "-")

  # --- Assemble df ---
  X <- t(clr_mat)  # samples × features
  colnames(X) <- rownames(counts_f)
  df <- data.frame(
    .outcome = y,
    X,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(df)[1] <- outcome_colname

  groups <- as.character(SummarizedExperiment::colData(tse_sub)$Patient)

  taxonomy <- as.data.frame(SummarizedExperiment::rowData(tse_sub)[keep_feat, , drop = FALSE])
  taxonomy$feature_id <- rownames(counts_f)

  list(
    df = df,
    groups = groups,
    taxonomy = taxonomy,
    meta = list(
      n_samples = n_samp,
      n_features_in = nrow(counts),
      n_features_out = nrow(counts_f),
      class_counts = as.list(table(y)),
      prevalence = prevalence,
      target_levels = target_levels,
      target_relabel = target_relabel,
      outcome_colname = outcome_colname,
      rank = if (is.null(aggregate_rank)) "asv" else aggregate_rank,
      pseudocount = pseudocount,
      assay_type = assay_type
    )
  )
}

#' Agglomerate to a rank and name features by their taxon at that rank
#'
#' Used both at training (`prepare_ml_data()`) and scoring (`score_ml_model()`) time, so
#' portable models reference biological names rather than ASV hashes and the two sides
#' derive identical feature names. Missing names fall back to the rowname; duplicates get
#' "_1", "_2", ... suffixes.
#'
#' @param tse A (Tree)SummarizedExperiment.
#' @param rank Taxonomic rank column.
#' @return The agglomerated object with renamed rows.
#' @keywords internal
agglomerate_for_ml <- function(tse, rank) {
  tse <- mia::agglomerateByRank(tse, rank = rank, na.rm = FALSE)
  rd <- SummarizedExperiment::rowData(tse)
  if (rank %in% names(rd)) {
    ids <- as.character(rd[[rank]])
    missing <- is.na(ids) | ids == ""
    ids[missing] <- rownames(tse)[missing]
    rownames(tse) <- make.unique(ids, sep = "_")
  }
  tse
}
