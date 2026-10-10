#' Aggregate a sample-level TSE to one observation per patient
#'
#' Collapses repeated biopsies/samples of a subject into a single patient-level profile so that
#' between-subject effects can be tested without pseudoreplication. Aggregation is by the **mean of
#' per-sample relative abundances** (depth-invariant — a deeply sequenced sample does not dominate),
#' not by summing counts.
#'
#' @param tse A (Tree)SummarizedExperiment with a `counts` assay; `rowTree` is preserved when present
#'   (needed for UniFrac on the aggregated object).
#' @param patient_var colData column identifying the subject (default "Patient").
#' @param warn Warn about colData columns that vary within patient (default TRUE). Set FALSE when
#'   the patient-level object is only plotted, not modelled.
#'
#' Note that the two assays are not two views of the same quantity: `relabundance` is a **mean of
#' proportions** (each sample weighted equally) while `counts` is a **sum** (deeply sequenced
#' samples dominate). Recomputing proportions from `counts` therefore does not reproduce
#' `relabundance` — use the `relabundance` assay for anything distance- or composition-based.
#'
#' @return A TreeSummarizedExperiment with one column per patient, a `relabundance` assay (columns sum
#'   to 1) and a `counts` assay (summed counts, for consumers that expect counts), colData carrying the
#'   first row per patient, and the original `rowData` and `rowTree`.
#' @export
aggregate_to_patient <- function(tse, patient_var = "Patient", warn = TRUE) {
  stopifnot(inherits(tse, "SummarizedExperiment"))
  cd <- SummarizedExperiment::colData(tse)
  # NB: check membership, not is.null(as.character(cd[[...]])) — a missing column yields
  # character(0), never NULL, so a coercion-based guard silently passes.
  if (!patient_var %in% colnames(cd)) {
    stop("patient_var '", patient_var, "' not found in colData. Available: ",
         paste(colnames(cd), collapse = ", "))
  }
  pv <- as.character(cd[[patient_var]])
  if (anyNA(pv)) {
    stop(sum(is.na(pv)), " sample(s) have a missing '", patient_var,
         "' and cannot be assigned to a patient.")
  }

  cnt <- as.matrix(SummarizedExperiment::assay(tse, "counts"))
  # Empty samples would contribute an all-zero column to the patient mean while still counting
  # in the divisor, leaving that patient's profile summing to < 1. Drop them instead.
  depth <- colSums(cnt)
  if (any(depth <= 0)) {
    warning(sum(depth <= 0), " sample(s) with zero total counts dropped before aggregation: ",
            paste(utils::head(colnames(cnt)[depth <= 0], 5), collapse = ", "),
            if (sum(depth <= 0) > 5) ", ..." else "", call. = FALSE)
    cnt <- cnt[, depth > 0, drop = FALSE]
    cd <- cd[depth > 0, , drop = FALSE]
    pv <- pv[depth > 0]
    if (ncol(cnt) == 0) stop("no samples left after dropping zero-count samples.")
  }

  # per-sample relative abundance, then mean within patient
  rel <- sweep(cnt, 2, colSums(cnt), "/")
  rel_sum <- t(rowsum(t(rel), group = pv))          # taxa x patient (sum of per-sample rel)
  n_per   <- as.integer(table(pv)[colnames(rel_sum)])
  rel_mean <- sweep(rel_sum, 2, n_per, "/")          # mean per-sample rel -> columns sum to 1
  cnt_sum <- t(rowsum(t(cnt), group = pv))[, colnames(rel_mean), drop = FALSE]

  # colData is carried over from the FIRST sample of each patient, which is only meaningful for
  # columns that are constant within a patient. Flag the rest so a caller cannot unknowingly
  # model e.g. Site on the patient-level object.
  varying <- check_between_subject(tse, colnames(cd), patient_var = patient_var,
                                   warn = FALSE)$varying
  if (warn && length(varying) > 0) {
    warning("colData column(s) vary within patient; the first sample's value was kept and is ",
            "not meaningful at the patient level: ", paste(varying, collapse = ", "),
            call. = FALSE)
  }

  pmd <- cd[!duplicated(pv), , drop = FALSE]
  rownames(pmd) <- pv[!duplicated(pv)]
  pmd <- pmd[colnames(rel_mean), , drop = FALSE]

  args <- list(
    assays = list(counts = cnt_sum, relabundance = rel_mean),
    rowData = SummarizedExperiment::rowData(tse),
    colData = pmd
  )
  rt <- tryCatch(TreeSummarizedExperiment::rowTree(tse), error = function(e) NULL)
  if (!is.null(rt)) {
    args$rowTree <- rt
    args$rowNodeLab <- TreeSummarizedExperiment::rowLinks(tse)$nodeLab
  }
  do.call(TreeSummarizedExperiment::TreeSummarizedExperiment, args)
}

#' Identify which variables are constant within patient (between-subject)
#'
#' A between-subject term (e.g. Diagnosis) takes one value per patient; a within-subject term (e.g.
#' Site) varies across a patient's samples. Tests that treat between-subject terms at the sample level
#' are pseudoreplicated, so callers use this to validate a "between" formula or to drop sample-level
#' terms before a patient-level analysis.
#'
#' @param tse A (Tree)SummarizedExperiment.
#' @param vars Character vector of colData column names to classify.
#' @param patient_var Subject identifier column (default "Patient").
#' @param warn Emit a warning naming any within-patient–varying variables (default TRUE).
#'
#' @return list(constant = <vars constant within patient>, varying = <vars that vary within patient>).
#' @export
check_between_subject <- function(tse, vars, patient_var = "Patient", warn = TRUE) {
  cd <- as.data.frame(SummarizedExperiment::colData(tse))
  if (!patient_var %in% colnames(cd)) {
    stop("patient_var '", patient_var, "' not found in colData. Available: ",
         paste(colnames(cd), collapse = ", "))
  }
  pv <- as.character(cd[[patient_var]])
  vars <- setdiff(intersect(vars, colnames(cd)), patient_var)
  varies <- vapply(vars, function(v) {
    n_distinct <- tapply(cd[[v]], pv, function(x) length(unique(x[!is.na(x)])))
    any(n_distinct > 1, na.rm = TRUE)
  }, logical(1))
  varying  <- vars[varies]
  constant <- vars[!varies]
  if (warn && length(varying) > 0) {
    warning("Variable(s) vary within patient (not between-subject): ",
            paste(varying, collapse = ", "), call. = FALSE)
  }
  list(constant = constant, varying = varying)
}
