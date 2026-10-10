#' Prepare MaAsLin3 Data for Plotting
#'
#' Internal helper that normalizes MaAsLin3 input for all plotting functions.
#' Accepts both in-memory result objects from `maaslin3::maaslin3()` and
#' data frames / TSV files matching the `all_results.tsv` structure.
#'
#' @param maaslin3_input Either:
#'   - A list returned by `maaslin3::maaslin3()` (contains `$fit_data_abundance`, etc.)
#'   - A data frame or tibble matching the `all_results.tsv` column structure
#'   - A file path (character) to a TSV file
#' @param coefficient_names Character vector of coefficient names to include
#'   (e.g., `"TypeStent"`, `c("GroupBenign Stricture", "GroupIntra Biliary Malignancy")`).
#'   If NULL, all coefficients are returned.
#' @param model_type One of `"abundance"`, `"prevalence"`, or `"joint"`.
#'   For `"joint"`: uses abundance `coef`/`stderr` but `pval_joint`/`qval_joint`.
#' @param taxa Optional character vector of taxa to filter to.
#'
#' @return A data frame with columns: `taxon`, `coefficient_name`, `coef`,
#'   `null_hypothesis`, `stderr`, `pval`, `qval`, `N`, `N_not_zero`, `error`.
#'   `null_hypothesis` is the value MaAsLin3 tested the coefficient against (0 for
#'   prevalence; the median-comparison reference for abundance); it is 0 when the
#'   source lacks the column.
#'
#' @keywords internal
prepare_maaslin3_data <- function(
  maaslin3_input,
  coefficient_names = NULL,
  model_type = c("abundance", "prevalence", "joint"),
  taxa = NULL
) {
  model_type <- match.arg(model_type)

  # --- Resolve input to a data frame ---
  df <- .resolve_maaslin3_input(maaslin3_input, model_type)

  # --- Filter to requested coefficients ---
  if (!is.null(coefficient_names)) {
    available <- unique(df$name)
    missing <- setdiff(coefficient_names, available)
    if (length(missing) > 0) {
      warning(
        "Coefficient(s) not found: ", paste(missing, collapse = ", "),
        "\nAvailable: ", paste(available, collapse = ", ")
      )
    }
    df <- df[df$name %in% coefficient_names, , drop = FALSE]
  }

  if (nrow(df) == 0) {
    stop("No rows remaining after filtering by coefficient_names.")
  }

  # --- Select p-value / q-value columns based on model_type ---
  if (model_type == "joint") {
    pval_col <- "pval_joint"
    qval_col <- "qval_joint"
  } else {
    pval_col <- "pval_individual"
    qval_col <- "qval_individual"
  }

  for (col in c(pval_col, qval_col)) {
    if (!col %in% names(df)) {
      stop("Required column '", col, "' not found in the data.")
    }
  }

  # --- Build standardized output ---
  # `null_hypothesis` is the value MaAsLin3 tested each coefficient against: 0 for
  # prevalence, and the median-comparison reference (median coefficient for that
  # metadatum across features) for abundance when median_comparison_abundance = TRUE.
  # Older outputs / hand-built TSVs may lack it, so it falls back to 0.
  out <- data.frame(
    taxon            = df$feature,
    coefficient_name = df$name,
    coef             = as.numeric(df$coef),
    null_hypothesis  = if ("null_hypothesis" %in% names(df)) as.numeric(df$null_hypothesis) else 0,
    stderr           = as.numeric(df$stderr),
    pval             = as.numeric(df[[pval_col]]),
    qval             = as.numeric(df[[qval_col]]),
    N                = if ("N" %in% names(df)) as.integer(df$N) else NA_integer_,
    N_not_zero       = if ("N_not_zero" %in% names(df)) as.integer(df$N_not_zero) else NA_integer_,
    error            = if ("error" %in% names(df)) df$error else NA_character_,
    stringsAsFactors = FALSE
  )

  # --- Filter to requested taxa ---
  if (!is.null(taxa)) {
    out <- out[out$taxon %in% taxa, , drop = FALSE]
    if (nrow(out) == 0) {
      warning("No rows remaining after filtering by taxa.")
    }
  }

  # Remove rows with errors if the error column exists and has values
  out <- out[is.na(out$error) | out$error == "" | out$error == "NA", , drop = FALSE]

  rownames(out) <- NULL
  out
}


#' Resolve MaAsLin3 input to a flat data frame
#'
#' @param maaslin3_input The user-provided input (list, data.frame, or file path).
#' @param model_type The model type requested.
#' @return A data frame with at least columns: feature, name, coef, stderr, plus p/q columns.
#' @keywords internal
.resolve_maaslin3_input <- function(maaslin3_input, model_type) {
  # Case 1: file path
  if (is.character(maaslin3_input) && length(maaslin3_input) == 1 && file.exists(maaslin3_input)) {
    df <- utils::read.delim(maaslin3_input, stringsAsFactors = FALSE)
    return(.validate_maaslin3_df(df, model_type))
  }

  # Case 2: data frame / tibble
  if (is.data.frame(maaslin3_input)) {
    return(.validate_maaslin3_df(maaslin3_input, model_type))
  }

  # Case 3: in-memory maaslin3 result object. Joint p/q-values live in the abundance results.
  if (is.list(maaslin3_input)) {
    slot_name <- if (model_type == "joint") "fit_data_abundance" else paste0("fit_data_", model_type)
    if (!slot_name %in% names(maaslin3_input)) {
      stop("'", slot_name, "' not found in maaslin3 result object.")
    }
    return(.validate_maaslin3_df(maaslin3_input[[slot_name]]$results, model_type))
  }

  stop("maaslin3_input must be a data frame, a file path to a TSV, or a maaslin3 result object.")
}


#' Validate that a MaAsLin3 data frame has the required columns
#' @keywords internal
.validate_maaslin3_df <- function(df, model_type) {
  required_base <- c("feature", "name", "coef", "stderr")

  if (model_type == "joint") {
    required <- c(required_base, "pval_joint", "qval_joint")
  } else {
    required <- c(required_base, "pval_individual", "qval_individual")
  }

  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(
      "Missing required columns: ", paste(missing, collapse = ", "),
      "\nAvailable columns: ", paste(names(df), collapse = ", ")
    )
  }

  # For combined TSV files (all_results.tsv), the model column contains
  # "abundance" or "prevalence". For in-memory objects, it contains the
  # model family name (e.g., "linear", "logistic"). Only filter when the
  # column values match the expected TSV format.
  if ("model" %in% names(df)) {
    model_vals <- unique(df$model)
    has_tsv_format <- any(c("abundance", "prevalence") %in% model_vals)

    if (has_tsv_format && model_type != "joint") {
      df <- df[df$model == model_type, , drop = FALSE]
      if (nrow(df) == 0) {
        stop("No rows with model == '", model_type, "' found in the data.")
      }
    } else if (has_tsv_format && model_type == "joint") {
      # For joint, use abundance rows (they contain pval_joint/qval_joint)
      df <- df[df$model == "abundance", , drop = FALSE]
      if (nrow(df) == 0) {
        stop("No abundance model rows found for joint analysis.")
      }
    }
    # If model column has in-memory format (e.g. "linear"), skip filtering —
    # the data was already extracted from the correct slot by .resolve_maaslin3_input()
  }

  df
}
