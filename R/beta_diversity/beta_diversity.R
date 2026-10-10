#' Build a permutation control object for a (possibly blocked) design
#'
#' Wraps `permutations` into a `permute::how()` control blocked by `strata`, so that
#' repeated-measures designs permute only within a subject. Passing the same `strata` to
#' `vegan::adonis2(strata = )` or `stats::anova(cca, strata = )` is equivalent — vegan
#' converts it into `how(blocks = )` internally (`vegan:::getPermuteMatrix`) — but building
#' the control explicitly keeps the design visible and lets callers pass a pre-built `how()`.
#'
#' @param permutations Either a permutation count, or an already-built `permute::how()` object
#'   (returned unchanged, so a caller-supplied design always wins).
#' @param strata Optional blocking vector (e.g. patient ID). `NULL` means free permutation.
#'
#' @return A `permute::how()` object, or `permutations` unchanged when `strata` is `NULL`.
#'
#' @export
build_perm_control <- function(permutations, strata = NULL) {
  # A caller-supplied how() already encodes its own design; overriding it would silently
  # discard their choice.
  if (inherits(permutations, "how")) return(permutations)
  if (is.null(strata)) return(permutations)
  permute::how(blocks = factor(strata), nperm = permutations)
}

#' Tidy a vegan ANOVA-type table (`adonis2()`, `anova.cca()`)
#'
#' Single definition used by every PERMANOVA / dbRDA table in shared. P-values stay numeric;
#' format them at display time.
#'
#' @param tab An `anova.cca` / `adonis2` result.
#' @param keep_terms Optional terms to keep (others dropped before the q-value adjustment).
#' @param p_adjust_method Adjustment across the kept terms (default "BH").
#' @return data.frame with `term`, `df`, `sum_of_sqs`, `r2` (adonis2 only), `f_statistic`,
#'   `p_value`, `q_value`.
#' @export
tidy_vegan_anova <- function(tab, keep_terms = NULL, p_adjust_method = "BH") {
  df <- as.data.frame(tab)
  df$term <- rownames(df)
  df <- df[!df$term %in% c("Residual", "Total"), , drop = FALSE]
  if (!is.null(keep_terms)) df <- df[df$term %in% keep_terms, , drop = FALSE]
  ren <- c(Df = "df", SumOfSqs = "sum_of_sqs", Variance = "sum_of_sqs", R2 = "r2",
           F = "f_statistic", `Pr(>F)` = "p_value")
  hit <- colnames(df) %in% names(ren)
  colnames(df)[hit] <- ren[colnames(df)[hit]]
  df$q_value <- if ("p_value" %in% colnames(df)) stats::p.adjust(df$p_value, method = p_adjust_method) else NA_real_
  rownames(df) <- NULL
  df[, intersect(c("term", "df", "sum_of_sqs", "r2", "f_statistic", "p_value", "q_value"), colnames(df)),
     drop = FALSE]
}

# Kruskal's stress-1 and Shepard-plot data for an ordination of a stored distance.
.shepard <- function(dist_matrix, coords) {
  d0 <- as.matrix(dist_matrix)
  dp <- as.matrix(stats::dist(coords))
  ord <- order(as.vector(d0))
  list(stress = sqrt(sum((dp - d0)^2) / sum(d0^2)),
       shepard_data = data.frame(original = as.vector(d0)[ord], ordination = as.vector(dp)[ord]))
}

# Percent variance explained per PCoA axis, from stored eigenvalues (`<distance>_eig`, written
# by calc_beta_diversity.R) or recomputed with cmdscale(). Negative eigenvalues of
# non-Euclidean distances are kept in the denominator via their absolute value.
.variance_explained <- function(tse, distance_label) {
  md <- S4Vectors::metadata(tse)
  eig <- md[[paste0(distance_label, "_eig")]]
  if (is.null(eig)) {
    eig <- tryCatch(stats::cmdscale(md[[distance_label]], k = 2, eig = TRUE)$eig, error = function(e) {
      warning("Eigenvalue calculation failed: ", e$message)
      NULL
    })
  }
  list(eigenvalues = eig,
       variance_explained = if (!is.null(eig) && sum(abs(eig)) > 0) eig / sum(abs(eig)) * 100)
}

#' Pairwise PERMANOVA using adonis2
#'
#' Performs pairwise adonis2 tests for all group combinations.
#' Tests community composition differences between each pair of groups.
#'
#' @param dist_matrix Distance matrix (class dist)
#' @param groups Factor or character vector of group labels
#' @param strata Optional strata for stratified permutations
#' @param permutations Number of permutations (default: 999)
#' @param p_adjust_method Method for p-value adjustment (default: "BH")
#' @param seed Optional RNG seed set before each adonis2 call, so permutation p-values are
#'   reproducible independently of how much RNG state upstream chunks consumed. `NULL`
#'   (default) preserves the previous behaviour of drawing from the ambient stream.
#'
#' @return data.frame with columns: comparison, df, sum_of_sqs, r2, f_statistic, p_value, q_value
#'
#' @export
pairwise_adonis2 <- function(
  dist_matrix,
  groups,
  strata = NULL,
  permutations = 999,
  p_adjust_method = "BH",
  seed = NULL
) {
  # Input validation
  if (!inherits(dist_matrix, "dist")) {
    stop("dist_matrix must be of class 'dist'")
  }

  groups <- as.character(groups)
  unique_groups <- unique(groups)

  if (length(unique_groups) < 2) {
    warning("Need at least 2 groups for pairwise comparisons")
    return(NULL)
  }

  # Generate all pairwise combinations
  combinations <- utils::combn(unique_groups, 2)
  n_comparisons <- ncol(combinations)

  # Initialize results vectors
  results <- data.frame(
    comparison = character(n_comparisons),
    df = numeric(n_comparisons),
    sum_of_sqs = numeric(n_comparisons),
    r2 = numeric(n_comparisons),
    f_statistic = numeric(n_comparisons),
    p_value = numeric(n_comparisons),
    stringsAsFactors = FALSE
  )

  # Run pairwise tests
  for (i in seq_len(n_comparisons)) {
    group1 <- combinations[1, i]
    group2 <- combinations[2, i]

    # Subset to these two groups
    keep_idx <- groups %in% c(group1, group2)
    dist_subset <- as.dist(as.matrix(dist_matrix)[keep_idx, keep_idx])
    groups_subset <- groups[keep_idx]

    # Prepare strata if provided
    strata_subset <- if (!is.null(strata)) strata[keep_idx] else NULL

    # Create model data
    model_df <- data.frame(Group = factor(groups_subset))

    # Run adonis2 with permutations restricted within `strata` (repeated measures).
    # The blocked design is written out as a permute::how() object rather than passed via
    # adonis2's `strata=`; the two are equivalent (vegan converts `strata=` into
    # how(blocks=) internally, see vegan:::getPermuteMatrix), but stating it here makes the
    # design visible at the call site. Restricting permutations only validates terms that
    # vary WITHIN a block; a term constant within a block cannot be tested this way at all
    # (see beta_permanova_split()).
    perm_ctrl <- build_perm_control(permutations, strata_subset)
    if (!is.null(seed)) set.seed(seed)
    ad_result <- tryCatch({
      vegan::adonis2(
        formula = dist_subset ~ Group,
        data = model_df,
        permutations = perm_ctrl,
        by = "terms"
      )
    }, error = function(e) {
      warning(paste("adonis2 failed for", group1, "vs", group2, ":", e$message))
      NULL
    })

    # Extract results
    results$comparison[i] <- paste(group1, "vs", group2)

    if (!is.null(ad_result)) {
      results$df[i] <- ad_result$Df[1]
      results$sum_of_sqs[i] <- ad_result$SumOfSqs[1]
      results$r2[i] <- ad_result$R2[1]
      results$f_statistic[i] <- ad_result$F[1]
      results$p_value[i] <- ad_result$`Pr(>F)`[1]
    } else {
      results$df[i] <- NA
      results$sum_of_sqs[i] <- NA
      results$r2[i] <- NA
      results$f_statistic[i] <- NA
      results$p_value[i] <- NA
    }
  }

  # Apply multiple testing correction
  results$q_value <- p.adjust(results$p_value, method = p_adjust_method)

  return(results)
}


#' Analyze Beta Diversity
#'
#' Performs PERMANOVA, betadisper, and calculates ordination stress.
#' Returns a list of results that can be passed to plot_beta_diversity().
#'
#' @param tse TreeSummarizedExperiment object
#' @param dimred_name Name of the MDS reduced dimension (e.g., "MDS_bray")
#' @param distance_label Name of the distance matrix in tse metadata (e.g., "bray")
#' @param group_var Column name in colData for the main grouping variable
#' @param formula Right-hand side for PERMANOVA (e.g., "Group + Type + Batch"; interactions
#'   allowed). If NULL, uses group_var only
#' @param strata_var Variable for stratification in PERMANOVA (e.g., "Patient"). Restricts
#'   permutations within each level, which validates *within-subject* terms only — see
#'   `beta_permanova_split()` for between-subject terms.
#' @param permutations Number of permutations for PERMANOVA and betadisper (default: 999)
#' @param p_adjust_method Method for multiple testing correction (default: "BH")
#' @param seed Optional RNG seed set before each permutation test, so results are reproducible
#'   independently of how much RNG state upstream chunks consumed. `NULL` (default) preserves
#'   the previous behaviour of drawing from the ambient stream.
#'
#' @return List with components (always all present; statistical slots are NULL when the
#'   test could not run):
#'   - stress: Numeric, Kruskal's stress-1 value
#'   - shepard_data: data.frame with original and ordination distances
#'   - eigenvalues, variance_explained: PCoA eigenvalues and percent variance per axis
#'   - permanova: data.frame with PERMANOVA results (`tidy_vegan_anova()`)
#'   - betadisper: data.frame with betadisper test results (homogeneity check)
#'   - pairwise_permanova: data.frame with pairwise adonis2 results (if > 2 groups)
#'
#' @export
analyze_beta_diversity <- function(
  tse,
  dimred_name,
  distance_label,
  group_var,
  formula = NULL,
  strata_var = NULL,
  permutations = 999,
  p_adjust_method = "BH",
  seed = NULL
) {
  stopifnot(inherits(tse, "SummarizedExperiment"))
  if (!dimred_name %in% SingleCellExperiment::reducedDimNames(tse)) {
    stop("Reduced dimension ", dimred_name, " not found in tse")
  }
  dist_matrix <- S4Vectors::metadata(tse)[[distance_label]]
  if (is.null(dist_matrix)) stop("Distance matrix ", distance_label, " not found in metadata(tse)")
  cd <- as.data.frame(SummarizedExperiment::colData(tse))
  if (!group_var %in% colnames(cd)) stop("Column ", group_var, " not found in colData")

  # Ordination diagnostics; the statistical slots are filled below when possible, so every
  # return has the same shape.
  out <- c(.shepard(dist_matrix, SingleCellExperiment::reducedDim(tse, dimred_name)),
           .variance_explained(tse, distance_label),
           list(permanova = NULL, betadisper = NULL, pairwise_permanova = NULL))

  # --- PERMANOVA ---
  # The formula is fitted as written (interactions and transformations included); its
  # variables only drive the complete-case filter.
  formula_obj <- stats::as.formula(paste("dist_subset ~", formula %||% group_var))
  model_vars <- unique(c(all.vars(formula_obj)[-1], group_var, strata_var))
  missing_vars <- setdiff(model_vars, colnames(cd))
  if (length(missing_vars) > 0) {
    warning("Variables ", paste(missing_vars, collapse = ", "), " not found in colData")
    return(out)
  }
  metadata_subset <- cd[stats::complete.cases(cd[, model_vars, drop = FALSE]), model_vars, drop = FALSE]
  if (nrow(metadata_subset) < 5) {
    warning("Not enough complete cases for PERMANOVA")
    return(out)
  }
  keep <- rownames(metadata_subset)
  dist_subset <- stats::as.dist(as.matrix(dist_matrix)[keep, keep])
  strata <- if (!is.null(strata_var)) metadata_subset[[strata_var]] else NULL

  # Permutations restricted within `strata_var` (written as an explicit permute::how(); vegan
  # converts adonis2's `strata=` into the same object). This validly tests *within-subject*
  # terms only: a term constant within a block is identical in every permutation, so its
  # p-value is meaningless. Test between-subject terms with beta_permanova_split().
  if (!is.null(seed)) set.seed(seed)
  permanova <- tryCatch(
    vegan::adonis2(formula_obj, data = metadata_subset,
                   permutations = build_perm_control(permutations, strata), by = "terms"),
    error = function(e) {
      warning("PERMANOVA failed: ", e$message)
      NULL
    })
  if (!is.null(permanova)) out$permanova <- tidy_vegan_anova(permanova, p_adjust_method = p_adjust_method)

  # --- Betadisper (homogeneity of dispersion) ---
  groups <- metadata_subset[[group_var]]
  out$betadisper <- tryCatch({
    disp <- vegan::betadisper(dist_subset, groups, type = "centroid")
    if (!is.null(seed)) set.seed(seed)
    disp_test <- vegan::permutest(disp, permutations = permutations)
    p <- disp_test$tab$`Pr(>F)`[1]
    data.frame(
      f_statistic = disp_test$statistic,
      p_value = p,
      interpretation = ifelse(p < 0.05, "Heterogeneous (caution: PERMANOVA may be confounded)",
                              "Homogeneous (valid)")
    )
  }, error = function(e) {
    warning("Betadisper failed: ", e$message)
    NULL
  })

  # --- Pairwise PERMANOVA (only with > 2 groups) ---
  if (length(unique(groups)) > 2) {
    out$pairwise_permanova <- tryCatch(
      pairwise_adonis2(dist_subset, groups, strata = strata, permutations = permutations,
                       p_adjust_method = p_adjust_method, seed = seed),
      error = function(e) {
        warning("Pairwise PERMANOVA failed: ", e$message)
        NULL
      })
  }

  out
}


#' Plot Beta Diversity Ordination
#'
#' Creates MDS ordination plot with optional Shepard stress diagram or variance explained barplot.
#'
#' @param tse TreeSummarizedExperiment object
#' @param dimred_name Name of the MDS reduced dimension (e.g., "MDS_bray")
#' @param color_by Column name in colData to color points by
#' @param shape_by Column name in colData to shape points by (optional)
#' @param color_label Label for color legend
#' @param shape_label Label for shape legend
#' @param distance_label Name of distance matrix in tse metadata (required if inset_type="shepard" and no analysis_result)
#' @param analysis_result Pre-computed results from analyze_beta_diversity() (optional)
#' @param inset_type Type of inset plot: "variance" (eigenvalue barplot), "shepard" (stress diagram), or "none"
#' @param annotate_stats Logical, add PERMANOVA annotation to plot (default: FALSE)
#' @param color_palette RColorBrewer palette name (default: "Set2")
#' @param shape_values Named vector of shape values for different categories (optional)
#' @param main_height Relative height for main plot in patchwork (default: 7)
#' @param inset_height Relative height for inset plot in patchwork (default: 1)
#' @param title Optional title for the main plot
#' @param label_by Column name in colData to use for point labels (optional)
#'
#' @return ggplot or patchwork object
#'
#' @export
plot_beta_diversity <- function(
  tse,
  dimred_name,
  color_by,
  shape_by = NULL,
  color_label = NULL,
  shape_label = NULL,
  distance_label = NULL,
  analysis_result = NULL,
  inset_type = c("variance", "shepard", "none"),
  annotate_stats = FALSE,
  color_palette = "Set2",
  shape_values = NULL,
  main_height = 7,
  inset_height = 1,
  title = NULL,
  label_by = NULL
) {
  # Match inset_type argument

  inset_type <- match.arg(inset_type)
  # Validate inputs
  if (!inherits(tse, "TreeSummarizedExperiment") && !inherits(tse, "SummarizedExperiment")) {
    stop("tse must be a TreeSummarizedExperiment or SummarizedExperiment object")
  }

  if (!(dimred_name %in% SingleCellExperiment::reducedDimNames(tse))) {
    stop(paste("Reduced dimension", dimred_name, "not found in tse"))
  }

  if (!(color_by %in% colnames(SummarizedExperiment::colData(tse)))) {
    stop(paste("Column", color_by, "not found in colData"))
  }

  if (!is.null(shape_by) && !(shape_by %in% colnames(SummarizedExperiment::colData(tse)))) {
    stop(paste("Column", shape_by, "not found in colData"))
  }

  # Set default labels if not provided
  if (is.null(color_label)) color_label <- color_by
  if (is.null(shape_label) && !is.null(shape_by)) shape_label <- shape_by

  # Determine title if not provided
  if (is.null(title)) {
    distance_name <- dplyr::case_when(
      stringr::str_detect(dimred_name, "jaccard") ~ "Jaccard",
      stringr::str_detect(dimred_name, "bray") ~ "Bray-Curtis",
      stringr::str_detect(dimred_name, "unifrac_weighted") ~ "Weighted UniFrac",
      stringr::str_detect(dimred_name, "unifrac") ~ "Unweighted UniFrac",
      stringr::str_detect(dimred_name, "rclr") ~ "(approx.) Aitchison with Robust CLR",
      TRUE ~ dimred_name
    )
    title <- paste("MDS Ordination:", distance_name)
  }

  # Extract reduced dimensions data
  dim_data <- SingleCellExperiment::reducedDim(tse, dimred_name)

  # Extract relevant colData
  coldata_vars <- c(color_by)
  if (!is.null(shape_by)) coldata_vars <- c(coldata_vars, shape_by)
  if (!is.null(label_by)) coldata_vars <- c(coldata_vars, label_by)
  coldata_subset <- SummarizedExperiment::colData(tse)[, coldata_vars, drop = FALSE]

  # Combine
  plot_df <- cbind(as.data.frame(dim_data), as.data.frame(coldata_subset))
  colnames(plot_df)[1:2] <- c("Dim1", "Dim2")

  # Create base MDS plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Dim1, y = Dim2)) +
    # Add ellipses first (behind points)
    ggplot2::stat_ellipse(
      ggplot2::aes(color = !!rlang::sym(color_by), fill = !!rlang::sym(color_by)),
      geom = "polygon",
      alpha = 0.1,
      level = 0.95,
      linewidth = 0.5
    ) +
    # Add points
    {
      if (!is.null(shape_by)) {
        ggplot2::geom_point(ggplot2::aes(color = !!rlang::sym(color_by), shape = !!rlang::sym(shape_by)), size = 2.5, alpha = 0.8)
      } else {
        ggplot2::geom_point(ggplot2::aes(color = !!rlang::sym(color_by)), size = 2.5, alpha = 0.8)
      }
    } +
    # Add labels if requested
    {
      if (!is.null(label_by)) {
        ggrepel::geom_text_repel(
          ggplot2::aes(label = !!rlang::sym(label_by), color = !!rlang::sym(color_by)),
          size = 2.5,
          max.overlaps = 20,
          box.padding = 0.3,
          point.padding = 0.3,
          segment.size = 0.2,
          segment.alpha = 0.6,
          show.legend = FALSE
        )
      }
    } +
    # Define scales
    ggplot2::scale_color_brewer(palette = color_palette, name = color_label) +
    ggplot2::scale_fill_brewer(palette = color_palette, name = color_label, guide = "none") +
    # Add shape scale if needed
    {
      if (!is.null(shape_by)) {
        if (!is.null(shape_values)) {
          ggplot2::scale_shape_manual(values = shape_values, name = shape_label)
        } else {
          ggplot2::scale_shape_discrete(name = shape_label)
        }
      }
    }

  # Determine axis labels with variance explained if available
  if (!is.null(analysis_result) && !is.null(analysis_result$variance_explained)) {
    var_exp <- analysis_result$variance_explained
    x_label <- sprintf("PCoA 1 (%.1f%%)", var_exp[1])
    y_label <- sprintf("PCoA 2 (%.1f%%)", var_exp[2])
  } else {
    x_label <- "PCoA 1"
    y_label <- "PCoA 2"
  }

  p <- p +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      title = title
    )

  # Add PERMANOVA annotation if requested
  if (annotate_stats && !is.null(analysis_result) && !is.null(analysis_result$permanova)) {
    first_term <- analysis_result$permanova$term[1]
    r2 <- analysis_result$permanova$r2[1]
    p_val <- analysis_result$permanova$p_value[1]
    p_val_formatted <- ifelse(p_val < 0.001, "< 0.001", format(round(p_val, 3), nsmall = 3))

    significance_text <- paste0(
      "PERMANOVA (", first_term, ")\n",
      "R² = ", format(round(r2, 3), nsmall = 3),
      ", p=", p_val_formatted
    )

    p <- p + ggplot2::annotate(
      "label",
      x = Inf, y = -Inf,
      label = significance_text,
      hjust = 1.1, vjust = -0.5,
      size = 2.8,
      color = "black",
      fill = ggplot2::alpha("white", 0.7),
      border.colour = NA,
      label.padding = ggplot2::unit(0.25, "lines")
    )
  }

  # If no inset requested, return main plot only
  if (inset_type == "none") {
    return(p)
  }

  # Create inset plot based on type
  inset_plot <- NULL

  if (inset_type == "variance") {
    # Create variance explained barplot
    if (!is.null(analysis_result) && !is.null(analysis_result$variance_explained)) {
      var_exp <- analysis_result$variance_explained
      n_axes <- min(10, length(var_exp))

      var_df <- data.frame(
        axis = factor(paste0("PC", 1:n_axes), levels = paste0("PC", 1:n_axes)),
        variance = var_exp[1:n_axes]
      )

      inset_plot <- ggplot2::ggplot(var_df, ggplot2::aes(x = axis, y = variance)) +
        ggplot2::geom_col(fill = "steelblue", alpha = 0.8) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", variance)),
                  vjust = -0.3, size = 2.5) +
        ggplot2::labs(
          x = NULL,
          y = "Variance (%)",
          subtitle = sprintf("PC1+PC2: %.1f%%", sum(var_exp[1:2]))
        )
    } else {
      warning("Variance explained data not available. Returning plot without inset.")
      return(p)
    }
  } else if (inset_type == "shepard") {
    # Create Shepard plot
    if (!is.null(analysis_result)) {
      shepard_df <- analysis_result$shepard_data
      stress <- analysis_result$stress
    } else {
      # Compute Shepard data on the fly
      if (is.null(distance_label)) {
        warning("distance_label required for Shepard plot when analysis_result not provided. Returning plot without inset.")
        return(p)
      }

      if (!(distance_label %in% names(tse@metadata))) {
        warning(paste("Distance matrix", distance_label, "not found in tse@metadata. Returning plot without inset."))
        return(p)
      }

      sh <- .shepard(tse@metadata[[distance_label]], SingleCellExperiment::reducedDim(tse, dimred_name))
      shepard_df <- sh$shepard_data
      stress <- sh$stress
    }

    inset_plot <- ggplot2::ggplot(shepard_df, ggplot2::aes(x = original, y = ordination)) +
      ggplot2::geom_smooth(method = "loess", se = TRUE, color = "steelblue") +
      ggplot2::geom_point(size = 0.5, alpha = 0.06) +
      ggplot2::labs(
        x = "Original distance",
        y = "MDS distance",
        subtitle = paste("Stress:", round(stress, 3))
      )
  }

  # Combine plots with patchwork
  combined <- p + inset_plot +
    patchwork::plot_layout(heights = c(main_height, inset_height), ncol = 1)

  return(combined)
}
