#' dbRDA (Bray-Curtis) ordination with constraint arrows labelled by PERMANOVA R²
#'
#' Factor constraints are drawn to their level centroids, numeric constraints as biplot
#' arrows. Each arrow is labelled with its term's sequential (by = "terms") PERMANOVA R² and
#' significance stars, so the label depends on the order of terms in the formula.
#'
#' Requires `tidy_vegan_anova()` (beta_diversity.R) and `classify_significance()`.
#'
#' @param tse TSE with a `relabundance` assay.
#' @param rda_formula_str Constraint formula, with or without the leading "~".
#' @param color_by,shape_by colData columns for point colour / shape (shape optional).
#' @param color_label,shape_label Legend titles (unused by the default scales; kept for callers).
#' @param arrow_scaling Arrow length relative to the largest site score (default 0.7).
#' @param strata_var Optional colData column restricting PERMANOVA permutations.
#' @param permutations Number of PERMANOVA permutations for the labels (default 999).
#' @return A ggplot object (a titled empty plot when the model cannot be fitted).
plot_rda_summary <- function(tse,
                             rda_formula_str,
                             color_by, shape_by = NULL,
                             color_label, shape_label,
                             arrow_scaling = 0.7,
                             strata_var = NULL,
                             permutations = 999) {
  error_plot <- function(subtitle) {
    ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::labs(title = "RDA Failed", subtitle = subtitle) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9, color = "grey40"))
  }

  metadata <- as.data.frame(SummarizedExperiment::colData(tse))
  rda_formula_str <- trimws(rda_formula_str)
  if (!startsWith(rda_formula_str, "~")) rda_formula_str <- paste("~", rda_formula_str)
  all_variables <- all.vars(stats::as.formula(rda_formula_str))
  if (!is.null(strata_var) && !strata_var %in% colnames(metadata)) {
    warning("Strata variable not found. Proceeding without strata.")
    strata_var <- NULL
  }

  # Complete cases over the constraints and strata
  complete <- stats::complete.cases(metadata[, unique(c(all_variables, strata_var)), drop = FALSE])
  if (!all(complete)) warning("Removed ", sum(!complete), " samples due to missing values")
  model_data <- metadata[complete, , drop = FALSE]
  rel <- SummarizedExperiment::assay(tse, "relabundance")[, complete, drop = FALSE]

  # Drop constant constraints
  constraining <- setdiff(all_variables, strata_var)
  is_constant <- vapply(constraining, function(v) {
    x <- model_data[[v]]
    if (is.numeric(x)) stats::var(x) < 1e-10 else length(unique(x)) <= 1
  }, logical(1))
  if (any(is_constant)) {
    warning("Excluding zero variance variables: ", paste(constraining[is_constant], collapse = ", "))
  }
  valid_vars <- constraining[!is_constant]
  if (length(valid_vars) == 0) return(error_plot("No valid constraining variables"))
  if (nrow(model_data) < max(5, length(unique(model_data[[valid_vars[1]]])) + 2)) {
    return(error_plot("Not enough samples"))
  }

  # --- Fit ---
  bray_dist <- vegan::vegdist(t(rel), method = "bray")
  fml <- stats::as.formula(paste("bray_dist ~", paste(valid_vars, collapse = " + ")))
  rda_object <- tryCatch(vegan::dbrda(fml, data = model_data, na.action = stats::na.fail),
                         error = function(e) {
                           warning("dbRDA failed: ", e$message)
                           NULL
                         })
  if (is.null(rda_object)) return(error_plot("dbRDA calculation error"))

  permanova <- tryCatch(
    tidy_vegan_anova(vegan::adonis2(fml, data = model_data, permutations = permutations, by = "terms",
                                    strata = if (!is.null(strata_var)) model_data[[strata_var]])),
    error = function(e) {
      warning("PERMANOVA failed: ", e$message)
      NULL
    })

  sample_scores <- vegan::scores(rda_object, display = "sites", choices = 1:2)
  importance <- summary(rda_object)$cont$importance
  var_expl <- if (!is.null(importance)) importance["Proportion Explained", c("dbRDA1", "dbRDA2")] * 100 else c(NA, NA)

  # --- Arrows: centroids for factor levels, biplot arrows for numeric constraints ---
  # term_of maps each score label to its term: "<var><level>" for factors, "<var>" for numerics.
  term_of <- unlist(lapply(valid_vars, function(v) {
    x <- model_data[[v]]
    labs <- if (is.numeric(x)) v else paste0(v, unique(as.character(x)))
    stats::setNames(rep(v, length(labs)), labs)
  }))
  numeric_vars <- valid_vars[vapply(valid_vars, function(v) is.numeric(model_data[[v]]), logical(1))]
  as_arrows <- function(sc) {
    if (is.null(sc) || nrow(sc) == 0) return(NULL)
    data.frame(label = rownames(sc), dbRDA1 = sc[, 1], dbRDA2 = sc[, 2])
  }
  cn <- as_arrows(tryCatch(vegan::scores(rda_object, display = "cn", choices = 1:2), error = function(e) NULL))
  bp <- as_arrows(tryCatch(vegan::scores(rda_object, display = "bp", choices = 1:2), error = function(e) NULL))
  if (!is.null(bp)) bp <- bp[unname(term_of[bp$label]) %in% numeric_vars, , drop = FALSE]
  arrows <- rbind(cn, bp)
  if (!is.null(arrows)) {
    arrows <- arrows[abs(arrows$dbRDA1) > 1e-10 | abs(arrows$dbRDA2) > 1e-10, , drop = FALSE]
  }

  if (!is.null(arrows) && nrow(arrows) > 0) {
    scale <- arrow_scaling * max(abs(sample_scores)) / max(abs(c(arrows$dbRDA1, arrows$dbRDA2)))
    arrows$dbRDA1 <- arrows$dbRDA1 * scale
    arrows$dbRDA2 <- arrows$dbRDA2 * scale
    term <- unname(term_of[arrows$label])
    hit <- if (!is.null(permanova)) match(term, permanova$term) else rep(NA_integer_, nrow(arrows))
    r2 <- permanova$r2[hit]
    # Same thresholds as everywhere else; a non-significant arrow gets no suffix rather than "ns".
    stars <- classify_significance(permanova$p_value[hit])
    stars[stars == "ns"] <- ""
    arrows$label_with_significance <- ifelse(
      is.na(r2), paste0(arrows$label, " (R² NA)"),
      paste0(arrows$label, stars, "\n(R²=", format(round(r2, 3), nsmall = 3), ")"))
  }

  # --- Plot ---
  plot_vars <- c(color_by, shape_by)
  plot_df <- cbind(as.data.frame(sample_scores), model_data[, plot_vars, drop = FALSE])
  for (v in plot_vars) {
    if (anyNA(plot_df[[v]])) {
      plot_df[[v]] <- forcats::fct_na_value_to_level(factor(plot_df[[v]]), "(NA)")
    }
  }

  point_aes <- if (!is.null(shape_by)) {
    ggplot2::aes(color = .data[[color_by]], shape = .data[[shape_by]])
  } else {
    ggplot2::aes(color = .data[[color_by]])
  }
  axis_label <- function(i) {
    if (is.na(var_expl[i])) paste0("dbRDA", i) else sprintf("dbRDA%d (%.1f%%)", i, var_expl[i])
  }
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = dbRDA1, y = dbRDA2)) +
    ggplot2::geom_point(point_aes, size = 2.5, alpha = 0.7) +
    ggplot2::labs(x = axis_label(1), y = axis_label(2))

  if (!is.null(arrows) && nrow(arrows) > 0) {
    p <- p +
      ggplot2::geom_segment(
        data = arrows, ggplot2::aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
        arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"), type = "closed"),
        color = "black", linewidth = 0.5, inherit.aes = FALSE
      ) +
      ggrepel::geom_text_repel(
        data = arrows, ggplot2::aes(x = dbRDA1, y = dbRDA2, label = label_with_significance),
        color = "black", size = 2.8,
        box.padding = ggplot2::unit(0.35, "lines"), point.padding = ggplot2::unit(0.5, "lines"),
        segment.color = "grey50", max.overlaps = Inf, force = 2, inherit.aes = FALSE
      )
  }
  p
}
