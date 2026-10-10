#' Plot Alpha Diversity with LMM Significance Testing
#'
#' Generates a boxplot for an alpha diversity index, overlaid with jittered points,
#' and annotates the model result. Inference comes from `calc_alpha_diversity_summary()`
#' (same model, same Satterthwaite post-hoc, same adjustment), so the brackets on the figure
#' match the post-hoc table.
#'
#' Requires `calc_alpha_diversity_summary()` and `describe_groups()` / `classify_significance()`.
#'
#' @param data A data frame containing the alpha diversity values and metadata (e.g., from colData).
#' @param x_var Character string, the name of the variable for the x-axis (grouping variable).
#' @param y_var Character string, the name of the alpha diversity index column.
#' @param color_var Character string, the name of the variable to map to color aesthetic. If NULL, uses x_var.
#' @param shape_var Character string, the name of the variable to map to shape aesthetic. If NULL, no shapes are used.
#' @param connect_subjects Logical, whether to draw dashed lines connecting points from the same subject
#'   (default: FALSE). Uses the random effect variable from the LMM formula; points are not jittered.
#' @param x_label Character string for the x-axis label.
#' @param y_label Character string for the y-axis label.
#' @param color_label Character string for the color legend title. Defaults to x_label when color_var is
#'   NULL, else color_var.
#' @param shape_label Character string for the shape legend title. Ignored if shape_var is NULL.
#' @param palette Character string, name of the RColorBrewer palette for color/fill (default: "Set2").
#' @param lmm_formula_str Character string defining the LMM formula (e.g., "Shannon ~ Group2 + Sex + (1|Patient)").
#' @param ref_group Reference level: brackets are drawn only for comparisons involving it. If NULL, all
#'   significant pairwise comparisons are drawn.
#' @param post_hoc_method P-value adjustment method across all pairwise comparisons (default: "tukey").
#' @param annotation_term ANOVA term shown in the top-right label. If NULL or "auto": the first
#'   interaction term if present, else x_var, else all significant terms.
#' @param show_posthoc_label Logical, whether to show the post-hoc family label at bottom-right (default: TRUE).
#' @param show_lmm_label Logical, whether to show the top-right ANOVA/Kruskal-Wallis label (default: TRUE).
#' @param test "anova" (default; linear (mixed) model, emmeans post-hoc) or "kruskal" (Kruskal-Wallis,
#'   Dunn's post-hoc; ignores random effects).
#' @param posthoc Logical, draw significance brackets (default: TRUE).
#' @param style Named list of cosmetic parameters from `alpha_plot_style()`.
#'
#' @return A ggplot object.
plot_alpha_diversity_summary_lmm <- function(data, x_var, y_var, color_var = NULL, shape_var = NULL,
                                             x_label, y_label, color_label = NULL, shape_label = NULL,
                                             palette = "Set2", connect_subjects = FALSE,
                                             lmm_formula_str,
                                             ref_group = NULL,
                                             post_hoc_method = "tukey",
                                             annotation_term = NULL,
                                             show_posthoc_label = TRUE,
                                             show_lmm_label = TRUE,
                                             test = "anova",
                                             posthoc = TRUE,
                                             style = alpha_plot_style()) {
  if (is.null(color_var)) {
    color_var <- x_var
    if (is.null(color_label)) color_label <- x_label
  } else if (is.null(color_label)) {
    color_label <- color_var
  }

  random_effect_var <- extract_random_effect_var(lmm_formula_str)
  if (connect_subjects && is.null(random_effect_var)) {
    stop("connect_subjects=TRUE but no random effect found in LMM formula. ",
         "Please include a random effect like '(1|Subject)' in your formula.")
  }

  # --- Complete cases over every modelled and plotted variable ---
  needed <- unique(c(all.vars(stats::as.formula(lmm_formula_str)),
                     y_var, x_var, color_var, shape_var, random_effect_var))
  missing_vars <- setdiff(needed, colnames(data))
  if (length(missing_vars) > 0) {
    stop("The following variables are missing from the input data: ", paste(missing_vars, collapse = ", "))
  }
  plot_data <- data[stats::complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
  if (nrow(plot_data) < nrow(data)) {
    message(sprintf("Filtered data: %d complete cases out of %d original rows (%.1f%% retained)",
                    nrow(plot_data), nrow(data), 100 * nrow(plot_data) / nrow(data)))
  }

  # --- Inference (shared with the tables) ---
  stats_res <- calc_alpha_diversity_summary(
    plot_data, y_var = y_var, group_var = x_var, lmm_formula_str = lmm_formula_str,
    ref_group = ref_group, test = test, post_hoc_method = post_hoc_method,
    statistics = "mean_sd", posthoc = posthoc
  )
  if (is.null(stats_res)) {
    stop("Insufficient complete cases (", nrow(plot_data), ") for the model on ", y_var,
         ". Need at least 10 complete observations.")
  }
  contrasts <- stats_res$posthoc
  if (!is.null(contrasts) && !is.null(ref_group) && !identical(test, "kruskal")) {
    contrasts <- contrasts[contrast_involves(contrasts$comparison, ref_group), , drop = FALSE]
  }

  # --- Base plot: subject lines, points, boxes ---
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]],
                                               color = .data[[color_var]]))

  if (connect_subjects) {
    p <- p + ggplot2::geom_line(ggplot2::aes(group = .data[[random_effect_var]]),
                                linetype = "dashed", alpha = style$connector_alpha,
                                linewidth = style$connector_linewidth,
                                color = style$connector_colour, na.rm = TRUE)
  }

  point_aes <- if (!is.null(shape_var)) ggplot2::aes(shape = .data[[shape_var]]) else ggplot2::aes()
  p <- p + if (connect_subjects) {
    ggplot2::geom_point(point_aes, alpha = style$point_alpha, size = style$point_size, na.rm = TRUE)
  } else {
    ggplot2::geom_jitter(point_aes, width = style$jitter_width, height = 0,
                         alpha = style$point_alpha, size = style$point_size, na.rm = TRUE)
  }

  # Boxes on top of the points (shape is mapped on the point layer only, so boxes ignore it).
  p <- p + ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[color_var]]),
                                 width = style$box_width, color = style$box_colour, outlier.shape = NA,
                                 alpha = style$box_alpha, linewidth = style$box_linewidth, na.rm = TRUE)

  # --- Annotations ---
  if (show_lmm_label) {
    p <- add_corner_label(p, model_label_text(stats_res$model_anova, x_var, test, annotation_term),
                          y = Inf, vjust = 1.5, style = style)
  }
  n_brackets <- 0
  if (posthoc && !is.null(contrasts) && nrow(contrasts) > 0) {
    bracket_res <- add_significance_brackets(p, contrasts, plot_data, x_var, y_var, style)
    p <- bracket_res$plot
    n_brackets <- bracket_res$n
  }
  if (show_posthoc_label && posthoc) {
    family <- if (!is.null(ref_group)) paste("Post-hoc vs", ref_group) else "All pairwise"
    status <- if (is.null(contrasts)) "not available" else if (n_brackets == 0) "none significant" else NULL
    label <- if (is.null(status)) paste0(family, " (", post_hoc_method, "-adj. p)") else paste0(family, ": ", status)
    p <- add_corner_label(p, label, y = -Inf, vjust = -0.5, style = style,
                          colour = if (is.null(status)) "black" else "grey50")
  }

  # --- Scales and theme ---
  p <- p +
    ggplot2::scale_y_continuous(expand = style$y_expand) +
    ggplot2::scale_color_brewer(palette = palette, name = color_label) +
    ggplot2::scale_fill_brewer(palette = palette, name = color_label) +
    ggplot2::labs(x = x_label, y = y_label) +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 10 + n_brackets * 10, r = 5, b = 5, l = 5,
                                                 unit = "pt")) +
    ggplot2::coord_cartesian(clip = "off")

  if (!is.null(shape_var) && !is.null(shape_label)) {
    p <- p + ggplot2::scale_shape_discrete(name = shape_label)
  }
  p
}

# Format p-values consistently
format_p_value <- function(p_val) {
  ifelse(p_val < 0.001, "p<0.001", paste0("p=", format(round(p_val, 3), nsmall = 3)))
}

# Grouping variable of the first random effect, e.g. "Patient" from "(1|Patient)"; NULL if none.
extract_random_effect_var <- function(lmm_formula_str) {
  m <- regmatches(lmm_formula_str, regexpr("\\([^|]*\\|\\s*([^)]+)\\)", lmm_formula_str))
  if (length(m) == 0) return(NULL)
  trimws(gsub("\\([^|]*\\|\\s*([^)]+)\\)", "\\1", m))
}

# Text of the top-right model label.
model_label_text <- function(model_anova, x_var, test, annotation_term = NULL) {
  if (is.null(model_anova)) return(if (identical(test, "kruskal")) "Kruskal-Wallis failed" else "ANOVA failed")
  if (identical(test, "kruskal")) {
    return(paste0("Kruskal-Wallis ", format_p_value(model_anova$p_value[1])))
  }
  p_of <- function(term) model_anova$p_value[model_anova$term == term]
  if (!is.null(annotation_term) && annotation_term != "auto") {
    if (!annotation_term %in% model_anova$term) return(paste("ANOVA", annotation_term, "not found"))
    return(paste("ANOVA", annotation_term, format_p_value(p_of(annotation_term))))
  }
  interaction <- grep(":", model_anova$term, value = TRUE, fixed = TRUE)
  if (length(interaction) > 0) {
    return(paste("ANOVA", interaction[1], format_p_value(p_of(interaction[1]))))
  }
  if (x_var %in% model_anova$term) {
    return(paste("ANOVA", x_var, format_p_value(p_of(x_var))))
  }
  sig <- model_anova[!is.na(model_anova$p_value) & model_anova$p_value < 0.05, , drop = FALSE]
  if (nrow(sig) == 0) return("ANOVA: no significant terms")
  paste0("ANOVA: ", paste0(sig$term, " (", format_p_value(sig$p_value), ")", collapse = ", "))
}

add_corner_label <- function(p, label, y, vjust, style, colour = "black") {
  p + ggplot2::annotate("label", x = Inf, y = y, label = label,
                        hjust = 1.1, vjust = vjust, size = style$label_size, color = colour,
                        fill = ggplot2::alpha("white", 0.7), border.colour = NA,
                        label.padding = ggplot2::unit(0.25, "lines"))
}

# Draw brackets for significant contrasts (columns `comparison`, `p_adj`). Returns the plot and
# the number of brackets drawn.
add_significance_brackets <- function(p, contrasts, plot_data, x_var, y_var, style = alpha_plot_style()) {
  x_levels <- levels(factor(plot_data[[x_var]]))
  g <- split_contrast(contrasts$comparison)

  # emmeans may prefix levels with the variable name (e.g. "Batch2023" for level "2023"):
  # fall back to the level the label ends with.
  to_position <- function(labels) {
    vapply(labels, function(l) {
      hit <- match(l, x_levels)
      if (is.na(hit)) hit <- which(endsWith(l, x_levels))[1]
      as.integer(hit)
    }, integer(1))
  }
  stat_data <- data.frame(xmin = to_position(g$group1), xmax = to_position(g$group2),
                          label = classify_significance(contrasts$p_adj))
  stat_data <- stat_data[!is.na(stat_data$xmin) & !is.na(stat_data$xmax) & stat_data$label != "ns", ,
                         drop = FALSE]
  if (nrow(stat_data) == 0) return(list(plot = p, n = 0))

  max_y <- max(plot_data[[y_var]], na.rm = TRUE)
  y_range <- max_y - min(plot_data[[y_var]], na.rm = TRUE)
  step <- if (y_range > 0) y_range * style$sig_step else 0.1
  tip <- y_range * style$sig_tip
  stat_data$y <- max_y + step * seq_len(nrow(stat_data))

  # Drawn with annotate() rather than ggpubr to avoid version-dependent behaviour.
  for (i in seq_len(nrow(stat_data))) {
    r <- stat_data[i, ]
    p <- p +
      ggplot2::annotate("segment", x = r$xmin, xend = r$xmax, y = r$y, yend = r$y, linewidth = 0.4) +
      ggplot2::annotate("segment", x = c(r$xmin, r$xmax), xend = c(r$xmin, r$xmax),
                        y = r$y, yend = r$y - tip, linewidth = 0.4) +
      ggplot2::annotate("text", x = (r$xmin + r$xmax) / 2, y = r$y, label = r$label,
                        size = style$sig_size, vjust = -0.3)
  }
  list(plot = p, n = nrow(stat_data))
}
