#' MaAsLin3 Volcano Plot
#'
#' Creates a volcano plot for MaAsLin3 differential abundance results.
#'
#' "Up"/"Down" are judged against the value MaAsLin3 tested each coefficient against
#' (`null_hypothesis`: the median coefficient for abundance with the default median comparison,
#' 0 for prevalence), not against 0; a dashed line marks that null when it is not 0.
#'
#' Requires `prepare_maaslin3_data()` (prepare_maaslin3_data.R).
#'
#' @param maaslin3_input MaAsLin3 input: result object, data frame, or TSV path.
#'   See `prepare_maaslin3_data()` for details.
#' @param coefficient_name Single character string specifying the coefficient
#'   (e.g., `"TypeStent"`, `"GroupBenign Stricture"`).
#' @param model_type One of `"abundance"`, `"prevalence"`, or `"joint"` (default: `"abundance"`).
#' @param q_threshold Numeric threshold for q-value significance (default: 0.05).
#' @param coef_threshold Minimum absolute distance of the coefficient from its null (default: 0).
#' @param label_significant Logical, whether to label significant taxa (default: TRUE).
#' @param max_labels Maximum number of labels to show (default: 20).
#' @param max_overlaps Passed to `ggrepel::geom_text_repel` (default: 25).
#' @param show_se Logical, whether to show SE as horizontal error bars (default: FALSE).
#' @param color_palette Named vector with colors for "Significant Up", "Significant Down",
#'   and "Not Significant".
#' @param title Optional plot title.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' plot_maaslin3_volcano(maaslin3_result, "TypeStent")
#' plot_maaslin3_volcano("all_results.tsv", "TypeStent", model_type = "joint")
#' }
#'
#' @export
plot_maaslin3_volcano <- function(
  maaslin3_input,
  coefficient_name,
  model_type = "abundance",
  q_threshold = 0.05,
  coef_threshold = 0,
  label_significant = TRUE,
  max_labels = 20,
  max_overlaps = 25,
  show_se = FALSE,
  color_palette = NULL,
  italicize_taxa = TRUE,
  title = NULL
) {

  # Parameter validation
  if (!is.character(coefficient_name) || length(coefficient_name) != 1) {
    stop("coefficient_name must be a single character string.")
  }
  if (!is.numeric(q_threshold) || q_threshold <= 0 || q_threshold > 1) {
    stop("q_threshold must be a numeric value between 0 and 1.")
  }
  if (!is.numeric(coef_threshold) || coef_threshold < 0) {
    stop("coef_threshold must be a non-negative numeric value.")
  }

  # Prepare data
  plot_data <- prepare_maaslin3_data(
    maaslin3_input,
    coefficient_names = coefficient_name,
    model_type = model_type
  )

  if (nrow(plot_data) == 0) {
    stop("No data available for the specified coefficient and model type.")
  }

  # Default color palette (same as ANCOM-BC2)
  if (is.null(color_palette)) {
    color_palette <- c(
      "Significant Down" = "#377eb8",
      "Significant Up"   = "#e41a1c",
      "Not Significant"  = "grey50"
    )
  }

  required_colors <- c("Significant Up", "Significant Down", "Not Significant")
  if (!all(required_colors %in% names(color_palette))) {
    missing_colors <- setdiff(required_colors, names(color_palette))
    stop("Color palette missing: ", paste(missing_colors, collapse = ", "))
  }

  # Compute plotting variables
  plot_data$neg_log10_q <- -log10(plot_data$qval)

  # Handle infinite values (q = 0)
  finite_vals <- plot_data$neg_log10_q[is.finite(plot_data$neg_log10_q)]
  if (length(finite_vals) > 0) {
    inf_replacement <- max(finite_vals, na.rm = TRUE) * 1.1
  } else {
    inf_replacement <- 10
  }
  plot_data$neg_log10_q[is.infinite(plot_data$neg_log10_q)] <- inf_replacement

  # Classify significance relative to the tested null
  effect <- plot_data$coef - plot_data$null_hypothesis
  plot_data$significant_flag <- !is.na(plot_data$qval) & plot_data$qval < q_threshold &
    abs(effect) > coef_threshold
  plot_data$significance_level <- dplyr::case_when(
    plot_data$significant_flag & effect > 0 ~ "Significant Up",
    plot_data$significant_flag & effect < 0 ~ "Significant Down",
    TRUE ~ "Not Significant"
  )
  plot_data$significance_level <- factor(
    plot_data$significance_level,
    levels = c("Significant Up", "Significant Down", "Not Significant")
  )

  # Determine labels: significant taxa first (by q-value), then remaining by |coef|
  plot_data$label <- ""
  if (label_significant && max_labels > 0) {
    sig_taxa <- plot_data[plot_data$significant_flag, ]
    sig_taxa <- sig_taxa[order(sig_taxa$qval), ]
    label_taxa <- utils::head(sig_taxa$taxon, max_labels)
    remaining_slots <- max_labels - length(label_taxa)
    if (remaining_slots > 0) {
      nonsig_taxa <- plot_data[!plot_data$significant_flag & !(plot_data$taxon %in% label_taxa), ]
      nonsig_taxa <- nonsig_taxa[order(-abs(nonsig_taxa$coef - nonsig_taxa$null_hypothesis)), ]
      label_taxa <- c(label_taxa, utils::head(nonsig_taxa$taxon, remaining_slots))
    }
    plot_data$label[plot_data$taxon %in% label_taxa] <- plot_data$taxon[plot_data$taxon %in% label_taxa]
  }

  # Default title
  if (is.null(title)) {
    model_label <- if (model_type == "joint") "Joint" else tools::toTitleCase(model_type)
    title <- paste0("MaAsLin3 Volcano Plot (", model_label, ")")
  }

  # Build plot
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = coef, y = neg_log10_q, color = significance_level, label = label)
  ) +
    ggplot2::geom_point(alpha = 0.7, size = 2.5) +
    ggplot2::geom_hline(
      yintercept = -log10(q_threshold),
      linetype = "dashed", color = "grey60", linewidth = 0.5
    ) +
    ggplot2::scale_color_manual(values = color_palette, name = "Significance") +
    ggplot2::labs(
      x = "Beta Coefficient",
      y = bquote(-Log[10] ~ "(Adjusted P-value)"),
      title = title,
      subtitle = paste("Coefficient:", coefficient_name)
    )

  # Null and coefficient-threshold lines (only where they differ from the default 0 / none)
  null_x <- unique(round(plot_data$null_hypothesis, 8))
  v_lines <- c(if (any(null_x != 0)) null_x,
               if (coef_threshold > 0) c(null_x - coef_threshold, null_x + coef_threshold))
  if (length(v_lines) > 0) {
    p <- p + ggplot2::geom_vline(xintercept = v_lines, linetype = "dashed", color = "grey60",
                                 linewidth = 0.5)
  }

  # Add SE error bars if requested
  if (show_se) {
    p <- p + ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = coef - stderr, xmax = coef + stderr),
      height = 0.02, alpha = 0.5, show.legend = FALSE
    )
  }

  # Add labels
  if (any(plot_data$label != "")) {
    p <- p + ggrepel::geom_text_repel(
      max.overlaps = max_overlaps,
      box.padding = 0.6,
      point.padding = 0.4,
      segment.color = "grey70",
      segment.size = 0.3,
      size = 3.5,
      min.segment.length = 0,
      fontface = if (italicize_taxa) "italic" else "plain",
      show.legend = FALSE
    )
  }

  p
}
