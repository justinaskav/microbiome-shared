#' Rarefaction curves: expected richness at increasing depth, per sample
#'
#' Uses `vegan::rarefy()`, i.e. the exact expected number of features at each depth
#' (no random subsampling), so the result is deterministic.
#'
#' @param tse A (Tree)SummarizedExperiment.
#' @param assay_type Count assay to rarefy (default "counts").
#' @param meta_cols colData columns to carry into the result (e.g. "Diagnosis").
#' @param n_steps Number of evenly spaced depths between 100 and the deepest sample.
#' @return A tibble with columns sample, depth, richness and `meta_cols`; the depth of the
#'   shallowest sample is attached as attribute `"min_depth"`.
#' @export
calc_rarefaction_curves <- function(tse, assay_type = "counts", meta_cols = character(),
                                    n_steps = 50) {
  count_mat <- t(SummarizedExperiment::assay(tse, assay_type))  # samples as rows (vegan)
  sample_meta <- as.data.frame(SummarizedExperiment::colData(tse))

  max_depth <- max(rowSums(count_mat))
  steps <- unique(c(1, round(seq(100, max_depth, length.out = n_steps)), max_depth))

  out <- purrr::map_dfr(seq_len(nrow(count_mat)), function(i) {
    sample_total <- sum(count_mat[i, ])
    valid_steps <- steps[steps <= sample_total]
    richness <- vegan::rarefy(count_mat[i, , drop = FALSE], valid_steps)
    row <- tibble::tibble(
      sample = rownames(count_mat)[i],
      depth = valid_steps,
      richness = as.numeric(richness)
    )
    for (mc in meta_cols) row[[mc]] <- sample_meta[[mc]][i]
    row
  })
  attr(out, "min_depth") <- min(rowSums(count_mat))
  out
}

#' Plot rarefaction curves, one line per sample
#'
#' @param data Output of `calc_rarefaction_curves()`.
#' @param color_by Column to colour lines by.
#' @param colors Named colour vector for `color_by` (NULL = ggplot default).
#' @param depth Depth to mark with a dashed line (default: the shallowest sample,
#'   `attr(data, "min_depth")`); NULL for no line.
#' @param line_alpha,line_width Line geom params.
#' @return A ggplot object.
#' @export
plot_rarefaction_curves <- function(data, color_by, colors = NULL,
                                    depth = attr(data, "min_depth"),
                                    line_alpha = 0.3, line_width = 0.65) {
  p <- ggplot2::ggplot(data, ggplot2::aes(x = depth, y = richness, group = sample,
                                          color = .data[[color_by]])) +
    ggplot2::geom_line(alpha = line_alpha, linewidth = line_width)
  if (!is.null(depth)) {
    p <- p +
      ggplot2::geom_vline(xintercept = depth, linetype = "dashed", color = "grey40") +
      ggplot2::annotate("text", x = depth, y = max(data$richness) * 0.95,
                        label = paste0("Rarefaction depth: ", scales::comma(depth)),
                        hjust = -0.05, size = 3, color = "grey30")
  }
  if (!is.null(colors)) p <- p + ggplot2::scale_color_manual(values = colors)
  p + ggplot2::labs(x = "Sequencing Depth", y = "Observed ASVs", color = color_by)
}
