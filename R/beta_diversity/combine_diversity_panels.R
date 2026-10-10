#' Lay out a beta-diversity panel's plots as one tagged patchwork figure
#'
#' @param panel Output of `create_beta_diversity_panel()` (uses `$plots`).
#' @param ncol,nrow Layout.
#' @param legend_position Position of the collected legend.
#' @return A patchwork object.
#' @export
wrap_beta_panel <- function(panel, ncol = 2, nrow = 1, legend_position = "right") {
  patchwork::wrap_plots(panel$plots, ncol = ncol, nrow = nrow, guides = "collect") &
    patchwork::plot_annotation(tag_levels = "a") &
    ggplot2::theme(legend.position = legend_position)
}

#' Stack a 2-plot alpha panel over a 2-plot beta panel as one 2x2 figure with one legend
#'
#' `guides = "collect"` merges guides only when their key glyphs match as well as their
#' titles and labels, so a legend drawn from boxplots (rect keys) never merges with the same
#' legend drawn from PCoA points (point keys); pinning the scales is not enough. The beta
#' panels therefore donate the legend (point keys read correctly for both rows) and the alpha
#' panels drop theirs.
#'
#' @param alpha_panel Patchwork from `create_alpha_diversity_panel()` (first two plots used).
#' @param beta_panel Output of `create_beta_diversity_panel()` (first two plots used).
#' @param color_values,shape_values Named colour / shape vectors shared by all four plots.
#' @param color_title,shape_title Legend titles.
#' @param legend_position Position of the collected legend.
#' @return A patchwork object.
#' @export
combine_alpha_beta_panels <- function(alpha_panel, beta_panel, color_values,
                                      shape_values, color_title, shape_title,
                                      legend_position = "right") {
  shared_scales <- list(
    ggplot2::scale_color_manual(values = color_values, limits = names(color_values),
                                drop = FALSE, na.value = "grey50", name = color_title),
    ggplot2::scale_fill_manual(values = color_values, limits = names(color_values),
                               drop = FALSE, na.value = "grey50", name = color_title),
    ggplot2::scale_shape_manual(values = shape_values, limits = names(shape_values),
                                drop = FALSE, name = shape_title)
  )
  alpha_plots <- lapply(list(alpha_panel[[1]], alpha_panel[[2]]), function(p) {
    p + shared_scales +
      ggplot2::guides(colour = "none", fill = "none", shape = "none")
  })
  beta_plots <- lapply(beta_panel$plots[1:2], function(p) p + shared_scales)

  patchwork::wrap_plots(c(alpha_plots, beta_plots), ncol = 2, nrow = 2) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "a") &
    ggplot2::theme(legend.position = legend_position)
}
