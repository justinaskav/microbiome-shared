#' Cosmetic spec for alpha-diversity panels
#'
#' Bundles the geom/theme constants that `plot_alpha_diversity_summary_lmm()` and
#' `create_alpha_diversity_panel()` apply to every panel (box width/colour, point
#' size, significance-label sizing, y-axis headroom, ...), so a caller can enforce
#' one look across several figures by building the list once and passing it
#' everywhere instead of hardcoding numbers inside the plotting functions.
#'
#' All defaults reproduce the values that were previously hardcoded in
#' `plot_alpha_diversity_summary_lmm()`, so existing callers that do not pass
#' `style` render identically to before this was split out.
#'
#' @param box_width,box_colour,box_alpha,box_linewidth Boxplot geom params.
#' @param point_size,point_alpha,jitter_width Point/jitter geom params (used for
#'   both `geom_jitter` and the non-jittered `geom_point` when `connect_subjects
#'   = TRUE`).
#' @param sig_size Text size of significance-bracket labels (e.g. "*", "**").
#' @param sig_step,sig_tip Bracket vertical spacing and tip length, each as a
#'   fraction of the y-range.
#' @param label_size Text size of the top-right LMM annotation and bottom-right
#'   post-hoc status label.
#' @param y_expand `expansion()` value for `scale_y_continuous()`, giving the
#'   brackets/labels headroom above the data.
#' @param connector_colour,connector_alpha,connector_linewidth `geom_line()`
#'   params for the `connect_subjects = TRUE` per-patient connectors.
#' @param zero_line Named list of `geom_hline()` params (used by delta-vs-reference plots,
#'   e.g. `plot_alpha_delta_vs_ref()`, to mark "no difference").
#' @param contrast_colour,contrast_size,contrast_shape Point params for the overlaid
#'   model-contrast estimate on delta-vs-reference plots.
#' @param errorbar_width,errorbar_linewidth `geom_errorbar()` params for the model-contrast
#'   95% CI on delta-vs-reference plots.
#'
#' @return A named list of style parameters.
#' @export
alpha_plot_style <- function(
  box_width = 0.75, box_colour = "black", box_alpha = 0.35, box_linewidth = 0.5,
  point_size = 2.5, point_alpha = 0.6, jitter_width = 0.2,
  sig_size = 3.5, sig_step = 0.08, sig_tip = 0.015,
  label_size = 2.8,
  y_expand = ggplot2::expansion(mult = c(0.05, 0.20)),
  connector_colour = "grey50", connector_alpha = 0.4, connector_linewidth = 0.5,
  # delta-vs-reference plots only (e.g. plot_alpha_delta_vs_ref())
  zero_line = list(linetype = "dashed", colour = "grey45", linewidth = 0.45),
  contrast_colour = "black", contrast_size = 2.3, contrast_shape = 18,
  errorbar_width = 0.09, errorbar_linewidth = 0.65
) {
  list(
    box_width = box_width, box_colour = box_colour, box_alpha = box_alpha,
    box_linewidth = box_linewidth,
    point_size = point_size, point_alpha = point_alpha, jitter_width = jitter_width,
    sig_size = sig_size, sig_step = sig_step, sig_tip = sig_tip,
    label_size = label_size,
    y_expand = y_expand,
    connector_colour = connector_colour, connector_alpha = connector_alpha,
    connector_linewidth = connector_linewidth,
    zero_line = zero_line,
    contrast_colour = contrast_colour, contrast_size = contrast_size,
    contrast_shape = contrast_shape,
    errorbar_width = errorbar_width, errorbar_linewidth = errorbar_linewidth
  )
}
