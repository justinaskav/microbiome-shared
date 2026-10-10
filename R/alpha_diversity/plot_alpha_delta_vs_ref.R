#' Per-subject change in an alpha index relative to the subject's own reference site
#'
#' A site LMM with a subject random intercept tests a WITHIN-subject contrast, which raw
#' values with one line per subject do not show (they show mostly between-subject spread).
#' This plots each subject's own difference against their own reference site, with a zero
#' line for "no difference", and overlays the model-adjusted contrast (point + 95% CI).
#'
#' `posthoc` is the `$posthoc` frame from `calc_alpha_diversity_summary()` run with
#' `ref_group = ref_site`. emmeans writes each comparison in factor-level order, so it reads
#' "<ref> - X" when the reference is the earlier level and "X - <ref>" otherwise. Estimates of
#' "<ref> - X" rows are negated (CI bounds swapped) to express "X relative to ref". Getting this
#' backwards would silently reverse the direction of the result, so exactly one matching
#' comparison per site is asserted rather than assumed.
#'
#' @param data Data frame with `subject_var`, `site_var` and the `index` column.
#' @param index Alpha-index column.
#' @param posthoc Post-hoc frame (see above).
#' @param ref_site Reference site level.
#' @param sites Sites to plot against the reference, in order.
#' @param colors Optional named colour vector for `sites`.
#' @param style Cosmetic spec from `alpha_plot_style()`.
#' @param subject_var,site_var Column names.
#'
#' Requires `split_contrast()` (calc_alpha_diversity_summary.R).
#' @return A ggplot object.
#' @export
plot_alpha_delta_vs_ref <- function(data, index, posthoc, ref_site, sites,
                                    colors = NULL,
                                    style = alpha_plot_style(),
                                    subject_var = "Patient", site_var = "Site") {
  wide <- data |>
    dplyr::select(Patient = dplyr::all_of(subject_var), Site = dplyr::all_of(site_var),
                  dplyr::all_of(index)) |>
    tidyr::pivot_wider(names_from = Site, values_from = dplyr::all_of(index),
                       values_fn = function(x) mean(x, na.rm = TRUE))
  stopifnot(ref_site %in% names(wide))

  deltas <- purrr::map_dfr(sites, function(s) {
    if (!s %in% names(wide)) return(NULL)
    tibble::tibble(Site = s, Patient = wide$Patient, delta = wide[[s]] - wide[[ref_site]])
  }) |>
    dplyr::filter(!is.na(delta)) |>
    dplyr::mutate(Site = factor(Site, levels = sites))

  # Model contrasts as "level - ref" (sign-flipped where emmeans wrote "ref - level").
  contr <- purrr::map_dfr(sites, function(s) {
    g <- split_contrast(posthoc$comparison)
    flip <- g$group1 == ref_site & g$group2 == s
    keep <- flip | (g$group1 == s & g$group2 == ref_site)
    if (sum(keep) != 1) {
      stop("expected exactly one posthoc row comparing '", ref_site, "' and '", s, "'; found ",
           sum(keep), ". Was calc_alpha_diversity_summary() run with ref_group = '", ref_site, "'?")
    }
    row <- posthoc[keep, , drop = FALSE]
    sgn <- if (flip[keep]) -1 else 1
    tibble::tibble(
      Site  = s,
      est   = sgn * row$estimate,
      lower = if (sgn < 0) -row$upper.CL else row$lower.CL,   # bounds swap when flipped
      upper = if (sgn < 0) -row$lower.CL else row$upper.CL,
      sig   = row$significance
    )
  }) |>
    dplyr::mutate(Site = factor(Site, levels = sites))

  # Same headroom rule as plot_alpha_diversity_summary_lmm()'s significance brackets: a step
  # above the max, sized as a fraction of the y-range rather than a flat multiplier, so it does
  # not collapse when the max sits near zero.
  max_y   <- max(deltas$delta, na.rm = TRUE)
  y_range <- max_y - min(deltas$delta, na.rm = TRUE)
  star_y  <- max_y + (if (y_range > 0) y_range * style$sig_step else 0.1)
  # The box spans +/- box_width/2 around the category centre, so the contrast diamond sits
  # just outside it and well short of the next category's jitter.
  nudge <- style$box_width / 2 + 0.10

  p <- ggplot2::ggplot(deltas, ggplot2::aes(x = Site, y = delta)) +
    do.call(ggplot2::geom_hline, c(list(yintercept = 0), style$zero_line)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = Site), width = style$box_width,
                          outlier.shape = NA, alpha = style$box_alpha,
                          colour = style$box_colour, linewidth = style$box_linewidth) +
    ggplot2::geom_jitter(ggplot2::aes(colour = Site), width = style$jitter_width, height = 0,
                         size = style$point_size, alpha = style$point_alpha) +
    ggplot2::geom_errorbar(data = contr,
                           ggplot2::aes(x = as.numeric(Site) + nudge, y = est,
                                        ymin = lower, ymax = upper),
                           width = style$errorbar_width, linewidth = style$errorbar_linewidth,
                           colour = style$contrast_colour, inherit.aes = FALSE) +
    ggplot2::geom_point(data = contr,
                        ggplot2::aes(x = as.numeric(Site) + nudge, y = est),
                        size = style$contrast_size, shape = style$contrast_shape,
                        colour = style$contrast_colour, inherit.aes = FALSE) +
    ggplot2::geom_text(data = contr[contr$sig != "ns", , drop = FALSE],
                       ggplot2::aes(x = as.numeric(Site) + nudge, y = star_y, label = sig),
                       size = style$sig_size, inherit.aes = FALSE) +
    # Right-hand axis spells out what the sign means, so "which way is up" does not depend on
    # the reader working out the subtraction order from the y label.
    ggplot2::scale_y_continuous(
      expand = style$y_expand,
      sec.axis = ggplot2::dup_axis(
        name = NULL,
        breaks = c(max_y * 0.6, 0, min(deltas$delta, na.rm = TRUE) * 0.6),
        labels = c(paste0("higher\nthan ", tolower(ref_site)),
                   paste0("= ", tolower(ref_site)),
                   paste0("lower\nthan ", tolower(ref_site))))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = paste0(site_var, ", relative to the same ", tolower(subject_var), "'s ",
                             tolower(ref_site)),
                  y = paste0("Δ ", alpha_index_label(index), " (site − ", tolower(ref_site), ")")) +
    ggplot2::theme(legend.position = "none",
                   axis.ticks.y.right = ggplot2::element_blank(),
                   axis.text.y.right = ggplot2::element_text(colour = "grey35", size = ggplot2::rel(0.8),
                                                             lineheight = 0.9, hjust = 0))

  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors, na.value = "grey70") +
      ggplot2::scale_colour_manual(values = colors, na.value = "grey70")
  }
  p
}
