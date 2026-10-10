#' Estimated marginal means for a two-factor interaction fit, without connecting lines
#'
#' Points + 95% CI per `group_var` x `x_var` cell, dodged. No lines join the x categories:
#' when they are not a continuum (e.g. mucosal sites followed by faeces) a line would draw a
#' trajectory that does not exist. A dotted rule can mark where one series ends
#' (`separate_before`).
#'
#' The secondary axis back-transforms the response. MaAsLin3's default "LOG" transform is
#' log2 of TSS relative abundance, so the default is `2^y * 100` (%). It needs EXPLICIT breaks:
#' left to itself ggplot picks breaks evenly spaced in percent, which the exponential
#' back-transform crams into the top of the axis. Decade breaks are evenly spaced in log2;
#' out-of-range breaks are dropped automatically. Note the back-transformed value is a
#' geometric mean conditional on detection (the abundance model fits non-zero observations
#' only), not the arithmetic mean relative abundance.
#'
#' @param model A fitted model emmeans can use (e.g. `maaslin3 result$fit_data_abundance$fits[[taxon]]`).
#' @param title Plot title (italic, for taxon names).
#' @param group_var,x_var Interaction factors: colour/shape groups and x axis.
#' @param group_levels,x_levels Factor levels (cells outside `x_levels` are dropped).
#' @param colors,shapes Named colour / shape vectors for `group_levels`.
#' @param separate_before Optional `x_var` level; a dotted rule is drawn before it.
#' @param sec_trans,sec_breaks,sec_labels,sec_name Secondary axis; `sec_trans = NULL` for none.
#' @param x_label,y_label Axis labels.
#' @return A ggplot object.
#' @export
plot_emm_interaction <- function(model, title,
                                 group_var = "Diagnosis", x_var = "Site",
                                 group_levels, x_levels,
                                 colors, shapes,
                                 separate_before = NULL,
                                 sec_trans = ~ 2^(.) * 100,
                                 sec_breaks = c(0.001, 0.01, 0.1, 1, 10),
                                 sec_labels = c("0.001", "0.01", "0.1", "1", "10"),
                                 sec_name = "Relative abundance (%)",
                                 x_label = x_var,
                                 y_label = expression(EMM ~ (log[2] ~ relative ~ abundance))) {
  emm <- emmeans::emmeans(model, stats::as.formula(paste("~", group_var, "*", x_var)))
  s <- as.data.frame(summary(emm))
  s[[x_var]] <- factor(s[[x_var]], levels = x_levels)
  s[[group_var]] <- factor(s[[group_var]], levels = group_levels)
  s <- s |> dplyr::filter(!is.na(.data[[x_var]]), !is.na(emmean))
  pd <- ggplot2::position_dodge(width = 0.45)

  p <- ggplot2::ggplot(s, ggplot2::aes(x = .data[[x_var]], y = emmean,
                                       color = .data[[group_var]], shape = .data[[group_var]],
                                       group = .data[[group_var]]))

  present <- levels(droplevels(s[[x_var]]))
  if (!is.null(separate_before) && separate_before %in% present && length(present) > 1) {
    p <- p + ggplot2::geom_vline(xintercept = which(present == separate_before) - 0.5,
                                 linetype = "dotted", colour = "grey55", linewidth = 0.4)
  }

  y_scale <- if (is.null(sec_trans)) {
    ggplot2::scale_y_continuous()
  } else {
    ggplot2::scale_y_continuous(
      sec.axis = ggplot2::sec_axis(sec_trans, name = sec_name,
                                   breaks = sec_breaks, labels = sec_labels)
    )
  }

  p +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower.CL, ymax = upper.CL), width = 0.15,
                           linewidth = 0.6, position = pd) +
    ggplot2::geom_point(size = 3, position = pd) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_shape_manual(values = shapes) +
    y_scale +
    ggplot2::labs(x = x_label, y = y_label, title = title,
                  color = group_var, shape = group_var) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "italic"))
}
