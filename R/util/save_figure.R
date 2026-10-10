#' Save a figure as <dir>/<name>.png on a white background
#'
#' Explicit-plot replacement for bare `ggsave()` calls: the plot is always passed in
#' (never "the last plot drawn") and the background is white even for themes whose
#' plot background is transparent.
#'
#' @param plot ggplot or patchwork object.
#' @param name File name without extension.
#' @param width,height Size in inches.
#' @param dpi Resolution (default 450).
#' @param dir Output directory, created if missing (default "figures").
#' @param bg Background colour (default "white").
#' @return The plot, invisibly, so the call can end a pipeline.
#' @export
save_figure <- function(plot, name, width, height, dpi = 450, dir = "figures", bg = "white") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file.path(dir, paste0(name, ".png")), plot,
                  width = width, height = height, dpi = dpi, bg = bg)
  invisible(plot)
}

#' Evaluate a plotting expression, reporting errors instead of raising them
#'
#' @param expr Expression to evaluate (lazily).
#' @return The value of `expr`, or `NULL` (after printing the message) on error.
#' @export
safe_render <- function(expr) {
  tryCatch(expr, error = function(e) cat("Error:", e$message, "\n"))
}
