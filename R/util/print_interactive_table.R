#' Print Interactive GT Table
#'
#' Renders a gt table as interactive HTML (with compact mode and page size selector)
#' when the output format is HTML, otherwise uses a smaller font size for static output.
#'
#' @param tbl A gt table object.
#'
#' @return A gt table object with interactive or static formatting applied.
print_interactive_table <- function(tbl) {
  if (knitr::is_html_output()) {
    gt::opt_interactive(tbl, use_compact_mode = TRUE, use_page_size_select = TRUE)
  } else {
    gt::tab_options(tbl, table.font.size = gt::px(10))
  }
}
