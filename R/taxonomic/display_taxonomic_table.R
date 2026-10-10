#' Display Taxonomic Summary as a GT Table
#'
#' Creates a formatted gt table from `calc_taxonomic_summary()` output,
#' with optional MaAsLin3 differential abundance statistics.
#'
#' @param summary_data Data frame from `calc_taxonomic_summary()` containing
#'   columns: taxon, group, formatted_median.
#' @param title Character string for the table title.
#' @param subtitle Character string for the table subtitle (default: "Median [IQR]").
#' @param taxon_label Character string for the taxon column label (default: "Taxon").
#' @param sort_descending Logical, whether to sort by first group's abundance descending (default: TRUE).
#' @param maaslin3_data Optional data frame with MaAsLin3 results containing columns:
#'   taxon, formatted_coef, formatted_pval, formatted_qval.
#'
#' @return A gt table object.
display_taxonomic_table <- function(
  summary_data,
  title,
  subtitle = "Median [IQR]",
  taxon_label = "Taxon",
  sort_descending = TRUE,
  maaslin3_data = NULL
) {
  wide_data <- summary_data |>
    dplyr::select(taxon, group, formatted_median) |>
    tidyr::pivot_wider(names_from = group, values_from = formatted_median)

  # Sort by first group column (extract numeric from "5.2% [2.1-8.3]")
  if (sort_descending) {
    first_group_col <- names(wide_data)[2]
    wide_data <- wide_data |>
      dplyr::mutate(sort_value = as.numeric(stringr::str_extract(.data[[first_group_col]], "^[0-9.]+"))) |>
      dplyr::arrange(dplyr::desc(sort_value)) |>
      dplyr::select(-sort_value)
  } else {
    wide_data <- wide_data |> dplyr::arrange(taxon)
  }

  # Join MaAsLin3 statistics if provided
  if (!is.null(maaslin3_data) && nrow(maaslin3_data) > 0) {
    wide_data <- wide_data |>
      dplyr::left_join(
        maaslin3_data |> dplyr::select(taxon, formatted_coef, formatted_pval, formatted_qval),
        by = "taxon"
      )
  }

  # Create gt table
  gt_table <- wide_data |>
    gt::gt() |>
    gt::tab_header(title = title, subtitle = subtitle) |>
    gt::cols_label(taxon = taxon_label)

  # Add labels for MaAsLin3 columns if present
  if (!is.null(maaslin3_data) && nrow(maaslin3_data) > 0) {
    gt_table <- gt_table |>
      gt::cols_label(
        formatted_coef = "Coefficient (\u03b2 \u00b1 SE)",
        formatted_pval = "P-value",
        formatted_qval = "Q-value (FDR)"
      )
  }

  return(gt_table)
}
