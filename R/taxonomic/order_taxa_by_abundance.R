order_taxa_by_abundance <- function(df, taxon_col = "FeatureID", abundance_col = "counts",
                                    remainder_string = "(Remainder)",
                                    unclassified_string = "(Unclassified)") {

  special_labels <- c(remainder_string, unclassified_string)

  # Calculate mean abundance for each taxon (excluding special labels)
  mean_abundance <- df |>
    dplyr::filter(!(.data[[taxon_col]] %in% special_labels)) |>
    dplyr::group_by(.data[[taxon_col]]) |>
    dplyr::summarise(mean_abund = mean(.data[[abundance_col]], na.rm = TRUE), .groups = "drop")

  # Check which special labels exist
  remainder_exists <- remainder_string %in% unique(df[[taxon_col]])
  unclassified_exists <- unclassified_string %in% unique(df[[taxon_col]])

  # Order levels by decreasing mean abundance - ensure character conversion
  ordered_levels <- as.character(mean_abundance[[taxon_col]][order(mean_abundance$mean_abund, decreasing = TRUE)])

  # Add remainder then unclassified at the end
  if (remainder_exists) {
    ordered_levels <- c(ordered_levels, remainder_string)
  }
  if (unclassified_exists) {
    ordered_levels <- c(ordered_levels, unclassified_string)
  }

  # Return df with ordered factor levels
  df |>
    dplyr::mutate(
      !!taxon_col := factor(as.character(.data[[taxon_col]]), levels = ordered_levels)
    )
}
