#' Format P-values as Significance Stars
#'
#' Single definition of the star thresholds used across the alpha-diversity, beta-diversity
#' and differential-abundance tables, so a threshold change cannot land in one table and not
#' another.
#'
#' @param p Numeric vector of p-values. NA yields "ns".
#'
#' @return Character vector of "***", "**", "*" or "ns".
#'
#' @export
classify_significance <- function(p) {
  dplyr::case_when(
    !is.na(p) & p < 0.001 ~ "***",
    !is.na(p) & p < 0.01 ~ "**",
    !is.na(p) & p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}
