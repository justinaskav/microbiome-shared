#' Subset a TSE's samples together with its stored beta-diversity results
#'
#' Distance matrices in `metadata()` do not follow column subsetting, and PCoA coordinates
#' and eigenvalues of the full object are wrong for a subset. This subsets each stored
#' distance to the kept samples and recomputes its `MDS_*` reducedDim (k = 2) and, when
#' present, its `<dist>_eig` eigenvalues (used for the axis variance-explained labels).
#'
#' @param tse TSE carrying distances in `metadata()` (as written by the workflows' `calc_beta_diversity.R`).
#' @param keep Logical or index vector of samples (columns) to keep.
#' @param distances Metadata keys of the distances to subset.
#' @return The subset TSE.
#' @export
subset_beta_metadata <- function(tse, keep,
                                 distances = c("jaccard", "bray", "unifrac", "unifrac_weighted")) {
  tse <- tse[, keep]
  keep_names <- colnames(tse)
  for (dn in intersect(distances, names(tse@metadata))) {
    d <- as.dist(as.matrix(tse@metadata[[dn]])[keep_names, keep_names])
    tse@metadata[[dn]] <- d
    mds_nm <- paste0("MDS_", dn)
    if (mds_nm %in% SingleCellExperiment::reducedDimNames(tse)) {
      mds <- stats::cmdscale(d, k = 2, eig = TRUE)
      SingleCellExperiment::reducedDim(tse, mds_nm) <- mds$points
      eig_key <- paste0(dn, "_eig")
      if (eig_key %in% names(tse@metadata)) tse@metadata[[eig_key]] <- mds$eig
    }
  }
  tse
}
