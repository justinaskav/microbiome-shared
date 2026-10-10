#' Calculate NetCoMi Network Comparison Between Two Groups
#'
#' Constructs two microbial co-occurrence networks (one per group) and compares them
#' using NetCoMi's network comparison framework.
#'
#' @param tse A TreeSummarizedExperiment object.
#' @param group_var Character string specifying the column name in colData for grouping.
#' @param group1 Character string specifying the first group (reference).
#' @param group2 Character string specifying the second group (comparison).
#' @param tax_rank Character string specifying the taxonomic rank (default: "genus").
#' @param filter_pattern Character pattern to filter out taxa (default: "Incertae_Sedis|Unknown|uncultured").
#' @param prevalence_threshold Numeric, minimum sample proportion for taxa filtering (default: 0.2).
#' @param measure Character, association measure for NetCoMi (default: "spring").
#' @param measure_par List of parameters for the measure.
#' @param spars_method Character, sparsification method (default: "none").
#' @param diss_func Character, dissimilarity function (default: "signed").
#' @param min_samples Integer, minimum samples per group (default: 10).
#' @param compute_diffnet Logical, whether to compute differential network via
#'   NetCoMi::diffnet() (default: FALSE).
#' @param diff_method Character, method for differential association testing
#'   (default: "fisherTest"). Passed to diffnet().
#' @param diff_adjust Character, p-value adjustment method (default: "lfdr").
#'   Falls back to "BH" if lfdr fails.
#' @param seed Numeric, random seed for reproducibility (default: 13075).
#'
#' @return A list containing:
#'   - netcomi_net: NetCoMi network object with two networks
#'   - netcomi_props: NetCoMi network properties for both networks
#'   - netcomi_comp: NetCoMi network comparison object
#'   - netcomi_diffnet: NetCoMi diffnet object (NULL if compute_diffnet = FALSE)
#'   - group1: Name of first group
#'   - group2: Name of second group
#'
#' @export
calc_netcomi_comparison <- function(tse,
                                    group_var,
                                    group1,
                                    group2,
                                    tax_rank = "genus",
                                    filter_pattern = "Incertae_Sedis|Unknown|uncultured",
                                    prevalence_threshold = 0.2,
                                    measure = "spring",
                                    measure_par = list(thresh = 0.05, Rmethod = "approx"),
                                    spars_method = "none",
                                    diss_func = "signed",
                                    min_samples = 10,
                                    compute_diffnet = FALSE,
                                    diff_method = "fisherTest",
                                    diff_adjust = "lfdr",
                                    seed = 13075) {

  # Check if NetCoMi is available
  if (!requireNamespace("NetCoMi", quietly = TRUE)) {
    stop("NetCoMi package is required but not installed.")
  }

  # --- 1. Data Preparation ---
  # Filter out unwanted taxa
  if (!is.null(filter_pattern) && nchar(filter_pattern) > 0) {
    rd <- SummarizedExperiment::rowData(tse)
    if (tax_rank %in% colnames(rd)) {
      selected_rows <- grepl(filter_pattern, rd[[tax_rank]], ignore.case = TRUE)
      tse <- tse[!selected_rows, ]
    }
  }

  # Subset to only the two groups of interest
  cd <- SummarizedExperiment::colData(tse)
  if (!group_var %in% colnames(cd)) {
    stop("group_var '", group_var, "' not found in colData")
  }

  group_vec <- as.character(cd[[group_var]])
  keep_samples <- group_vec %in% c(group1, group2)

  if (sum(keep_samples) < 6) {
    stop("Insufficient samples after filtering. Need at least 6 samples total, got ", sum(keep_samples))
  }

  tse_subset <- tse[, keep_samples]

  # Check sample counts per group
  cd_subset <- SummarizedExperiment::colData(tse_subset)
  group_counts <- table(cd_subset[[group_var]])

  if (!group1 %in% names(group_counts) || group_counts[group1] < min_samples) {
    stop("Insufficient samples in group1 '", group1, "'. Need at least ", min_samples, " samples, got ",
         ifelse(group1 %in% names(group_counts), group_counts[group1], 0), ".")
  }
  if (!group2 %in% names(group_counts) || group_counts[group2] < min_samples) {
    stop("Insufficient samples in group2 '", group2, "'. Need at least ", min_samples, " samples, got ",
         ifelse(group2 %in% names(group_counts), group_counts[group2], 0), ".")
  }

  # Warn about SPRING reliability with small samples
  for (g in c(group1, group2)) {
    if (group_counts[g] < 20) {
      warning("Group '", g, "' has only ", group_counts[g],
              " samples. SPRING partial correlation estimation is unreliable with <20-30 samples.",
              call. = FALSE)
    }
  }

  message("Group counts: ", group1, " = ", group_counts[group1], ", ", group2, " = ", group_counts[group2])

  # Agglomerate to specified rank
  tse_agg <- mia::agglomerateByRank(tse_subset, rank = tax_rank)

  # Store rowData for later use in plotting
  rd_original <- SummarizedExperiment::rowData(tse_agg)

  # --- 2. Prepare data for NetCoMi ---
  # NetCoMi expects: samples as rows, taxa as columns
  count_mat <- t(SummarizedExperiment::assay(tse_agg, "counts"))

  # Extract the group vector for samples

  cd_agg <- SummarizedExperiment::colData(tse_agg)
  group_vec_final <- as.character(cd_agg[[group_var]])

  # Get taxa names for column labels
  taxa_names <- rownames(SummarizedExperiment::assay(tse_agg, "counts"))
  colnames(count_mat) <- taxa_names

  message("Network construction: ", nrow(count_mat), " samples x ", ncol(count_mat), " taxa")
  message("Using method: ", measure, " (this may take several minutes...)")

  # --- 3. Network Construction (Two Networks) ---
  netcomi_net <- NetCoMi::netConstruct(
    data = count_mat,
    dataType = "counts",
    group = group_vec_final,
    filtTax = "numbSamp",
    filtTaxPar = list(numbSamp = prevalence_threshold),
    filtSamp = "none",
    measure = measure,
    measurePar = measure_par,
    sparsMethod = spars_method,
    dissFunc = diss_func,
    seed = seed,
    verbose = 2
  )

  # --- 3. Network Analysis ---
  netcomi_props <- NetCoMi::netAnalyze(
    netcomi_net,
    clustMethod = "cluster_fast_greedy",
    hubPar = "eigenvector",
    normDeg = FALSE
  )

  # --- 4. Network Comparison ---
  netcomi_comp <- NetCoMi::netCompare(
    netcomi_props,
    permTest = FALSE,  # Skip permutation test for speed
    verbose = FALSE
  )

  # --- 5. Differential Network (optional) ---
  netcomi_diffnet <- NULL
  if (compute_diffnet) {
    message("Computing differential network (", diff_method, ", adjust = ", diff_adjust, ")...")
    netcomi_diffnet <- tryCatch(
      NetCoMi::diffnet(netcomi_net, diffMethod = diff_method, adjust = diff_adjust),
      error = function(e) {
        warning("diffnet with '", diff_adjust, "' failed: ", e$message,
                ". Falling back to 'BH'.", call. = FALSE)
        NetCoMi::diffnet(netcomi_net, diffMethod = diff_method, adjust = "BH")
      }
    )
    n_diff <- sum(netcomi_diffnet$diffAdjustMat[upper.tri(netcomi_diffnet$diffAdjustMat)] != 0)
    message("Differential associations found: ", n_diff)
  }

  return(list(
    netcomi_net = netcomi_net,
    netcomi_props = netcomi_props,
    netcomi_comp = netcomi_comp,
    netcomi_diffnet = netcomi_diffnet,
    group1 = group1,
    group2 = group2,
    tse_agg = tse_agg,
    tax_rank = tax_rank,
    rowData = rd_original
  ))
}


#' Plot NetCoMi Network Comparison
#'
#' Plots a side-by-side comparison of two microbial networks from calc_netcomi_comparison.
#'
#' @param netcomi_result A list returned by calc_netcomi_comparison.
#' @param cex_nodes Numeric, node size scaling factor (default: 1.5).
#' @param node_transp Numeric, node transparency 0-100 (default: 20).
#' @param edge_transp_high Numeric, positive edge transparency (default: 20).
#' @param edge_transp_low Numeric, negative edge transparency (default: 20).
#' @param repulsion Numeric, layout repulsion factor (default: 0.95).
#' @param show_legends Logical, whether to show legends (default: TRUE).
#'
#' @return Invisibly returns NULL. Called for side effect of plotting.
#'
#' @export
plot_netcomi_comparison <- function(netcomi_result,
                                    cex_nodes = 1.5,
                                    node_transp = 20,
                                    edge_transp_high = 20,
                                    edge_transp_low = 20,
                                    repulsion = 0.95,
                                    show_legends = TRUE) {

  if (!requireNamespace("NetCoMi", quietly = TRUE)) {
    stop("NetCoMi package is required but not installed.")
  }

  netcomi_props <- netcomi_result$netcomi_props
  rd <- netcomi_result$rowData
  tax_rank <- netcomi_result$tax_rank
  group1 <- netcomi_result$group1
  group2 <- netcomi_result$group2

  # --- Prepare Node Coloring by Phylum ---
  # Get taxa names that are in the network
  network_taxa <- netcomi_props$input$assoMat1 |> rownames()

  # Match to rowData to get phylum info
  if ("phylum" %in% colnames(rd) && tax_rank %in% colnames(rd)) {
    # Build a lookup table from rowData
    taxa_to_phylum <- setNames(as.character(rd$phylum), rownames(rd))

    # Extract phyla for network taxa
    phyla_vec <- taxa_to_phylum[network_taxa]
    phyla_vec[is.na(phyla_vec)] <- "Unknown"
    phyla <- as.factor(phyla_vec)
    names(phyla) <- network_taxa

    # Create color palette for phyla
    unique_phyla <- levels(phyla)
    n_phyla <- length(unique_phyla)
    colvec <- colorRampPalette(RColorBrewer::brewer.pal(min(n_phyla, 9), "Set1"))(n_phyla)
    names(colvec) <- unique_phyla

    use_phylum_colors <- TRUE
  } else {
    use_phylum_colors <- FALSE
    n_phyla <- 0
  }

  # --- Plot Comparison ---
  if (use_phylum_colors) {
    plot(netcomi_props,
      sameLayout = TRUE,
      layoutGroup = "union",
      repulsion = repulsion,
      rmSingles = "inboth",
      labelScale = FALSE,
      nodeSize = "mclr",
      nodeColor = "feature",
      featVecCol = phyla,
      colorVec = colvec,
      nodeTransp = node_transp,
      highlightHubs = FALSE,
      cexNodes = cex_nodes,
      edgeTranspHigh = edge_transp_high,
      edgeTranspLow = edge_transp_low,
      title1 = group1,
      title2 = group2,
      showTitle = TRUE,
      cexTitle = 1.3,
      mar = c(1, 3, 3, 8)
    )
  } else {
    plot(netcomi_props,
      sameLayout = TRUE,
      layoutGroup = "union",
      repulsion = repulsion,
      rmSingles = "inboth",
      labelScale = FALSE,
      nodeSize = "mclr",
      nodeTransp = node_transp,
      highlightHubs = FALSE,
      cexNodes = cex_nodes,
      edgeTranspHigh = edge_transp_high,
      edgeTranspLow = edge_transp_low,
      title1 = group1,
      title2 = group2,
      showTitle = TRUE,
      cexTitle = 1.3,
      mar = c(1, 3, 3, 8)
    )
  }

  # --- Add Legends ---
  if (show_legends) {
    # Correlation sign legend
    legend(
      "bottomright",
      inset = c(-0.02, 0),
      cex = 0.7,
      title = "Correlation:",
      legend = c("Positive", "Negative"),
      lty = 1,
      lwd = 2,
      col = c("#009900", "red"),
      bty = "n",
      horiz = FALSE
    )

    # Phylum legend (only if using phylum colors and not too many)
    if (use_phylum_colors && n_phyla <= 12) {
      col_transp <- NetCoMi::colToTransp(colvec, node_transp)
      legend(
        "topright",
        inset = c(-0.02, 0),
        cex = 0.6,
        pt.cex = 1.2,
        title = "Phylum:",
        legend = names(colvec),
        fill = col_transp,
        bty = "n"
      )
    }
  }

  return(invisible(NULL))
}
