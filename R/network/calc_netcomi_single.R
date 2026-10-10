#' Calculate Single NetCoMi Network
#'
#' Constructs and analyzes a single microbial co-occurrence network using NetCoMi
#' with SPRING associations. Returns computation results without plotting.
#'
#' @param tse A TreeSummarizedExperiment object.
#' @param tax_rank Character string specifying the taxonomic rank (default: "genus").
#' @param filter_pattern Character pattern to filter out taxa (default: "Incertae_Sedis|Unknown|uncultured").
#' @param prevalence_threshold Numeric, minimum prevalence for taxa filtering (default: 0.3).
#' @param measure Character, association measure for NetCoMi (default: "spring").
#' @param measure_par List of parameters for the measure.
#' @param spars_method Character, sparsification method (default: "none").
#' @param diss_func Character, dissimilarity function (default: "signed").
#' @param clust_method Character, clustering method for netAnalyze (default: "cluster_fast_greedy").
#' @param hub_par Character, hub definition parameter (default: "eigenvector").
#' @param seed Numeric, random seed for reproducibility (default: 13075).
#'
#' @return A list containing:
#'   - netcomi_net: NetCoMi network object
#'   - netcomi_props: NetCoMi network properties
#'   - tse_agg: Aggregated TreeSummarizedExperiment used for the network
#'
#' @export
calc_netcomi_single <- function(tse,
                                tax_rank = "genus",
                                filter_pattern = "Incertae_Sedis|Unknown|uncultured",
                                prevalence_threshold = 0.3,
                                measure = "spring",
                                measure_par = list(thresh = 0.05, Rmethod = "approx"),
                                spars_method = "none",
                                diss_func = "signed",
                                clust_method = "cluster_fast_greedy",
                                hub_par = "eigenvector",
                                seed = 13075) {

  if (!requireNamespace("NetCoMi", quietly = TRUE)) {
    stop("NetCoMi package is required but not installed.")
  }

  # --- 1. Data Preparation ---
  if (!is.null(filter_pattern) && nchar(filter_pattern) > 0) {
    rd <- SummarizedExperiment::rowData(tse)
    if (tax_rank %in% colnames(rd)) {
      selected_rows <- grepl(filter_pattern, rd[[tax_rank]], ignore.case = TRUE)
      tse <- tse[!selected_rows, ]
    }
  }

  # Agglomerate to specified rank
  tse_agg <- mia::agglomerateByRank(tse, rank = tax_rank)

  # Rename taxa for cleaner display
  tse_agg <- NetCoMi::renameTaxa(tse_agg, pat = "<name>")

  message("Network construction: ", ncol(tse_agg), " samples x ", nrow(tse_agg), " taxa")
  message("Using method: ", measure, " with sparsification: ", spars_method)

  # --- 2. Network Construction ---
  netcomi_net <- NetCoMi::netConstruct(
    tse_agg,
    taxRank = tax_rank,
    filtTax = "numbSamp",
    filtTaxPar = list(numbSamp = prevalence_threshold),
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
    clustMethod = clust_method,
    hubPar = hub_par,
    normDeg = FALSE
  )

  return(list(
    netcomi_net = netcomi_net,
    netcomi_props = netcomi_props,
    tse_agg = tse_agg
  ))
}
