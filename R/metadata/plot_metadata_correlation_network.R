#' Plot Metadata Correlation Network
#'
#' Calculates and plots correlations between *numeric* metadata variables as a network graph.
#'
#' @param tse A TreeSummarizedExperiment object.
#' @param variables Character vector of *numeric* variable names from colData(tse).
#' @param pretty_labels Optional named character vector to rename variables for node labels.
#' @param method Correlation method (passed to Hmisc::rcorr, default "spearman").
#' @param p_adjust_method P-value adjustment method (passed to p.adjust, default "BH").
#' @param p_threshold Significance threshold for filtering edges (default 0.05).
#' @param cor_threshold Minimum absolute correlation strength for filtering edges (default 0.2).
#' @param layout Network layout algorithm for ggraph (default "fr").
#' @param seed Optional random seed for layout reproducibility.
#' @param network_title Title for the network plot (default: NULL).
#' @param base_fontsize Base font size for network elements (default: 10).
#' @param ... Additional arguments passed to ggraph theme.
#'
#' @return A ggraph object.
#'
#' @examples
#' \dontrun{
#' network_plot <- plot_metadata_correlation_network(
#'   tse = tse_object,
#'   variables = c("Age", "BMI", "DiseaseScore"),
#'   method = "spearman",
#'   p_threshold = 0.05,
#'   cor_threshold = 0.3
#' )
#' print(network_plot)
#' }
#'
#' @export
plot_metadata_correlation_network <- function(tse,
                                              variables,
                                              pretty_labels = NULL,
                                              method = "spearman",
                                              p_adjust_method = "BH",
                                              p_threshold = 0.05,
                                              cor_threshold = 0.2,
                                              layout = "fr",
                                              seed = NULL,
                                              network_title = NULL,
                                              base_fontsize = 10,
                                              ...) {

  # --- Data Preparation ---
  # Variables must already be numeric: as.numeric() on a factor silently returns its level
  # codes, which would then be rank-correlated as if they were ordered quantities.
  metadata_subset <- as.data.frame(SummarizedExperiment::colData(tse)[, variables, drop = FALSE])
  non_numeric <- variables[!vapply(metadata_subset, is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("Non-numeric variable(s): ", paste(non_numeric, collapse = ", "),
         ". Encode them explicitly before calling.")
  }

  # --- Calculate Correlations ---
  corr_results <- Hmisc::rcorr(as.matrix(metadata_subset), type = method)

  # --- Create Edge List ---
  # One row per unordered pair, with r, p and p_adj all read from the same (i, j) index.
  idx <- which(lower.tri(corr_results$r), arr.ind = TRUE)
  edge_list <- data.frame(
    from = rownames(corr_results$r)[idx[, "row"]],
    to = colnames(corr_results$r)[idx[, "col"]],
    correlation = corr_results$r[idx],
    p = corr_results$P[idx]
  )
  edge_list$p_adj <- stats::p.adjust(edge_list$p, method = p_adjust_method)

  # Filter edges based on significance and correlation strength
  filtered_edges <- edge_list |>
    dplyr::filter(!is.na(p_adj) & p_adj < p_threshold & abs(correlation) >= cor_threshold) |>
    dplyr::mutate(
      sign = ifelse(correlation > 0, "Positive", "Negative"),
      abs_cor = abs(correlation)
    )

  # --- Create Node List (only nodes with edges) ---
  if (nrow(filtered_edges) > 0) {
    connected_nodes <- unique(c(as.character(filtered_edges$from), as.character(filtered_edges$to)))
    nodes <- data.frame(name = connected_nodes)
  } else {
    nodes <- data.frame(name = character(0))
  }

  # --- Create Graph Object ---
  graph_obj <- igraph::graph_from_data_frame(d = filtered_edges, vertices = nodes, directed = FALSE)

  # --- Plot Network ---
  if (!is.null(seed)) {
    set.seed(seed)
  }

  final_network_title <- if (is.null(network_title)) "Metadata Correlation Network" else network_title

  # Use pretty labels if provided
  if (nrow(nodes) > 0) {
    node_labels <- if (!is.null(pretty_labels)) pretty_labels[nodes$name] else nodes$name
    node_labels[is.na(node_labels)] <- nodes$name[is.na(node_labels)]
  } else {
    node_labels <- character(0)
  }

  # Handle case with no significant edges
  if (nrow(filtered_edges) == 0) {
    warning(paste("No significant correlations found (p.adj <", p_threshold,
      ", |", method, "| >=", cor_threshold, "). Plotting nodes only."))
    network_ggraph <- ggraph::ggraph(graph_obj, layout = layout) +
      ggraph::geom_node_point(size = 4, color = "grey50") +
      ggraph::geom_node_text(ggplot2::aes(label = node_labels), repel = TRUE,
        size = 3, max.overlaps = 15) +
      ggraph::theme_graph(base_family = "sans", base_size = base_fontsize) +
      ggplot2::labs(title = NULL, subtitle = NULL) +
      ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_blank(),
        plot.subtitle = ggplot2::element_blank())
  } else {
    network_ggraph <- ggraph::ggraph(graph_obj, layout = layout) +
      ggraph::geom_edge_link(ggplot2::aes(edge_width = abs_cor, color = sign),
        alpha = 0.7, show.legend = TRUE) +
      ggraph::scale_edge_width(range = c(0.5, 3), name = "Abs. Correlation") +
      ggraph::scale_edge_color_manual(values = c("Positive" = "firebrick", "Negative" = "steelblue"),
        name = "Correlation Sign") +
      ggraph::geom_node_point(size = 4, color = "grey30") +
      ggraph::geom_node_text(ggplot2::aes(label = node_labels), repel = TRUE,
        size = 3, max.overlaps = 15) +
      ggraph::theme_graph(base_family = "sans", base_size = base_fontsize) +
      ggplot2::labs(title = NULL, subtitle = NULL) +
      ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_blank(),
        plot.subtitle = ggplot2::element_blank())
  }

  return(network_ggraph)
}
