#' Default Oral Microbiota Genera
#'
#' Returns a character vector of genera commonly associated with oral microbiota
#' and relevant to CRC research via the oral-gut axis. Based on literature
#' consensus from multiple CRC microbiome studies.
#'
#' @return Character vector of oral microbiota genus names.
#'
#' @export
default_oral_genera <- function() {
  c("Fusobacterium", "Porphyromonas", "Prevotella",
    "Peptostreptococcus", "Parvimonas", "Streptococcus", "Gemella",
    "Haemophilus", "Neisseria", "Veillonella", "Actinomyces", "Rothia",
    "Capnocytophaga", "Granulicatella", "Treponema", "Campylobacter",
    "Leptotrichia", "Selenomonas")
}


#' Default Gut Anaerobic Commensal Genera
#'
#' Returns a character vector of genera commonly associated with healthy gut
#' anaerobic commensals (SCFA producers, fiber fermenters).
#'
#' @return Character vector of commensal genus names.
#'
#' @export
default_commensal_genera <- function() {
  c("Faecalibacterium", "Roseburia", "Eubacterium",
    "Ruminococcus", "Coprococcus", "Blautia", "Dorea", "Anaerostipes",
    "Butyrivibrio", "Subdoligranulum", "Lachnospira", "Bifidobacterium",
    "Akkermansia", "Alistipes", "Bacteroides", "Parabacteroides")
}


#' Clean Taxon Names for Display
#'
#' Removes SILVA database suffixes and brackets from taxon names for cleaner
#' plot labels. Handles collisions by appending numeric suffixes.
#'
#' @param names Character vector of taxon names.
#'
#' @return Character vector of cleaned display names.
#'
#' @export
clean_taxon_names <- function(names) {
  cleaned <- names
  # Remove brackets: [Eubacterium] -> Eubacterium, (Ruminococcus) -> Ruminococcus
  cleaned <- gsub("^\\[(.+?)\\]", "\\1", cleaned)
  cleaned <- gsub("^\\((.+?)\\)", "\\1", cleaned)
  # Remove SILVA suffixes: _UCG_NNN, _UCG-NNN, _NK4ANNN_group, _AD3011_group
  cleaned <- gsub("_UCG[_-]\\d+", "", cleaned)
  cleaned <- gsub("_NK\\w+", "", cleaned)
  cleaned <- gsub("_AD\\d+", "", cleaned)
  cleaned <- gsub("_R-7", "", cleaned)
  # Remove trailing _group
  cleaned <- gsub("_group$", "", cleaned)
  # Remove _NC2004 and similar numeric codes
  cleaned <- gsub("_NC\\d+", "", cleaned)
  # Clean up trailing underscores/hyphens
  cleaned <- gsub("[_-]+$", "", cleaned)
  # Replace remaining underscores with spaces for display
  cleaned <- gsub("_", " ", cleaned)
  # Trim whitespace
  cleaned <- trimws(cleaned)

  # Deduplicate: append suffix for collisions
  dups <- duplicated(cleaned) | duplicated(cleaned, fromLast = TRUE)
  if (any(dups)) {
    dup_names <- unique(cleaned[dups])
    for (dn in dup_names) {
      idx <- which(cleaned == dn)
      cleaned[idx] <- paste0(cleaned[idx], " (", seq_along(idx), ")")
    }
  }

  cleaned
}


# Node attributes shared by the association and differential graphs: taxonomy colour group
# (rare levels lumped into "Other" when max_other_levels > 0), mean CLR abundance for node
# size, annotation group, display name and the taxonomy columns.
.network_node_data <- function(node_names, coords, tse_agg, color_tax_rank, max_other_levels,
                               annotation_groups, min_group_nodes, clean_names) {
  tax_data <- as.data.frame(SummarizedExperiment::rowData(tse_agg))[node_names, , drop = FALSE]

  color_col <- if (color_tax_rank %in% colnames(tax_data)) {
    vals <- as.character(tax_data[[color_tax_rank]])
    vals[is.na(vals) | grepl("Incertae_Sedis|Unknown|uncultured", vals, ignore.case = TRUE)] <- "Unclassified"
    vals
  } else {
    rep("Unknown", length(node_names))
  }
  color_col <- lump_rare_levels(list(color_col), max_other_levels)[[1]]

  # Mean CLR (pseudocount 1) over all samples of tse_agg, shifted to be positive for sizing.
  mclr <- tryCatch({
    clr <- SummarizedExperiment::assay(
      mia::transformAssay(tse_agg, assay.type = "counts", method = "clr", pseudocount = 1), "clr")
    rowMeans(clr)[node_names]
  }, error = function(e) {
    message("Could not compute CLR, using log counts: ", e$message)
    log1p(rowMeans(SummarizedExperiment::assay(tse_agg, "counts")))[node_names]
  })
  mclr <- mclr - min(mclr, na.rm = TRUE) + 0.1

  # Annotation groups (first match wins; groups with < min_group_nodes members are not drawn)
  annotation_group <- rep(NA_character_, length(node_names))
  for (grp_name in names(annotation_groups)) {
    grp_match <- node_names %in% annotation_groups[[grp_name]] & is.na(annotation_group)
    if (sum(grp_match) >= min_group_nodes) annotation_group[grp_match] <- grp_name
  }

  node_data <- tibble::tibble(
    name = node_names,
    display_name = if (clean_names) clean_taxon_names(node_names) else node_names,
    x = coords[, 1],
    y = coords[, 2],
    color_tax = color_col,
    mclr = as.numeric(mclr),
    annotation_group = annotation_group
  )
  for (col in intersect(c("kingdom", "phylum", "class", "order", "family", "genus"), colnames(tax_data))) {
    node_data[[col]] <- as.character(tax_data[[col]])
  }
  node_data
}

#' Lump taxonomy levels with fewer than `min_n` nodes (counted across all vectors) into "Other"
#'
#' @param values_list List of character vectors (one per network panel).
#' @param min_n Minimum count to keep a level; 0 disables lumping.
#' @return The list with rare levels replaced by "Other" ("Unclassified"/"Unknown" are kept).
#' @keywords internal
lump_rare_levels <- function(values_list, min_n = 2) {
  counts <- table(unlist(lapply(values_list, as.character)))
  rare <- setdiff(names(counts)[counts < min_n], c("Unclassified", "Unknown"))
  lapply(values_list, function(v) {
    v <- as.character(v)
    v[v %in% rare] <- "Other"
    v
  })
}

#' Qualitative palette for taxonomy levels (Set2, Set3, then interpolated Set3), with fixed
#' greys for "Other" and "Unclassified".
#' @keywords internal
network_level_palette <- function(levels) {
  levels <- sort(unique(levels))
  n <- length(levels)
  cols <- if (n <= 8) {
    RColorBrewer::brewer.pal(max(3, n), "Set2")[seq_len(n)]
  } else if (n <= 12) {
    RColorBrewer::brewer.pal(max(3, n), "Set3")[seq_len(n)]
  } else {
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n)
  }
  pal <- stats::setNames(cols, levels)
  pal[names(pal) == "Other"] <- "#B0B0B0"
  pal[names(pal) == "Unclassified"] <- "#C0C0C0"
  pal
}

# Re-level color_tax on several graph-data lists to one shared palette and rebuild their graphs.
.apply_shared_palette <- function(gd_list) {
  palette <- network_level_palette(unlist(lapply(gd_list, function(g) as.character(g$node_data$color_tax))))
  gd_list <- lapply(gd_list, function(gd) {
    gd$node_data$color_tax <- factor(gd$node_data$color_tax, levels = names(palette))
    gd$graph <- tidygraph::tbl_graph(nodes = gd$node_data, edges = gd$edge_data, directed = FALSE)
    gd
  })
  list(graph_data = gd_list, palette = palette)
}


#' Extract NetCoMi Results into ggraph-Compatible Data
#'
#' Converts NetCoMi network objects into a tidygraph/igraph structure with
#' taxonomy metadata, annotation groups, and layout coordinates for ggraph plotting.
#' Removes singleton nodes (no edges) by default.
#'
#' @param netcomi_result A list from calc_netcomi_single() containing netcomi_net,
#'   netcomi_props, and tse_agg.
#' @param color_tax_rank Character, taxonomic rank for node coloring (default: "order").
#' @param annotation_groups Named list of character vectors defining annotation groups
#'   for polygon highlighting. Names are group labels, values are genus name vectors.
#'   Example: list("Oral microbiota" = c("Fusobacterium", "Prevotella"),
#'                 "Anaerobic commensals" = c("Faecalibacterium", "Roseburia"))
#'   (default: NULL, no annotations).
#' @param min_group_nodes Integer, minimum nodes in an annotation group to draw hull (default: 2).
#' @param max_other_levels Integer, taxonomy levels with fewer than this many genera
#'   are grouped as "Other". Set to 0 to disable (default: 2).
#' @param remove_singletons Logical, remove nodes with no edges (default: TRUE).
#' @param layout_coords Optional matrix of pre-computed layout coordinates (2 columns: x, y).
#'   Row names must match node names. If provided, layout is not recomputed.
#' @param network_index Integer, which network's properties to use from netcomi_props
#'   (1 or 2). Relevant for comparison networks where netcomi_props contains
#'   properties for both networks (default: 1).
#' @param clean_names Logical, whether to create cleaned display names by removing
#'   SILVA suffixes and brackets (default: TRUE).
#' @param seed Numeric, random seed for layout (default: 13075).
#'
#' @return A list containing:
#'   - graph: tidygraph tbl_graph object
#'   - node_data: tibble with node attributes (includes display_name if clean_names=TRUE)
#'   - edge_data: tibble with edge attributes
#'   - layout_coords: matrix of x, y coordinates with node names as rownames
#'
#' @export
extract_netcomi_graph_data <- function(netcomi_result,
                                       color_tax_rank = "order",
                                       annotation_groups = NULL,
                                       min_group_nodes = 2,
                                       max_other_levels = 2,
                                       remove_singletons = TRUE,
                                       layout_coords = NULL,
                                       network_index = 1,
                                       clean_names = TRUE,
                                       seed = 13075) {

  netcomi_net <- netcomi_result$netcomi_net
  netcomi_props <- netcomi_result$netcomi_props
  tse_agg <- netcomi_result$tse_agg

  adj_mat <- netcomi_net$adjaMat1
  # assoMat1 has the original signed association values;
  # adjaMat1 is non-negative after dissimilarity transformation
  asso_mat <- netcomi_net$assoMat1

  # --- 1. Build igraph from adjacency matrix ---
  ig <- igraph::graph_from_adjacency_matrix(
    adj_mat,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  # Determine edge sign from the signed association matrix
  edge_list <- igraph::as_edgelist(ig)
  if (!is.null(asso_mat) && nrow(edge_list) > 0) {
    edge_signs <- vapply(seq_len(nrow(edge_list)), function(i) {
      val <- asso_mat[edge_list[i, 1], edge_list[i, 2]]
      if (val >= 0) "Positive" else "Negative"
    }, character(1))
  } else {
    edge_signs <- rep("Positive", nrow(edge_list))
  }
  igraph::E(ig)$sign <- edge_signs

  # --- 2. Remove singletons (nodes with degree 0) ---
  if (remove_singletons) {
    deg <- igraph::degree(ig)
    keep <- deg > 0
    if (sum(keep) > 0) {
      ig <- igraph::induced_subgraph(ig, which(keep))
    }
  }

  node_names <- igraph::V(ig)$name

  # --- 3. Compute layout ---
  if (!is.null(layout_coords)) {
    available <- intersect(node_names, rownames(layout_coords))
    coords <- matrix(NA_real_, nrow = length(node_names), ncol = 2,
                     dimnames = list(node_names, c("x", "y")))
    coords[available, ] <- layout_coords[available, , drop = FALSE]
    missing <- setdiff(node_names, available)
    if (length(missing) > 0) {
      set.seed(seed)
      cx <- mean(coords[available, 1], na.rm = TRUE)
      cy <- mean(coords[available, 2], na.rm = TRUE)
      coords[missing, ] <- cbind(
        stats::rnorm(length(missing), mean = cx, sd = 0.5),
        stats::rnorm(length(missing), mean = cy, sd = 0.5)
      )
    }
  } else {
    set.seed(seed)
    coords <- igraph::layout_with_fr(ig, niter = 1000)
    rownames(coords) <- node_names
    colnames(coords) <- c("x", "y")
  }

  # --- 4. Node attributes ---
  node_data <- .network_node_data(node_names, coords, tse_agg, color_tax_rank, max_other_levels,
                                  annotation_groups, min_group_nodes, clean_names)

  # Hub scores and clusters of the requested network (comparison objects hold two)
  props_vec <- function(slot, name) {
    v <- netcomi_props[[slot]][[paste0(name, network_index)]]
    if (is.null(v)) rep(NA, length(node_names)) else v[node_names]
  }
  node_data$hub_score <- as.numeric(props_vec("centralities", "eigenvecConnect"))
  node_data$cluster <- as.integer(props_vec("clustering", "clust"))

  # --- 5. Assemble edge data ---
  edge_df <- igraph::as_data_frame(ig, what = "edges")
  edge_data <- tibble::tibble(
    from = edge_df$from,
    to = edge_df$to,
    weight = edge_df$weight,
    sign = edge_df$sign
  )

  # --- 6. Build tidygraph ---
  tbl_graph <- tidygraph::tbl_graph(
    nodes = node_data,
    edges = edge_data,
    directed = FALSE
  )

  return(list(
    graph = tbl_graph,
    node_data = node_data,
    edge_data = edge_data,
    layout_coords = coords
  ))
}


#' Compute Union Layout for Multiple Networks
#'
#' Computes a shared Fruchterman-Reingold layout from the union of multiple
#' adjacency matrices, ensuring consistent node positions across panels.
#' Only includes nodes that have at least one edge in any network.
#'
#' @param netcomi_results A list of results from calc_netcomi_single().
#' @param seed Numeric, random seed (default: 13075).
#'
#' @return A matrix with x, y columns and node names as rownames.
#'
#' @export
compute_union_layout <- function(netcomi_results, seed = 13075) {
  adj_mats <- lapply(netcomi_results, function(res) res$netcomi_net$adjaMat1)

  all_names <- unique(unlist(lapply(adj_mats, rownames)))

  # Build union matrix (element-wise max of absolute values)
  union_mat <- matrix(0, nrow = length(all_names), ncol = length(all_names),
                      dimnames = list(all_names, all_names))

  for (mat in adj_mats) {
    shared <- intersect(rownames(mat), all_names)
    union_mat[shared, shared] <- pmax(union_mat[shared, shared], abs(mat[shared, shared]))
  }

  # Remove singletons from union
  row_sums <- rowSums(union_mat)
  connected <- row_sums > 0
  if (sum(connected) > 0) {
    union_mat <- union_mat[connected, connected]
    all_names <- all_names[connected]
  }

  union_ig <- igraph::graph_from_adjacency_matrix(
    union_mat, mode = "undirected", weighted = TRUE, diag = FALSE
  )

  set.seed(seed)
  coords <- igraph::layout_with_fr(union_ig, niter = 1000)
  rownames(coords) <- all_names
  colnames(coords) <- c("x", "y")

  coords
}


#' Extract Graph Data from a NetCoMi Comparison Result
#'
#' Convenience wrapper that properly extracts ggraph-compatible data for both
#' networks from a calc_netcomi_comparison() result, using the correct network
#' index for centralities and clustering, and computing a shared union layout.
#'
#' @param netcomi_comparison A list from calc_netcomi_comparison().
#' @param color_tax_rank Character, taxonomic rank for node coloring (default: "order").
#' @param annotation_groups Named list of character vectors for polygon annotations (default: NULL).
#' @param clean_names Logical, whether to clean taxon names (default: TRUE).
#' @param seed Numeric, random seed (default: 13075).
#'
#' @return A list containing:
#'   - group1_data: graph data for group 1 (from extract_netcomi_graph_data)
#'   - group2_data: graph data for group 2
#'   - union_coords: shared layout coordinates
#'   - shared_palette: named color vector for consistent coloring across panels
#'   - group1: name of group 1
#'   - group2: name of group 2
#'
#' @export
extract_comparison_graph_data <- function(netcomi_comparison,
                                          color_tax_rank = "order",
                                          annotation_groups = NULL,
                                          clean_names = TRUE,
                                          seed = 13075) {

  net <- netcomi_comparison$netcomi_net
  props <- netcomi_comparison$netcomi_props

  # Unpack into single-network format for each group
  result1 <- list(
    netcomi_net = list(adjaMat1 = net$adjaMat1, assoMat1 = net$assoMat1),
    netcomi_props = props,
    tse_agg = netcomi_comparison$tse_agg
  )
  result2 <- list(
    netcomi_net = list(adjaMat1 = net$adjaMat2, assoMat1 = net$assoMat2),
    netcomi_props = props,
    tse_agg = netcomi_comparison$tse_agg
  )

  # Compute shared union layout
  union_coords <- compute_union_layout(
    list(g1 = result1, g2 = result2),
    seed = seed
  )

  # Disable per-network "Other" grouping; will be done globally below
  gd1 <- extract_netcomi_graph_data(
    result1,
    color_tax_rank = color_tax_rank,
    annotation_groups = annotation_groups,
    max_other_levels = 0,
    layout_coords = union_coords,
    network_index = 1,
    clean_names = clean_names,
    seed = seed
  )

  gd2 <- extract_netcomi_graph_data(
    result2,
    color_tax_rank = color_tax_rank,
    annotation_groups = annotation_groups,
    max_other_levels = 0,
    layout_coords = union_coords,
    network_index = 2,
    clean_names = clean_names,
    seed = seed
  )

  # Global "Other" grouping (levels with < 2 genera across BOTH networks), then one palette
  lumped <- lump_rare_levels(list(gd1$node_data$color_tax, gd2$node_data$color_tax))
  gd1$node_data$color_tax <- lumped[[1]]
  gd2$node_data$color_tax <- lumped[[2]]
  shared <- .apply_shared_palette(list(gd1, gd2))
  gd1 <- shared$graph_data[[1]]
  gd2 <- shared$graph_data[[2]]
  shared_palette <- shared$palette

  list(
    group1_data = gd1,
    group2_data = gd2,
    union_coords = union_coords,
    shared_palette = shared_palette,
    group1 = netcomi_comparison$group1,
    group2 = netcomi_comparison$group2
  )
}


#' Extract Differential Network Graph Data for ggraph
#'
#' Converts a NetCoMi diffnet result into ggraph-compatible data. Edges represent
#' significantly different associations between two groups (Fisher's z-test).
#' Edge signs indicate which group has the stronger association.
#'
#' @param netcomi_comparison A list from calc_netcomi_comparison() with compute_diffnet = TRUE.
#' @param color_tax_rank Character, taxonomic rank for node coloring (default: "order").
#' @param annotation_groups Named list of character vectors for polygon annotations (default: NULL).
#' @param min_group_nodes Integer, minimum annotated nodes to form a polygon (default: 2).
#' @param max_other_levels Integer, minimum genera for a taxonomy level to be shown
#'   individually; levels with fewer genera become "Other" (default: 2).
#' @param clean_names Logical, whether to clean taxon names for display (default: TRUE).
#' @param seed Numeric, random seed for reproducibility (default: 13075).
#'
#' @return A list with the same structure as extract_netcomi_graph_data():
#'   - graph: tidygraph object
#'   - node_data: tibble with node attributes
#'   - edge_data: tibble with edge attributes (sign = "Stronger in {group}" direction)
#'   - layout_coords: matrix of x, y coordinates
#'
#' @export
extract_diffnet_graph_data <- function(netcomi_comparison,
                                       color_tax_rank = "order",
                                       annotation_groups = NULL,
                                       min_group_nodes = 2,
                                       max_other_levels = 2,
                                       clean_names = TRUE,
                                       seed = 13075) {

  diffnet_obj <- netcomi_comparison$netcomi_diffnet
  if (is.null(diffnet_obj)) {
    stop("No diffnet result found. Run calc_netcomi_comparison() with compute_diffnet = TRUE.")
  }

  net <- netcomi_comparison$netcomi_net
  tse_agg <- netcomi_comparison$tse_agg
  group1 <- netcomi_comparison$group1
  group2 <- netcomi_comparison$group2

  diff_adj <- diffnet_obj$diffAdjustMat

  # --- 1. Build igraph from differential adjacency matrix ---
  # diffAdjustMat is binary (0/1) after adjustment — use association difference as weight
  asso_diff <- abs(net$assoMat1 - net$assoMat2)

  # Zero out non-significant entries
  weighted_diff <- asso_diff * (diff_adj != 0)

  ig <- igraph::graph_from_adjacency_matrix(
    weighted_diff,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  # Determine edge direction: which group has the stronger association
  edge_list <- igraph::as_edgelist(ig)
  if (nrow(edge_list) > 0) {
    edge_signs <- vapply(seq_len(nrow(edge_list)), function(i) {
      v1 <- net$assoMat1[edge_list[i, 1], edge_list[i, 2]]
      v2 <- net$assoMat2[edge_list[i, 1], edge_list[i, 2]]
      if (abs(v1) >= abs(v2)) {
        paste0("Stronger in ", group1)
      } else {
        paste0("Stronger in ", group2)
      }
    }, character(1))
    igraph::E(ig)$sign <- edge_signs
  }

  # --- 2. Remove singletons ---
  deg <- igraph::degree(ig)
  keep <- deg > 0
  if (sum(keep) == 0) {
    # Return empty graph data
    return(list(
      graph = tidygraph::tbl_graph(
        nodes = tibble::tibble(name = character(0)),
        edges = tibble::tibble(from = integer(0), to = integer(0)),
        directed = FALSE
      ),
      node_data = tibble::tibble(name = character(0)),
      edge_data = tibble::tibble(from = character(0), to = character(0),
                                 weight = numeric(0), sign = character(0)),
      layout_coords = matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("x", "y")))
    ))
  }
  ig <- igraph::induced_subgraph(ig, which(keep))
  node_names <- igraph::V(ig)$name

  # --- 3. Compute layout ---
  set.seed(seed)
  coords <- igraph::layout_with_fr(ig, niter = 1000)
  rownames(coords) <- node_names
  colnames(coords) <- c("x", "y")

  # --- 4. Node attributes ---
  node_data <- .network_node_data(node_names, coords, tse_agg, color_tax_rank, max_other_levels,
                                  annotation_groups, min_group_nodes, clean_names)

  # --- 5. Assemble edge data ---
  edge_df <- igraph::as_data_frame(ig, what = "edges")
  edge_data <- tibble::tibble(
    from = edge_df$from,
    to = edge_df$to,
    weight = edge_df$weight,
    sign = edge_df$sign
  )

  # --- 6. Build tidygraph ---
  tbl_graph <- tidygraph::tbl_graph(
    nodes = node_data,
    edges = edge_data,
    directed = FALSE
  )

  return(list(
    graph = tbl_graph,
    node_data = node_data,
    edge_data = edge_data,
    layout_coords = coords
  ))
}
