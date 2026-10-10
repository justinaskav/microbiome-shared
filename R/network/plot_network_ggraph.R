#' Default Annotation Group Colors
#'
#' Returns a named character vector of colors for annotation groups. Colors are
#' assigned deterministically by hashing the group name so that the same group
#' always gets the same color regardless of which other groups are present.
#'
#' @param group_names Character vector of annotation group names.
#' @return Named character vector of colors.
#' @export
default_annotation_colors <- function(group_names) {
  ann_pal <- c("#E8927C", "#7CB5E8", "#A8D5A2", "#D4A0D4",
               "#F0C987", "#87CEEB", "#FFB6C1", "#B0E0E6")
  unique_names <- unique(group_names)
  # Use character code sum mod palette length for stable per-name assignment
  indices <- vapply(unique_names, function(nm) {
    (sum(utf8ToInt(nm)) %% length(ann_pal)) + 1L
  }, integer(1))
  # Resolve collisions by incrementing
  used <- integer(0)
  for (i in seq_along(indices)) {
    while (indices[i] %in% used) {
      indices[i] <- (indices[i] %% length(ann_pal)) + 1L
    }
    used <- c(used, indices[i])
  }
  setNames(ann_pal[indices], unique_names)
}


#' Plot Microbial Network with ggraph
#'
#' Creates a ggraph-based network visualization with polygon hull annotations
#' for oral microbiota and gut anaerobic commensals. Nodes are colored by a
#' configurable taxonomic rank and sized by mean CLR abundance.
#'
#' Distant nodes within the same annotation group are split into separate
#' sub-clusters (using hierarchical clustering on coordinates) so that
#' polygons remain compact.
#'
#' @param graph_data A list from extract_netcomi_graph_data().
#' @param color_tax_rank Character, label for the color legend (default: "Order").
#' @param show_polygons Logical, whether to show annotation hull polygons (default: TRUE).
#' @param polygon_alpha Numeric, fill alpha for polygons (default: 0.10).
#' @param polygon_expand Numeric, hull expansion in mm (default: 2).
#' @param polygon_concavity Numeric, concavity parameter for geom_mark_hull (default: 3).
#'   Lower values produce more convex hulls; higher values allow more concave shapes.
#' @param polygon_radius Numeric, corner rounding radius in mm (default: 2).
#' @param polygon_max_dist Numeric, maximum Euclidean distance between nodes in the
#'   same polygon sub-cluster; nodes farther apart are split into separate polygons
#'   (default: NULL, auto-computed as 18% of layout diameter).
#' @param node_size_range Numeric vector of length 2 for node size scale (default: c(2, 7)).
#' @param edge_width_range Numeric vector of length 2 for edge width scale (default: c(0.3, 1.5)).
#' @param edge_alpha Numeric, edge transparency (default: 0.6).
#' @param adaptive_edge_alpha Logical, whether to auto-scale edge alpha based on
#'   edge count (default: FALSE). When TRUE, overrides edge_alpha.
#' @param label_size Numeric, genus label text size (default: 2.5).
#' @param show_labels Logical, whether to show node labels (default: TRUE).
#' @param label_annotated_fill Logical, whether to add colored backgrounds to labels
#'   of annotated genera (oral/commensal) matching polygon colors (default: TRUE).
#' @param label_fill_alpha Numeric, alpha for label background fill (default: 0.4).
#' @param max_overlaps Integer, max label overlaps for ggrepel (default: 20).
#' @param color_palette Optional named character vector of colors for taxonomy levels.
#' @param annotation_colors Optional named character vector of fill colors for annotation
#'   groups. Names must match the annotation group labels used in extract_netcomi_graph_data().
#'   If NULL, auto-generated from a qualitative palette (default: NULL).
#' @param edge_type Character, type of edge geom: "link" for straight lines,
#'   "arc" for curved arcs (default: "link").
#' @param edge_strength Numeric, curvature strength for arc edges (default: 0.3).
#' @param edge_colors Named character vector for edge sign colors (default: NULL,
#'   uses c("Positive" = "#2E8B57", "Negative" = "#CD4F39")). Override for
#'   differential networks, e.g. c("Stronger in Ileum" = "#4A8DC8", "Stronger in Lesion" = "#8B5CB8").
#' @param title Character, plot title (default: NULL).
#' @param base_size Numeric, base font size (default: 10).
#'
#' @return A ggplot object.
#'
#' @export
plot_network_ggraph <- function(graph_data,
                                color_tax_rank = "Order",
                                show_polygons = TRUE,
                                polygon_alpha = 0.10,
                                polygon_expand = 2,
                                polygon_concavity = 3,
                                polygon_radius = 2,
                                polygon_max_dist = NULL,
                                node_size_range = c(2, 7),
                                edge_width_range = c(0.3, 1.5),
                                edge_alpha = 0.6,
                                adaptive_edge_alpha = FALSE,
                                label_size = 2.5,
                                show_labels = TRUE,
                                label_annotated_fill = TRUE,
                                label_fill_alpha = 0.4,
                                max_overlaps = 20,
                                color_palette = NULL,
                                annotation_colors = NULL,
                                edge_type = "link",
                                edge_strength = 0.3,
                                edge_colors = NULL,
                                title = NULL,
                                base_size = 10) {

  node_data <- graph_data$node_data
  graph <- graph_data$graph

  # Determine label column: use display_name if available, else name
  label_col <- if ("display_name" %in% colnames(node_data)) "display_name" else "name"

  # Adaptive edge alpha based on edge count
  if (adaptive_edge_alpha) {
    n_edges <- nrow(graph_data$edge_data)
    edge_alpha <- max(0.12, min(0.5, 0.7 - 0.25 * log10(max(n_edges, 1))))
  }

  # Create manual layout
  layout_df <- ggraph::create_layout(graph, layout = "manual", x = node_data$x, y = node_data$y)

  if (is.null(color_palette)) color_palette <- network_level_palette(as.character(node_data$color_tax))

  # Annotation group colors — deterministic by group name for consistency across plots
  ann_groups_present <- unique(na.omit(node_data$annotation_group))
  if (is.null(annotation_colors) && length(ann_groups_present) > 0) {
    annotation_colors <- default_annotation_colors(ann_groups_present)
  }

  # --- Prepare polygon data with sub-clustering for distant nodes ---
  has_annotations <- show_polygons && any(!is.na(layout_df$annotation_group))
  polygon_data <- NULL

  if (has_annotations) {
    ann_rows <- layout_df[!is.na(layout_df$annotation_group), , drop = FALSE]

    # Auto-compute max distance as 18% of layout diameter
    if (is.null(polygon_max_dist)) {
      x_range <- diff(range(layout_df$x, na.rm = TRUE))
      y_range <- diff(range(layout_df$y, na.rm = TRUE))
      polygon_max_dist <- 0.18 * sqrt(x_range^2 + y_range^2)
    }

    # Split each annotation group into sub-clusters
    polygon_rows <- list()
    for (grp in unique(ann_rows$annotation_group)) {
      grp_data <- ann_rows[ann_rows$annotation_group == grp, , drop = FALSE]
      if (nrow(grp_data) < 2) {
        grp_data$polygon_id <- paste0(grp, "_1")
        polygon_rows[[length(polygon_rows) + 1]] <- grp_data
        next
      }
      # Hierarchical clustering on spatial coordinates
      d <- stats::dist(grp_data[, c("x", "y")])
      hc <- stats::hclust(d, method = "complete")
      sub_clusters <- stats::cutree(hc, h = polygon_max_dist)
      grp_data$polygon_id <- paste0(grp, "_", sub_clusters)

      # Only keep sub-clusters with >= 2 nodes (can't draw hull around 1 point)
      cluster_sizes <- table(grp_data$polygon_id)
      grp_data <- grp_data[grp_data$polygon_id %in% names(cluster_sizes[cluster_sizes >= 2]), , drop = FALSE]
      if (nrow(grp_data) > 0) {
        polygon_rows[[length(polygon_rows) + 1]] <- grp_data
      }
    }

    if (length(polygon_rows) > 0) {
      polygon_data <- do.call(rbind, polygon_rows)
    }
  }

  # --- Build plot ---
  p <- ggraph::ggraph(layout_df)

  # Layer 1: Hull polygons (fill aesthetic)
  if (show_polygons && !is.null(polygon_data) && nrow(polygon_data) > 0) {
    p <- p +
      ggforce::geom_mark_hull(
        data = polygon_data,
        ggplot2::aes(x = x, y = y, fill = annotation_group, group = polygon_id),
        alpha = polygon_alpha,
        expand = ggplot2::unit(polygon_expand, "mm"),
        concavity = polygon_concavity,
        radius = ggplot2::unit(polygon_radius, "mm"),
        show.legend = TRUE
      ) +
      ggplot2::scale_fill_manual(
        values = annotation_colors,
        name = "Annotation",
        drop = FALSE
      )
  }

  # Layer 2: Edges
  if (edge_type == "arc") {
    p <- p +
      ggraph::geom_edge_arc(
        ggplot2::aes(edge_width = weight, edge_colour = sign),
        alpha = edge_alpha,
        strength = edge_strength,
        show.legend = TRUE
      )
  } else {
    p <- p +
      ggraph::geom_edge_link(
        ggplot2::aes(edge_width = weight, edge_colour = sign),
        alpha = edge_alpha,
        show.legend = TRUE
      )
  }
  p <- p +
    ggraph::scale_edge_width_continuous(
      range = edge_width_range,
      name = "Association\nstrength",
      guide = "none"
    ) +
    ggraph::scale_edge_color_manual(
      values = if (!is.null(edge_colors)) edge_colors else c("Positive" = "#2E8B57", "Negative" = "#CD4F39"),
      name = "Association"
    )

  # Layer 3: Nodes (use colour aesthetic to avoid fill conflict with polygons)
  layout_df$color_tax <- factor(layout_df$color_tax, levels = names(color_palette))

  p <- p +
    ggraph::geom_node_point(
      ggplot2::aes(size = mclr, colour = color_tax),
      shape = 16
    ) +
    ggplot2::scale_size_continuous(
      range = node_size_range,
      name = "Mean CLR\nabundance",
      guide = "none"
    ) +
    ggplot2::scale_colour_manual(
      values = color_palette,
      name = color_tax_rank,
      drop = FALSE,
      guide = ggplot2::guide_legend(override.aes = list(size = 3))
    )

  # Layer 4: Labels — split by annotation group (colored bg) and unannotated (white halo)
  if (show_labels) {
    has_annotation_col <- "annotation_group" %in% colnames(layout_df)
    unannotated_idx <- if (has_annotation_col) which(is.na(layout_df$annotation_group)) else seq_len(nrow(layout_df))

    # Add a label layer per annotation group with matching background color
    if (label_annotated_fill && has_annotation_col && !is.null(annotation_colors)) {
      for (grp_name in names(annotation_colors)) {
        grp_idx <- which(layout_df$annotation_group == grp_name)
        if (length(grp_idx) > 0) {
          grp_data <- layout_df[grp_idx, , drop = FALSE]
          p <- p +
            ggraph::geom_node_text(
              data = grp_data,
              ggplot2::aes(x = x, y = y, label = .data[[label_col]]),
              repel = TRUE,
              size = label_size,
              max.overlaps = max_overlaps,
              segment.color = "grey60",
              segment.size = 0.2,
              bg.colour = scales::alpha(annotation_colors[[grp_name]], label_fill_alpha),
              bg.r = 0.15
            )
        }
      }
    }

    if (length(unannotated_idx) > 0 || !label_annotated_fill) {
      # Unannotated nodes (or all nodes if label_annotated_fill is FALSE)
      if (label_annotated_fill) {
        unannotated_data <- layout_df[unannotated_idx, , drop = FALSE]
      } else {
        unannotated_data <- layout_df
      }
      p <- p +
        ggraph::geom_node_text(
          data = unannotated_data,
          ggplot2::aes(x = x, y = y, label = .data[[label_col]]),
          repel = TRUE,
          size = label_size,
          max.overlaps = max_overlaps,
          segment.color = "grey60",
          segment.size = 0.2,
          bg.colour = scales::alpha("white", 0.7),
          bg.r = 0.08
        )
    }
  }

  # Theme
  p <- p +
    ggraph::theme_graph(base_family = "sans", base_size = base_size) +
    ggplot2::theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.key.size = ggplot2::unit(0.4, "cm"),
      legend.text = ggplot2::element_text(size = base_size - 2),
      legend.title = ggplot2::element_text(size = base_size - 1),
      plot.title = ggplot2::element_text(size = base_size + 2, face = "bold"),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )

  if (!is.null(title)) {
    p <- p + ggplot2::labs(title = title)
  }

  p
}


#' Create Multi-Panel Network Figure
#'
#' Constructs separate co-occurrence networks for each group level and combines
#' them into a multi-panel figure with shared layout coordinates.
#'
#' @param tse A TreeSummarizedExperiment object.
#' @param group_var Character, column name in colData for grouping.
#' @param group_levels Character vector of group levels to include.
#' @param tax_rank Character, taxonomic rank for network (default: "genus").
#' @param color_tax_rank Character, taxonomic rank for node coloring (default: "order").
#' @param annotation_groups Named list of character vectors for polygon annotations (default: NULL).
#' @param prevalence_threshold Numeric, prevalence filter (default: 0.3).
#' @param spars_method Character, sparsification method (default: "none").
#' @param measure Character, association measure (default: "spring").
#' @param measure_par List of measure parameters.
#' @param min_samples Integer, minimum samples per group (default: 10).
#' @param shared_layout Logical, use shared union layout (default: TRUE).
#' @param ncol Integer, number of columns in panel (default: 2).
#' @param label_size Numeric, text label size (default: 2.2).
#' @param annotation_colors Optional named fill colours for annotation groups.
#' @param seed Numeric, random seed (default: 13075).
#' @param ... Additional arguments passed to plot_network_ggraph().
#'
#' @return A list containing:
#'   - plot: patchwork object
#'   - netcomi_results: list of calc_netcomi_single results per group
#'   - graph_data: list of extract_netcomi_graph_data results per group
#'
#' @export
create_network_panel <- function(tse,
                                 group_var,
                                 group_levels,
                                 tax_rank = "genus",
                                 color_tax_rank = "order",
                                 annotation_groups = NULL,
                                 annotation_colors = NULL,
                                 prevalence_threshold = 0.3,
                                 spars_method = "none",
                                 measure = "spring",
                                 measure_par = list(thresh = 0.05, Rmethod = "approx"),
                                 min_samples = 10,
                                 shared_layout = TRUE,
                                 ncol = 2,
                                 label_size = 2.2,
                                 seed = 13075,
                                 ...) {

  cd <- SummarizedExperiment::colData(tse)
  if (!group_var %in% colnames(cd)) {
    stop("group_var '", group_var, "' not found in colData")
  }

  # --- 1. Compute networks per group ---
  netcomi_results <- list()
  for (grp in group_levels) {
    grp_idx <- as.character(cd[[group_var]]) == grp
    n_samples <- sum(grp_idx)

    if (n_samples < min_samples) {
      message("Skipping group '", grp, "': only ", n_samples, " samples (need ", min_samples, ")")
      next
    }

    message("\n--- Computing network for: ", grp, " (n=", n_samples, ") ---")
    tse_grp <- tse[, grp_idx]

    netcomi_results[[grp]] <- calc_netcomi_single(
      tse_grp,
      tax_rank = tax_rank,
      prevalence_threshold = prevalence_threshold,
      measure = measure,
      measure_par = measure_par,
      spars_method = spars_method,
      seed = seed
    )
  }

  if (length(netcomi_results) == 0) {
    stop("No groups had sufficient samples for network construction.")
  }

  # --- 2. Compute shared layout ---
  union_coords <- if (shared_layout && length(netcomi_results) > 1) {
    compute_union_layout(netcomi_results, seed = seed)
  } else {
    NULL
  }

  # --- 3. Extract graph data (disable per-panel "Other" grouping) ---
  graph_data_list <- lapply(netcomi_results, extract_netcomi_graph_data,
                            color_tax_rank = color_tax_rank,
                            annotation_groups = annotation_groups,
                            max_other_levels = 0,
                            layout_coords = union_coords,
                            seed = seed)

  # --- 4. Global "Other" grouping (levels with < 2 genera across ALL panels), one palette ---
  lumped <- lump_rare_levels(lapply(graph_data_list, function(gd) gd$node_data$color_tax))
  for (grp in names(graph_data_list)) graph_data_list[[grp]]$node_data$color_tax <- lumped[[grp]]
  shared <- .apply_shared_palette(graph_data_list)
  graph_data_list <- shared$graph_data
  shared_palette <- shared$palette

  # --- 5. Build individual plots ---
  # Use provided annotation_colors or compute from annotation_groups definition
  shared_annotation_colors <- if (!is.null(annotation_colors)) {
    annotation_colors
  } else if (!is.null(annotation_groups)) {
    default_annotation_colors(names(annotation_groups))
  } else {
    NULL
  }

  plots <- list()
  color_label <- tools::toTitleCase(color_tax_rank)

  for (grp in names(graph_data_list)) {
    n_samples <- sum(as.character(cd[[group_var]]) == grp)
    n_edges <- nrow(graph_data_list[[grp]]$edge_data)
    plot_title <- paste0(grp, " (n=", n_samples, ", edges=", n_edges, ")")

    plots[[grp]] <- plot_network_ggraph(
      graph_data_list[[grp]],
      color_tax_rank = color_label,
      color_palette = shared_palette,
      annotation_colors = shared_annotation_colors,
      title = plot_title,
      label_size = label_size,
      ...
    )
  }

  # --- 6. Combine with patchwork and shared legend ---
  # Use the panel with most Order levels as the legend source
  legend_idx <- which.max(sapply(graph_data_list, function(g) nrow(g$node_data)))
  legend_plot <- plots[[legend_idx]]

  # Remove legends from all panels
  plots_no_legend <- lapply(plots, function(p) {
    p + ggplot2::theme(legend.position = "none")
  })

  # Build grid: panels on the left, shared legend on the right
  legend_grob <- cowplot::get_legend(legend_plot)

  combined <- cowplot::plot_grid(
    patchwork::wrap_plots(plots_no_legend, ncol = ncol),
    legend_grob,
    ncol = 2,
    rel_widths = c(1, 0.22)
  )

  return(list(
    plot = combined,
    netcomi_results = netcomi_results,
    graph_data = graph_data_list
  ))
}


#' Extract Network Topology Summary
#'
#' Creates a summary table of network topology metrics across groups. `Taxa` counts every taxon
#' that entered the network; `Nodes` only those with at least one edge (singletons are removed
#' from the plotted graph). Density is computed over all taxa, so it is comparable across groups
#' with different numbers of isolated taxa.
#'
#' @param panel_result A list from create_network_panel().
#'
#' @return A tibble with topology metrics per group.
#'
#' @export
extract_network_topology <- function(panel_result) {
  metrics <- list()

  for (grp in names(panel_result$netcomi_results)) {
    props <- panel_result$netcomi_results[[grp]]$netcomi_props
    gd <- panel_result$graph_data[[grp]]

    n_taxa <- nrow(panel_result$netcomi_results[[grp]]$netcomi_net$adjaMat1)
    n_nodes <- nrow(gd$node_data)
    n_edges <- nrow(gd$edge_data)
    n_pos <- sum(gd$edge_data$sign == "Positive")
    n_neg <- sum(gd$edge_data$sign == "Negative")
    density <- if (n_taxa > 1) 2 * n_edges / (n_taxa * (n_taxa - 1)) else 0

    modularity <- tryCatch(
      as.numeric(props$clustering$modularity1),
      error = function(e) NA_real_
    )
    if (is.null(modularity) || length(modularity) == 0) modularity <- NA_real_

    avg_path <- tryCatch(
      as.numeric(props$globalProps$avPath1),
      error = function(e) NA_real_
    )
    if (is.null(avg_path) || length(avg_path) == 0) avg_path <- NA_real_

    metrics[[grp]] <- tibble::tibble(
      Group = grp,
      Taxa = n_taxa,
      Nodes = n_nodes,
      Edges = n_edges,
      `Positive edges` = n_pos,
      `Negative edges` = n_neg,
      Density = round(density, 3),
      Modularity = round(modularity, 3),
      `Avg. path length` = round(avg_path, 2)
    )
  }

  dplyr::bind_rows(metrics)
}
