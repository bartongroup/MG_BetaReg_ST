okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#AAAAAA", "#000000")

th <- ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank()
  )

thm <- th +
  ggplot2::theme(
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  )

gs <- function(gg, name, sample = NULL, width, height) {
  if(is.null(sample)) {
    fname <- file.path("fig", paste0(name, ".png"))
  } else {
    fname <- file.path("fig", paste0(name, "_", sample, ".png"))
  }
  ggsave(filename = fname, plot = gg, device = "png",
         width = width, height = height, dpi = 300)
}

scientific_10 <- function(x) {
  parse(
    text = scales::scientific_format()(x) |>
      str_replace("e\\+*", " %*% 10^") |> 
      str_replace("^1 \\%\\*\\% 10\\^", "10^")
  )
}





make_table_cuts <- function(params) {
  tribble(
    ~Group, ~Quantity, ~Limit,
    "Per cell", "Number of detected genes per cell", str_c(">", params$min_detected_genes_per_cell),
    "Per cell", "Number of counts per cell", str_c(">", params$min_counts_per_cell),
    "Per cell", "Mitochondrial count fraction", str_c("<", params$max_mt_fraction),
    "Per gene", "Number of cells a gene is detected in", str_c(">", params$min_detections_per_gene)
  )
}



plot_cell_marginals <- function(cell_data, params, with_vlines = TRUE) {
  pl_ <- function(what, xlab, log_scale = TRUE) {
    pl <- cell_data |> 
      mutate(val = get(what)) |> 
    ggplot(aes(x = val)) +
      th +
      theme(legend.position = "none") +
      #geom_density_ridges() +
      geom_density(fill = "grey70") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
      #scale_fill_manual(values = okabe_ito_palette) +
      labs(x = xlab, y = NULL)
    if(log_scale)
      pl <- pl + scale_x_log10(labels = scales::comma)
    pl
  }
  
  pl1 <- pl_("total_count", "Total count per cell")
  if(with_vlines) pl1 <- pl1 + geom_vline(xintercept = params$min_counts_per_cell, linetype = "dotted")
  pl2 <- pl_("n_genes", "Genes detected per cell")
  if(with_vlines) pl2 <-pl2 + geom_vline(xintercept = params$min_detected_genes_per_cell, linetype = "dotted")
  pl3 <- pl_("mito_fraction", "Mitochondrial fraction", log_scale = FALSE) +
  if(with_vlines) pl3 <- pl3 + geom_vline(xintercept = params$max_mt_fraction, linetype = "dotted")
  
  #cowplot::plot_grid(pl1, pl2, pl3, nrow = 1, align = "h", labels = c("A", "B", "C"))
  
  # cowplot throws cryptic errors while in targets; using patchwork instead
  pl1 + pl2 + pl3
  
}


plot_gene_detections <- function(gene_detection, max_det = 8) {
  gene_detection |> 
    pivot_longer(-gene_id, names_to = "group", values_to = "detections") |> 
    group_by(group, detections) |> 
    tally() |> 
    mutate(
      det = as.character(detections),
      det = if_else(detections > max_det, sprintf(">%d", max_det), det) |> as_factor()
    ) |> 
    group_by(group, det) |>
    summarise(tot = sum(n)) |> 
  ggplot(aes(x = det, y = tot, fill = group)) +
    th +
    geom_col() +
    facet_wrap(~group) +
    scale_fill_manual(values = okabe_ito_palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Number of detections", y = "Number of genes")
  
}

plot_count_det <- function(cell_data, params, with_lines = TRUE) {
  pl <- cell_data |> 
    ggplot(aes(x = total_count, y = n_genes, colour = mito_fraction)) +
    th +
    geom_point(size = 0.5) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    scale_colour_viridis_c(option = "cividis") +
    labs(x = "Total count per cell", y = "Genes detected per cell", colour = "MT fraction")
  if(with_lines) {
    pl <- pl +
      geom_hline(yintercept = params$min_detected_genes_per_cell, linetype = "dotted") +
      geom_vline(xintercept = params$min_counts_per_cell, linetype = "dotted")
  }
    
  return(pl)
}

plot_all_counts <- function(cell_data_all, gene_detection_all) {
  pl_ <- function(d, ylab) {
    ggplot(d, aes(x = sample, y = n)) +
    th +
    geom_col() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Sample", y = ylab)
  }

  pl1 <- cell_data_all |> 
    count(sample) |> 
    pl_("Cell count")

  pl3 <- gene_detection_all |> 
    summarise(n = median(blank), .by = "sample") |> 
    pl_("Median gene count per cell")

  pl2 <- cell_data_all |> 
    summarise(n = sum(total_count) / 1e6, .by = "sample") |> 
    pl_("Total read count (millions)")


  plot_grid(pl1, pl2, pl3, nrow = 1)
}

#'
#' Generates a scatter plot using either UMAP or spatial coordinates from a Seurat object.
#' Points can be colored by cluster identity or gene expression values.
#'
#' @param seu A Seurat object containing the data to plot.
#' @param resolution Character or numeric value indicating the clustering resolution.
#' @param point_size Numeric value specifying the size of points in the plot. Default is 0.2.
#' @param what Character value specifying the type of coordinates to plot. Choices are "umap" or "spatial". Default is "umap".
#' @param gene Optional character specifying the gene names whose expression will be used for coloring. Default is NULL.
#' @param with_legend Logical indicating whether to include the legend in the plot. Default is FALSE.
#'
#' @return A ggplot2 object.
plot_one_xy <- function(seu, resolution = NULL, point_size = 0.2, what = c("umap", "spatial"),
                        genes = NULL, with_legend = FALSE, label_size = 6) {
  what <- match.arg(what)

  # Convert resolutions (e.g. "0.1") into metadata column name ("Spatial.008um_snn_res.0.1")
  if(!is.null(resolution))
    res <- get_resolution_columns(resolution)

  # Create x-y coordinates for umap or spatial
  if(what == "umap") {
    xy <- FetchData(seu, vars = c("umap_1", "umap_2")) |> 
      as_tibble(rownames = "cell") |> 
      rename(x = umap_1, y = umap_2)
    point_shape <- 16
  } else {
    xy <- GetTissueCoordinates(seu) |> 
      as_tibble() |> 
      # Rearrange cooridnates to look like Newcastle report
      rename(
        t = x,
        x = y
      ) |> 
      rename(y = t) |> 
      mutate(y = -y) |> 
      relocate(cell, .before = 1) 
    point_shape <- 15
  }

  # if gene(s) provided, make an expression plot
  if(!is.null(genes)) {

    ex <- get_gene_expression(seu, genes, norm_per_gene = TRUE)
    d <- inner_join(xy, ex, by = "cell") 
    # Plot zeroes first, then expressed genes on top
    d0 <- d |> filter(expression == 0)
    d1 <- d |> filter(expression > 0)

    pl <- d0 |> 
      ggplot(aes(x = x, y = y, colour = expression)) +
        thm +
        geom_point(data = d0, size = point_size, shape = point_shape, stroke = 0.1) +
        geom_point(data = d1, size = point_size, shape = point_shape, stroke = 0.1) +
        scale_colour_distiller(palette = "Blues", direction = 1) +
        labs(x = NULL, y = NULL)
   if(length(genes) > 1)
      pl <- pl + facet_wrap(~ gene_symbol)
    
    
  # Otherwise make clustering plot
  } else {

    cl <- get_cluster_cells(seu, res)
    xy_cl <- inner_join(xy, cl, by = "cell")
    pl <- xy_cl |>
      ggplot(aes(x = x, y = y)) +
      thm +
      geom_point(aes(colour = cluster), size = point_size, shape = point_shape, stroke = 0.1) +
      scale_colour_brewer(palette = "Paired") +
      labs(x = NULL, y = NULL)

    if(what == "umap") {
      xy_cl_mean <- xy_cl |>
        group_by(cluster) |> 
        summarise(mx = mean(x), my = mean(y))
      pl <- pl +  
        geom_text(data = xy_cl_mean, aes(label = cluster, x = mx, y = my), colour = "black", size = label_size)
    }

  }

  if(!with_legend)
      pl <- pl + theme(legend.position = "none")
  
  return(pl)
}


plot_resolutions <- function(seu, resolutions, umap_point_size = 0.2, spatial_point_size = 0.4, label_size = 5) {

  pl <- map(resolutions, function(res) {

    pl1 <- plot_one_xy(seu, resolution = res, point_size = umap_point_size, what = "umap")
    pl2 <- plot_one_xy(seu, resolution = res, point_size = spatial_point_size, what = "spatial")
    plot_grid(pl1, pl2, nrow = 1, rel_widths = c(1, 1))
  })

  plot_grid(plotlist = pl, ncol = 1, labels = resolutions)

}


plot_expression <- function(seu, genes, umap_point_size = 0.1, spatial_point_size = 0.2, label_size = 5) {

  pl <- map(genes, function(gene) {

    pl1 <- plot_one_xy(seu, genes = gene, point_size = umap_point_size, what = "umap")
    pl2 <- plot_one_xy(seu, genes = gene, point_size = spatial_point_size, what = "spatial")
    plot_grid(pl1, pl2, nrow = 1, rel_widths = c(1, 1))
  })

  plot_grid(plotlist = pl, ncol = 1, labels = genes)

}


plot_expression_vs_clusters <- function(seu, genes, resolution, point_size = 0.5) {
  res <- get_resolution_columns(resolution)
  ex <- get_gene_expression(seu, genes) |> 
    filter(expression > 0)
  cl <- get_cluster_cells(seu, res)
  
  inner_join(ex, cl, by = "cell") |> 
    ggplot(aes(x = cluster, y = expression)) + 
    th +
    geom_quasirandom(size = point_size) +
    facet_wrap(~gene_symbol)
}




plot_top_cluster_markers_heatmap <- function(seu, top_cluster_markers, resolution) {
  env <- new.env(parent = globalenv())
  
  gns <- unique(top_cluster_markers$gene_symbol)

  res <- get_resolution_columns(resolution)
  Idents(seu) <- res
  
  d <- AverageExpression(seu, layer = "data", features = gns, group.by = "ident") |> 
    pluck("Spatial.008um") |> 
    as_tibble(rownames = "gene_symbol") |> 
    pivot_longer(-gene_symbol, names_to = "cluster") |> 
    group_by(gene_symbol) |> 
    mutate(
      z = value / mean(value)
    ) |> 
    ungroup()
  
  env$tab <- d |>
    pivot_wider(id_cols = gene_symbol, names_from = cluster, values_from = z) |>
    column_to_rownames("gene_symbol")

  with(env, {
    ggheatmap(tab, order.col = FALSE, with.y.text = TRUE, palette = "viridis",
              legend.name = "rel expr")  
  })
}


plot_top_cluster_markers_dots <- function(seu, top_cluster_markers, resolution, extra_genes = NULL) {
  env <- new.env(parent = globalenv())
  
  gns <- unique(top_cluster_markers$gene_symbol)
  
  if(!is.null(extra_genes))
    gns <- c(gns, extra_genes)
  
  res <- get_resolution_columns(resolution)
  Idents(seu) <- res
  
  
  dp <- DotPlot(seu, features = gns)
  env$d <- dp$data |> 
    as_tibble(rownames = "gene_symbol")
   
  with(env, {
    ggplot(d, aes(x = id, y = features.plot, colour = avg.exp.scaled, size = pct.exp)) +
      th +
      geom_point() +
      scale_colour_gradient(low = "grey90", high = "blue") +
      scale_radius(range = c(1, 6)) +
      labs(x = "Cluster", y = "Gene", size = "Percent expressed", colour = "Average expresion")
  })
}


plot_png_expression_genes <- function(seu, genes, file_name) {
  plot_expression(seu, genes, spatial_point_size = 0.5, umap_point_size = 0.5) |> 
    gs(file_name, width = 8, height = 4 * length(genes))
}


plot_key_markers_umap_spatial <- function(seu, key_markers, file_prefix) {
  all_feat <- Features(seu)

  km <- key_markers |> 
    filter(gene_symbol %in% all_feat)
  nums <- key_markers |> 
    count(cell_type)

  map(nums$cell_type, function(ct) {
    gname <- str_glue("{file_prefix}_{ct}")
    gns <- km |> 
      filter(cell_type == ct) |> 
      pull(gene_symbol)
    plot_png_expression_genes(seu, gns, gname)
  }) 
}

plot_top_cluster_umap_spatial <- function(seu, top_cluster_markers, file_prefix = "top_clust_gene") {
  
  clusters <- top_cluster_markers$cluster |> 
    unique() |> 
    sort()
  
  map(clusters, function(cl) {
    gname <- str_glue("{file_prefix}_cl{cl}")
    gns <- top_cluster_markers |> 
      filter(cluster == cl) |> 
      pull(gene_symbol)
    plot_png_expression_genes(seu, gns, gname)
  })
}
