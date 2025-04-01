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
    pivot_longer(-c(gene_symbol, sample)) |> 
    summarise(tot = sum(value), .by = c("gene_symbol", "sample")) |> 
    summarise(n = median(tot), .by = "sample") |> 
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
#' @param combine_genes If true, combine normalised-to-max epxression across genes
#'
#' @return A ggplot2 object.
plot_one_xy <- function(seu, resolution = NULL, point_size = 0.2, what = c("umap", "spatial"), meta_column = "cluster",
                        genes = NULL, with_legend = FALSE, label_size = 6, combine_genes = FALSE, discrete_palette = NULL) {
  what <- match.arg(what)

  # Convert resolutions (e.g. "0.1") into metadata column name ("Spatial.008um_snn_res.0.1")
  if(!is.null(resolution))
    res <- get_resolution_columns(resolution)

  # Create x-y coordinates for umap or spatial
  if(what == "umap") {
    xy <- FetchData(seu, vars = c("umap_1", "umap_2")) |> 
      as_tibble(rownames = "cell_id") |> 
      rename(x = umap_1, y = umap_2)
    point_shape <- 16
  } else {
    xy <- GetTissueCoordinates(seu) |> 
      as_tibble(rownames = "cell_id") |>
      select(-cell) |> 
      # Rearrange cooridnates to look like Newcastle report
      rename(
        t = x,
        x = y
      ) |> 
      rename(y = t) |> 
      mutate(y = -y) |> 
      relocate(cell_id, .before = 1) 
    point_shape <- 15
  }

  # if gene(s) provided, make an expression plot
  if(!is.null(genes)) {

    ex <- get_gene_expression(seu, genes, norm_per_gene = TRUE)
    if(combine_genes) {
      ex <- ex |> 
        group_by(cell_id) |> 
        summarise(expression = mean(expression)) |> 
        add_column(gene_symbol = "Combined")
    }
    d <- inner_join(xy, ex, by = "cell_id") 
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
   if(length(genes) > 1 & !combine_genes)
      pl <- pl + facet_wrap(~ gene_symbol)
    
    
  # Otherwise make clustering or cell type plot
  } else {

    if(meta_column == "cluster") {
      col <- res
    } else {
      col <- meta_column
    }

    cl <- get_meta_column(seu, col)
    xy_cl <- inner_join(xy, cl, by = "cell_id")
    pl <- xy_cl |>
      ggplot(aes(x = x, y = y)) +
      thm +
      geom_point(aes(colour = value), size = point_size, shape = point_shape, stroke = 0.1) +
      labs(x = NULL, y = NULL)
    if(is.null(discrete_palette)) {
      pl <- pl + scale_colour_brewer(palette = "Paired")
    } else {
      pl <- pl + scale_colour_manual(values = discrete_palette)
    }

    if(what == "umap") {
      xy_cl_mean <- xy_cl |>
        group_by(value) |> 
        summarise(mx = mean(x), my = mean(y))
      pl <- pl +  
        geom_text(data = xy_cl_mean, aes(label = value, x = mx, y = my), colour = "black", size = label_size)
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


plot_cell_type <- function(seu, umap_point_size = 0.2, spatial_point_size = 0.4, label_size = 5) {

  pal <- c(okabe_ito_palette[c(3, 6, 4)], "grey90")
  pl1 <- plot_one_xy(seu, point_size = umap_point_size, what = "umap", meta_column = "cell_type", discrete_palette = pal)
  pl2 <- plot_one_xy(seu, point_size = spatial_point_size, what = "spatial", meta_column = "cell_type", discrete_palette = pal)
  plot_grid(pl1, pl2, nrow = 1, rel_widths = c(1, 1))

}


plot_expression <- function(seu, genes, umap_point_size = 0.1, spatial_point_size = 0.2, label_size = 5, combine_genes = FALSE) {

  plot_side_by_side <- function(gn) {
    pl1 <- plot_one_xy(seu, genes = gn, point_size = umap_point_size, what = "umap", combine_genes = combine_genes)
    pl2 <- plot_one_xy(seu, genes = gn, point_size = spatial_point_size, what = "spatial", combine_genes = combine_genes)
    plot_grid(pl1, pl2, nrow = 1, rel_widths = c(1, 1))
  }

  if(combine_genes) {
    plot_side_by_side(genes)
  } else {
    pl <- map(genes, function(gene) {
      plot_side_by_side(gene)
    })
    plot_grid(plotlist = pl, ncol = 1, labels = genes)
  }

}


plot_gene_expression <- function(seu, gids, layer = "data", width = 0.3, ncol = NULL) {
  
  # Necessary to avoid keeping large (~1 GB) environment within ggplot object
  env <- new.env(parent = globalenv())
  env$width <- width
  env$ncol <- ncol
  
  mt <- seu@meta.data |> 
    as_tibble(rownames = "cell_id")
  
  env$d <- GetAssayData(seu, layer = layer)[gids, , drop = FALSE] |> 
    t() |> 
    as.data.frame() |> 
    rownames_to_column("cell_id") |> 
    as_tibble() |> 
    pivot_longer(-cell_id, names_to = "gene_symbol") |>
    left_join(mt, by = "cell_id") |> 
    mutate(
      grp = cluster,
      igrp = as.integer(grp)
    )
  env$groups <- levels(env$d$grp)
  env$m <- env$d |> 
    group_by(gene_symbol, igrp) |> 
    summarise(
      mn = mean(value[value > 0], na.rm = TRUE),
      md = median(value[value > 0], na.rm = TRUE),
      p = length(which(value > 0)) / n()
    )
  
  mrg <- 0.2
  env$max_y <- max(env$d$value)
  env$min_y <- -mrg * env$max_y

  env$m <- env$m |>
    mutate(y = (1 - p) * env$min_y)
  
  env$d <- env$d |> 
    filter(value > 0)
  
  
  
  env$wid <- 0.3
  with(env, {
    ggplot(d, aes(x = igrp, y = value)) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      geom_quasirandom(size = 0.5, colour = "grey60", width = width) +
      geom_segment(data = m, aes(x = igrp - wid, xend = igrp + wid, y = md, yend = md), linewidth = 1) +
      geom_col(data = m, aes(y = min_y), fill = "brown", width = wid) +
      geom_col(data = m, aes(y = y), fill = "grey80", width = wid) +
      facet_wrap( ~ gene_symbol, ncol = ncol) +
      geom_hline(yintercept = 0, colour = "grey70") +
      scale_x_continuous(breaks = seq_along(groups), labels = groups) +
      scale_y_continuous(limits = c(min_y, max_y), expand = expansion(mult = c(0, 0.03))) +
      labs(x = NULL, y = "Value")
  })
}



plot_gene_heatmap <- function(seu, gns) {
  env <- new.env(parent = globalenv())
  
  if(!is.null(resolution)) {
    res <- get_resolution_columns(resolution)
    Idents(seu) <- res  
  }

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


plot_gene_dots <- function(seu, gns, resolution = NULL) {
  env <- new.env(parent = globalenv())

  if(!is.null(resolution)) {
    res <- get_resolution_columns(resolution)
    Idents(seu) <- res  
  }
  
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


plot_gene_cor <- function(tbl_cor, min_cor = 0.3) {
  good_genes <- tbl_cor |> 
    filter(correlation >= min_cor) |> 
    pull(gene_symbol) |> 
    unique()
  if(length(good_genes) == 0)
    return(NULL)
  tab <- tbl_cor |> 
    filter(gene_symbol %in% good_genes)  |> 
    pivot_wider(id_cols = gene_symbol, names_from = sel_gene, values_from = correlation) |> 
    column_to_rownames("gene_symbol")
  ggheatmap(tab, palette = "viridis", with.y.text = TRUE, legend.name = "Correlation")
}

plot_gene_vs_gene_cluster <- function(seu, genes1, genes2, clust_sel = NULL, point_size = 0.5) {

  mt <- get_cell_data(seu) |> 
    select(cell_id, cluster)

  palette <- okabe_ito_palette
  if(!is.null(clust_sel)) {
    mt <- mt |> 
      mutate(
        cluster = as.character(cluster),
        cluster = if_else(cluster %in% clust_sel, cluster, "0")
      ) 
    palette <- c("grey90", okabe_ito_palette)
  }

  d1 <- get_gene_expression(seu, genes1, norm_per_gene = TRUE) |> 
    summarise(gr1 = mean(expression), .by = cell_id)
  d2 <- get_gene_expression(seu, genes2, norm_per_gene = TRUE) |> 
    summarise(gr2 = mean(expression), .by = cell_id)

  inner_join(d1, d2, by = "cell_id") |> 
    left_join(mt, by = "cell_id") |> 
  ggplot(aes(x = gr1, y = gr2, colour = cluster)) +
    th +
    geom_point(size = point_size) +
    scale_colour_manual(values = palette) +
    labs(x = str_c(genes1, collapse = ","), y = str_c(genes2, collapse = ","))
}