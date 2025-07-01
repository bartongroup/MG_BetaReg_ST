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

gs <- function(gg, name, sample = NULL, width, height, dpi = 240) {
  if(is.null(sample)) {
    fname <- file.path("doc", "images", paste0(name, ".png"))
  } else {
    fname <- file.path("doc", "images", paste0(name, "_", sample, ".png"))
  }

  ggsave(filename = fname, plot = gg, device = "png", limitsize = FALSE,
         width = width, height = height, dpi = dpi)
}

rs <- function(gg, name, sample = NULL) {
  if(is.null(sample)) {
    fname <- file.path("doc", "images", paste0(name, "rds"))
  } else {
    fname <- file.path("doc", "images", paste0(name, "_", sample, ".rds"))
  }

  readr::write_rds(gg, fname)
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

plot_cell_marginals <- function(cell_data, params = NULL, with_vlines = TRUE) {
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

plot_count_det <- function(cell_data, params = NULL, with_lines = TRUE) {
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
  env <- new.env(parent = globalenv())

  env$d1 <- cell_data_all |> 
    count(sample)
  env$d3 <- gene_detection_all |> 
    pivot_longer(-c(gene_symbol, sample)) |> 
    summarise(tot = sum(value), .by = c("gene_symbol", "sample")) |> 
    summarise(n = median(tot), .by = "sample")
  env$d2 <- cell_data_all |> 
    summarise(n = sum(total_count) / 1e6, .by = "sample")

  with(env, {
    pl_ <- function(d, ylab) {
      ggplot(d, aes(x = sample, y = n)) +
      th +
      geom_col() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(x = "Sample", y = ylab)
    }

    pl1 <- pl_(d1, "Cell count")
    pl3 <- pl_(d3, "Median gene count per cell")
    pl2 <- pl_(d2, "Total read count (millions)")
  
    plot_grid(pl1, pl2, pl3, nrow = 1)
  })
}

plot_spatial_meta <- function(seu, meta_column = "pseudotime", point_size = 0.6, palette = "viridis") {
  d <- get_spatial_coordinates(seu) |> 
    mutate(val = seu@meta.data[[meta_column]])
  mx <- max(abs(d$val))
  pl <-  d |> 
  ggplot(aes(x = x, y = y, colour = val)) +
    thm +
    geom_point(size = point_size, shape = 15, stroke = 0.1) +
    labs(x = NULL, y = NULL, colour = meta_column)
  if(palette == "viridis") {
    pl <- pl + scale_colour_viridis_c(option = "cividis")
  } else if(palette == "divergent") {
    pl <- pl + scale_colour_distiller(type = "div", palette = "RdBu", limits = c(-mx, mx))
  }
  pl
}

percentiles_no_zeroes <- function(df, percentile) {
  df |> 
    mutate(
      prc = if_else(expression == 0, NA, expression),
      prc = rank(prc, na.last = "keep") / length(na.omit(prc)),
      prc = replace_na(prc, 0), 
      val = case_when(
        prc == 0 ~ 0,
        prc < percentile ~ 1,
        prc >= percentile ~ 2
      ) |> as.factor()
    )
}

plot_spatial_selection_genes <- function(seu, genes, box = list(x1 = -Inf, x2 = Inf, y1 = -Inf, y2 = Inf), percentile = NULL, spots = NULL, ncol = 2) {
  ex <- get_gene_expression(seu, genes) |> 
    mutate(val = expression)
  if(!is.null(percentile)) {
    # percentiles, but ignore zeroes
    ex <- ex |> 
      group_by(gene_symbol) |>
      percentiles_no_zeroes(percentile) |> 
      ungroup()
  }

  xy <- get_spatial_coordinates(seu) |> 
    filter(x >= box$x1 & x <= box$x2 & y >= box$y1 & y <= box$y2)
  corners <- get_resolution_grid(xy)
  d <- inner_join(corners, ex, by = "cell_id", relationship = "many-to-many")

  pl <- d |> 
    ggplot() +
    thm +
    geom_polygon(aes(x = x_corner, y = y_corner, group = cell_id, fill = val), colour = NA, linewidth = 0) +
    labs(x = NULL, y = NULL, fill = "Expression") +
    facet_wrap(~gene_symbol, ncol = ncol)

  if(is.null(percentile)) {
    pl <- pl + scale_fill_viridis_c(option = "cividis")
  } else {
    pl <- pl + scale_fill_manual(values = c("beige", "darkseagreen", "black")) +
      theme(legend.position = "none")
  }

  r <- 20
  a <- 2 * pi * seq(1, 30) / 29
  if(!is.null(spots)) {
    sp <- spots |> 
      mutate(idx = row_number()) |> 
      group_split(idx) |> 
      map(function(w) {
        tibble(
          idx = w$idx,
          fx = w$x + r * cos(a),
          fy = w$y + r * sin(a)
        )
      }) |> 
      list_rbind()
    pl <- pl +
      geom_path(data = sp, aes(x = fx, y = fy, group = idx), colour = "red", linewidth = 0.5, alpha = 0.7)
  }
  pl
}

animate_spatial_genes <- function(seu, genes, box = list(x1 = -Inf, x2 = Inf, y1 = -Inf, y2 = Inf), percentile = NULL,
                                      transition_length = 1) {
  ex <- get_gene_expression(seu, genes) |> 
    mutate(val = expression)
  if(!is.null(percentile)) {
    # percentiles, but ignore zeroes
    ex <- ex |> 
      group_by(gene_symbol) |>
      percentiles_no_zeroes(percentile) |> 
      ungroup()
  }

  xy <- get_spatial_coordinates(seu) |> 
    filter(x >= box$x1 & x <= box$x2 & y >= box$y1 & y <= box$y2)
  corners <- get_resolution_grid(xy)
  d <- inner_join(corners, ex, by = "cell_id", relationship = "many-to-many") |> 
    mutate(gene_symbol = factor(gene_symbol, levels = genes))

  pl <- d |> 
    ggplot(aes(x = x_corner, y = y_corner, group = cell_id, fill = val)) +
    thm +
    geom_polygon(colour = NA, linewidth = 0) +
    labs(x = NULL, y = NULL, fill = "Expression") +
    theme(legend.position = "none") +
    labs(title = "{closest_state}") +
    gganimate::transition_states(gene_symbol, transition_length = transition_length)

    if(is.null(percentile)) {
    pl <- pl + scale_fill_viridis_c(option = "cividis")
  } else {
    pl <- pl + scale_fill_manual(values = c("beige", "darkseagreen", "black")) +
      theme(legend.position = "none")
  }

  sgenes <- str_c(genes, collapse = "_")
  fname <- str_glue("spatial_{sgenes}.gif")
  gganimate::anim_save(fname, pl, path = file.path("doc", "images"))
}

  



plot_spatial_selection_meta <- function(seu, meta_column = "cell_type", box = list(x1 = -Inf, x2 = Inf, y1 = -Inf, y2 = Inf)) {
  mt <- get_cell_data(seu) |> 
    mutate(val = seu@meta.data[[meta_column]]) |> 
    select(cell_id, val)
  xy <- get_spatial_coordinates(seu) |> 
    filter(x >= box$x1 & x <= box$x2 & y >= box$y1 & y <= box$y2)
  corners <- get_resolution_grid(xy)
  d <- inner_join(corners, mt, by = "cell_id", relationship = "many-to-many") 

  pal <- c(okabe_ito_palette[c(4, 1, 8, 5, 2, 3, 6, 7, 9)], "darkgreen", "lightgreen", "red", "brown", "magenta", grey(0.95))
  d |> 
    ggplot(aes(x = x_corner, y = y_corner, group = cell_id, fill = val)) +
    thm +
    geom_polygon(colour = NA, linewidth = 0) +
    labs(x = NULL, y = NULL, fill = "Expression") +
    scale_fill_manual(values = pal)
}


plot_reduction_meta <- function(seu, reduction = "pca", meta_column = "pseudotime", point_size = 0.6) {
  Embeddings(seu, reduction = reduction) |>
    as_tibble(rownames = "cell_id") |>
    select(1:3) |> 
    set_names(c("cell_id", "x", "y")) |> 
    mutate(val = seu@meta.data[[meta_column]]) |> 
  ggplot(aes(x = x, y = y, colour = val)) +
    thm +
    geom_point(size = point_size, shape = 15, stroke = 0.1) +
    scale_colour_viridis_c(option = "cividis") +
    labs(x = NULL, y = NULL, colour = meta_column)
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
    xy <- get_spatial_coordinates(seu)
    corners <- get_resolution_grid(xy)
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

    if(what == "umap") {
      d <- inner_join(xy, ex, by = "cell_id") 
    # Plot zeroes first, then expressed genes on top
      d0 <- d |> filter(expression == 0)
      d1 <- d |> filter(expression > 0)
      pl <- d |> 
        ggplot(aes(x = x, y = y, colour = expression)) +
          thm +
          geom_point(data = d0, size = point_size, shape = point_shape, stroke = 0.1) +
          geom_point(data = d1, size = point_size, shape = point_shape, stroke = 0.1) +
          scale_colour_distiller(palette = "Purples", direction = 1) +
          labs(x = NULL, y = NULL)
    } else {
      d <- inner_join(corners, ex, by = "cell_id") 
      pl <- d |> 
        ggplot(aes(x = x_corner, y = y_corner, group = cell_id, fill = expression)) +
          thm +
          geom_polygon(colour = NA, linewidth = 0) +
          scale_fill_distiller(palette = "Purples", direction = 1) +
          labs(x = NULL, y = NULL)  
    }


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
    

    if(what == "umap") {
      xy_cl <- inner_join(xy, cl, by = "cell_id")
      pl <- xy_cl |>
      ggplot(aes(x = x, y = y)) +
      thm +
      geom_point(aes(colour = value), size = point_size, shape = point_shape, stroke = 0.1) +
      guides(color = guide_legend(override.aes = list(size = 3))) +
      labs(x = NULL, y = NULL)
      if(is.null(discrete_palette)) {
        pl <- pl + scale_colour_brewer(palette = "Paired")
      } else {
        pl <- pl + scale_colour_manual(values = discrete_palette)
      }
      xy_cl_mean <- xy_cl |>
        group_by(value) |> 
        summarise(mx = mean(x), my = mean(y))
      pl <- pl +  
        geom_text(data = xy_cl_mean, aes(label = value, x = mx, y = my), colour = "black", size = label_size)

    } else {
      xy_cl <- inner_join(corners, cl, by = "cell_id")
      pl <- xy_cl |>
        ggplot(aes(x = x_corner, y = y_corner, group = centre_id, fill = value)) +
        thm +
        geom_polygon(colour = NA, linewidth = 0) +
        guides(fill = guide_legend(override.aes = list(size = 3))) +
        labs(x = NULL, y = NULL)  
        if(is.null(discrete_palette)) {
          pl <- pl + scale_fill_brewer(palette = "Paired")
        } else {
          pl <- pl + scale_fill_manual(values = discrete_palette)
        }
  
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


plot_cell_type <- function(seu, umap_point_size = 0.3, spatial_point_size = 0.6, label_size = 3) {

  pal <- c(okabe_ito_palette[c(4, 1, 8, 5, 2, 3, 6, 7, 9)], "darkgreen", "lightgreen", "red", "brown", "magenta", grey(0.95))
  pl1 <- plot_one_xy(seu, point_size = umap_point_size, label_size = label_size, what = "umap", meta_column = "cell_type", discrete_palette = pal, with_legend = FALSE)
  pl2 <- plot_one_xy(seu, point_size = spatial_point_size, label_size = label_size, what = "spatial", meta_column = "cell_type", discrete_palette = pal, with_legend = TRUE) + guides(colour = guide_legend(override.aes = list(size = 3)))
  plot_grid(pl1, pl2, nrow = 1, rel_widths = c(1, 1.4))

}


plot_expression <- function(seu, genes, umap_point_size = 0.1, spatial_point_size = 0.2, label_size = 5, grid_label_size = 10, combine_genes = FALSE) {

  plot_side_by_side <- function(gn) {
    pl1 <- plot_one_xy(seu, genes = gn, point_size = umap_point_size, what = "umap", combine_genes = combine_genes, label_size = label_size)
    pl2 <- plot_one_xy(seu, genes = gn, point_size = spatial_point_size, what = "spatial", combine_genes = combine_genes, label_size = label_size)
    plot_grid(pl1, pl2, nrow = 1, rel_widths = c(1, 1))
  }

  if(combine_genes) {
    plot_side_by_side(genes)
  } else {
    pl <- map(genes, function(gene) {
      plot_side_by_side(gene)
    })
    plot_grid(plotlist = pl, ncol = 1, labels = genes, label_size = grid_label_size, label_y = 0.95)
  }

}


plot_gene_expression <- function(seu, gids, layer = "data", meta_column = "cluster", width = 0.3, ncol = NULL) {
  
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
      grp = get(meta_column),
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
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      ) +
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



plot_gene_heatmap <- function(seu, gns, resolution = NULL) {
  env <- new.env(parent = globalenv())
  
  if(!is.null(resolution)) {
    res <- get_resolution_columns(resolution)
    Idents(seu) <- res  
  }

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


plot_gene_dots <- function(seu, gns, meta_column = "cluster", with_legend = TRUE, title = NULL) {
  env <- new.env(parent = globalenv())
  env$meta_column <- meta_column
  env$with_legend <- with_legend
  env$title <- title

  Idents(seu) <- meta_column
  
  dp <- DotPlot(seu, features = gns)
  env$d <- dp$data |> 
    as_tibble(rownames = "gene_symbol")
   
  with(env, {
    pl <- ggplot(d, aes(x = id, y = fct_rev(features.plot), colour = avg.exp.scaled, size = pct.exp)) +
      th +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      ) +
      geom_point() +
      scale_colour_gradient(low = "grey90", high = "blue") +
      scale_radius(range = c(1, 6)) +
      labs(x = meta_column, y = "Gene", size = "Percent expressed", colour = "Average expresion", title = title)
    if(!with_legend)
      pl <- pl + theme(legend.position = "none")
    pl
  })
}

plot_cell_type_marker_dots <- function(seu, cell_type_markers, n_top = 30) {
  cell_types <- cell_type_markers$cluster |> unique() |> as.character()
  pl <- map(cell_types, function(ct) {
    gns <- cell_type_markers |> 
      filter(cluster == ct) |> 
      slice_max(avg_log2FC, n = n_top) |> 
      arrange(desc(avg_log2FC)) |> 
      pull(gene_symbol)
    plot_gene_dots(seu, gns, meta_column = "cell_type", with_legend = FALSE, title = ct)
  })
  plot_grid(plotlist = pl, ncol = 2)
}



plot_png_expression_genes <- function(seu, genes, file_name, width = 8) {
  plot_expression(seu, genes, spatial_point_size = 0.4, umap_point_size = 0.4, grid_label_size = 10) |> 
    gs(file_name, width = width, height = 4 * length(genes))
}


plot_key_markers_umap_spatial <- function(seu, key_markers, file_prefix, width = 8) {
  all_feat <- Features(seu)

  km <- key_markers |> 
    filter(gene_symbol %in% all_feat)
  nums <- key_markers |> 
    count(cell_type)

  map(nums$cell_type, function(ct) {
    gname <- str_glue("{file_prefix}_{ct}")
    gns <- km |> 
      filter(cell_type == ct) |> 
      pull(gene_symbol) |> 
      as.character()
    plot_png_expression_genes(seu, gns, gname, width)
  }) 
}

plot_top_cluster_umap_spatial <- function(seu, top_cluster_markers, file_prefix = "top_clust_gene") {
  
  clusters <- top_cluster_markers$cluster |> 
    unique() |> 
    sort()
  
  map(clusters, function(cl) {
    gname <- sprintf("%s_cl%02d", file_prefix, cl)
    gns <- top_cluster_markers |> 
      filter(cluster == cl) |> 
      pull(gene_symbol)
    plot_png_expression_genes(seu, gns, gname)
  })
}


plot_cor_heat <- function(tbl_cor, fdr_limit = 0.01, min_cor = NULL, n_top = NULL, cluster_rows = TRUE, cluster_columns = TRUE) {
  env <- new.env(parent = globalenv())

  top_genes <- find_top_cor_genes(tbl_cor, fdr_limit, min_cor, n_top)

  # Make matrix symmetric for the heatmap plot
  tb <- symmetrise_table(tbl_cor)

  if(length(top_genes) == 0)
    return(NULL)
  tab <- tb |> 
    filter(gene_2 %in% top_genes & p_adj < fdr_limit) |> 
    pivot_wider(id_cols = gene_2, names_from = gene_1, values_from = correlation, values_fill = NA) |> 
    column_to_rownames("gene_2") |> 
    as.matrix()
  
  env$tab <- tab[, intersect(levels(tbl_cor$gene_1), colnames(tab))]
  env$palet <- circlize::colorRamp2(c(-1, 0, 1), c("#2166AC", "#FFFFFF", "#B2182B"), space = "LAB")
  env$cluster_columns <- cluster_columns
  env$cluster_rows = cluster_rows

  with(env, {
    ComplexHeatmap::Heatmap(tab, col = palet, name = "Correlation", cluster_columns = cluster_columns, cluster_rows = cluster_rows)
  })
  #ggheatmap(tab, palette = "distiller", with.y.text = TRUE, legend.name = "Correlation")
}

plot_gene_vs_gene <- function(seu, genes1, genes2, group_var = "cluster", group_sel = NULL, point_size = 0.5) {

  mt <- get_cell_data(seu) |> 
    mutate(group = get(group_var)) |> 
    select(cell_id, group)

  palette <- okabe_ito_palette
  if(!is.null(group_sel)) {
    mt <- mt |> 
      mutate(
        group = as.character(group),
        group = if_else(group %in% group_sel, group, "Other"),
        group = fct_relevel(group, "Other")
      ) 
    palette <- c("grey90", okabe_ito_palette)
  }

  d1 <- get_gene_expression(seu, genes1, norm_per_gene = FALSE) |> 
    summarise(gr1 = mean(expression), .by = cell_id)
  d2 <- get_gene_expression(seu, genes2, norm_per_gene = FALSE) |> 
    summarise(gr2 = mean(expression), .by = cell_id)

  inner_join(d1, d2, by = "cell_id") |> 
    left_join(mt, by = "cell_id") |> 
  ggplot(aes(x = gr1, y = gr2, colour = group)) +
    th +
    geom_point(size = point_size) +
    scale_colour_manual(values = palette) +
    labs(x = str_c(genes1, collapse = ","), y = str_c(genes2, collapse = ","), colour = group_var)
}


plot_cell_marker_pca <- function(seu, select_markers) {

  env <- new.env(parent = globalenv())

  tab <- get_group_gene_expression(seu, select_markers, group_var = "cell_type") |> 
    select(-cell_id) |> 
    as.matrix() |> 
    t()

  tab <- tab[apply(tab, 1, function(v) sum(v) > 0), ]
  tab <- tab[apply(tab, 1, sd) > 0, ]

  pca <- prcomp(t(tab), scale. = TRUE, center = TRUE)
  var.perc <- 100 * (pca$sdev)^2 / sum((pca$sdev)^2)
  
  env$pca1 <- sprintf("PCA1 (%5.1f%%)", var.perc[1])
  env$pca2 <- sprintf("PCA2 (%5.1f%%)", var.perc[2])

  env$xy <- pca$x |>
    as_tibble(rownames = "id") |>
    select(x = PC1, y = PC2, id)

  with(env, {
    xy |> 
      ggplot(aes(x = x, y = y)) +
      th +
      geom_point() +
      labs(x = pca1, y = pca2)
  })
}


make_list_for_upset <- function(mex) {
  d <- mex |> 
    pivot_longer(-cell_id) |> 
    filter(value > 0.5)
  groups <- unique(d$name)

  map(groups, function(gr) {
    d |> 
      filter(name == gr) |> 
      pull(cell_id)
  }) |> 
    set_names(groups)
}

make_combmat_for_upset <- function(mex, ...) {
  mex |> 
    make_list_for_upset() |> 
    ComplexHeatmap::make_comb_mat(...)
}


plot_cluster_cell_type_jaccard <- function(seu) {
  env <- new.env(parent = globalenv())

  env$d <- seu |>
    jaccard_clusters_cell_types()
  with(env, {
    d |> 
      ggplot(aes(x = cluster, y = cell_type, fill = j)) +
      th +
      geom_tile() +
      scale_fill_viridis_c(option = "cividis") +
      labs(fill = "Jaccard")
  })
}


plot_expression_pseudotime <- function(seu, gids, pt_column = "pseudotime_beta", point_size = 0.8) {
  ex <- get_gene_expression(seu, gids) |> 
    filter(expression > 0)
  pt <- seu@meta.data |> 
    as_tibble(rownames = "cell_id") |> 
    select(cell_id, pseudo_time = all_of(pt_column)) |> 
    drop_na()
  inner_join(ex, pt, by = "cell_id") |> 
    ggplot(aes(x = pseudo_time, y = expression)) +
    th +
    geom_point(size = point_size) +
    facet_wrap(~ gene_symbol, scales = "free_y") +
    labs(x = "Pseudo-time", y = "Normalised expression")
}




plot_volma <- function(res, group_var, point_size, point_alpha) {
  r <- res |>
    mutate(
      group = get(group_var)
    ) |>
    select(x, y, sel, group)
  r_sel <- r |> filter(sel)
  r_nsel <- r |> filter(!sel)

  rm(res, r)  # Minimise environment for serialisation

  g <- ggplot() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_point(data = r_nsel, aes(x = x, y = y), colour = "grey70",
               size = point_size, alpha = point_alpha) +
    facet_wrap( ~ group)

  if (nrow(r_sel) > 0) {
    g <- g + geom_point(data = r_sel, aes(x = x, y = y), colour = "black",
                        size = point_size, alpha = point_alpha)
  }
  g
}

plot_ma <- function(res, a = "AveExpr", fc = "logFC", p = "PValue", fdr = "FDR", group_var = "contrast",
                    fdr_limit = 0.05, point_size = 0.5, point_alpha = 0.5) {
  res |>
    mutate(
      x = get(a),
      y = get(fc),
      sel = get(fdr) < fdr_limit
    ) |>
    plot_volma(group_var, point_size, point_alpha) +
    geom_hline(yintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(log[10]~Intensity), y = expression(log[2]~FC))
}

plot_volcano <- function(res, fc = "avg_log2FC", p = "p_val", fdr = "p_val_adj", group_var = "contrast",
                         fdr_limit = 0.05, point_size = 0.5, point_alpha = 0.5) {
  res |>
    mutate(
      pval = get(p),
      pval = if_else(pval == 0, .Machine$double.xmin, pval),
      x = get(fc),
      y = -log10(pval),
      sel = get(fdr) < fdr_limit
    ) |>
    plot_volma(group_var, point_size, point_alpha) +
    geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
}

plot_pdist <- function(res, p = "PValue", group_var = "contrast", n_bins = 50) {
  brks <- seq(0, 1, length.out = n_bins)
  r <- res |>
    mutate(p = get(p), grp = get(group_var)) |>
    select(p, grp)
  rm(res)
  ggplot(r, aes(x = p, y = after_stat(density))) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_histogram(breaks = brks) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    facet_wrap(~grp)
}

plot_vmp <- function(res, a = "AveExpr", fc = "logFC", p = "PValue", fdr = "FDR", group = "contrast",
                     fdr_limit = 0.05, point_size = 0.5, point_alpha = 0.5, n_bins = 50){
  list(
    volcano = plot_volcano(res, fc, p, fdr, group, fdr_limit, point_size, point_alpha),
    ma = plot_ma(res, a, fc, p, fdr, group, fdr_limit, point_size, point_alpha),
    pdist = plot_pdist(res, p, group, n_bins)
  )
}

plot_cell_type_scores <- function(seu, tissue_markers, point_size = 0.4) {
  cell_types <- tissue_markers$cell_type |> unique()
  pls <- map(cell_types, function(ct) {
    ct <- str_glue("ES_{ct}")
    plot_spatial_meta(seu, ct, point_size)
  })
  plot_grid(plotlist = pls, ncol = 2)
}


plot_gene_overlap <- function(seu, gene_1, gene_2, section = NULL, point_size = 0.6) {
  only_1 <- str_glue("{gene_1} only")
  only_2 <- str_glue("{gene_2} only")

  xy <- get_spatial_coordinates(seu)
  ex <- get_gene_expression(seu, c(gene_1, gene_2), layer = "count") |> 
    pivot_wider(id_cols = cell_id, names_from = gene_symbol, values_from = expression) |> 
    rename(g1 := !!gene_1, g2 := !!gene_2) |> 
    mutate(
      state = case_when(
        g1 > 0 & g2 > 0 ~ "Both",
        g1 == 0 & g2 > 0 ~ only_2,
        g1 > 0 & g2 == 0 ~ only_1,
        g1 == 0 & g2 == 0 ~ "None"
      ) |> 
        factor(levels = c(only_1, only_2, "Both", "None"))
    )
  
  pl <- inner_join(xy, ex, by = "cell_id") |> 
    ggplot(aes(x = x, y = y, colour = state)) +
    thm +
    geom_point(size = point_size, shape = 15, stroke = 0.1) +
    scale_colour_manual(values = c(okabe_ito_palette[c(2, 3, 1)], "grey90")) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(x = NULL, y = NULL, colour = "Gene expression")

  if(!is.null(section)) {
    pl <- pl + 
      geom_path(data = section, aes(x = x, y = y), colour = "black", linewidth = 1)
  }

  pl
}


dist_to_segment <- function(xy, segment) {
  x1 <- segment$x[1]
  x2 <- segment$x[2]
  y1 <- segment$y[1]
  y2 <- segment$y[2]
  dx <- x2 - x1
  dy <- y2 - y1

  xy |> 
    mutate(
      # projection parameter t
      t = ((x - x1) * dx + (y - y1) * dy) / (dx^2 + dy^2),
      # projection point
      proj_x = x1 + t * dx,
      proj_y = y1 + t * dy,
      # Euclidean distance
      dst = sqrt((x - proj_x)^2 + (y - proj_y)^2)
    )
}



plot_expression_section <- function(seu, gene_1, gene_2, section, layer = "data", slice = "slice1.008um") {
  
  rad <- seu@images[[slice]]@boundaries$centroids@radius
  ex <- get_gene_expression(seu, c(gene_1, gene_2), layer = layer) |> 
    mutate(gene_symbol = factor(gene_symbol, levels = c(gene_1, gene_2)))
  ylab <- ifelse(layer == "count", "UMI count", "log1p normalised count")

  pl <- seu |> 
    get_spatial_coordinates() |> 
    dist_to_segment(section) |> 
    filter(dst < rad / 2 & t >=0 & t <= 1) |> 
    arrange(t) |> 
    left_join(ex, by = "cell_id") |> 
  ggplot(aes(x = t, y = expression, colour = gene_symbol)) +
    th +
    geom_line(, linewidth = 1.3) +
    scale_colour_manual(values = okabe_ito_palette[c(2, 3)]) +
    labs(x = "t", y = ylab, colour = "Gene")

  if(layer == "count")
    pl <- pl + geom_hline(yintercept = 0, colour = "brown", alpha = 0.3)
  pl
}

plot_spatial_section <- function(seu, gene_1, gene_2, section, point_size = 0.8, layer = "count", slice = "slice1.008um") {
  pl1 <- plot_gene_overlap(seu, gene_1, gene_2, section, point_size)
  pl2 <- plot_expression_section(seu, gene_1, gene_2, section, layer = layer, slice = slice)

  plot_grid(pl1, pl2, ncol = 1, rel_heights = c(2, 1), align = "v")
}

plot_spatial_difference <- function(seu, type_1, type_2, point_size = 0.6) {
  newcol <- str_glue("{type_1} - {type_2}")
  seu[[newcol]] <- seu[[type_1]] - seu[[type_2]]
  plot_spatial_meta(seu, meta_column = newcol, point_size = point_size, palette = "divergent") +
    theme(
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7),
      legend.key.size = unit(0.4, 'cm')
    )
}

# pairs is a list of vectors
plot_spatial_differences <- function(seu, pairs, point_size = 0.6) {
  pl <- map(pairs, function(pair) {
    ct1 <- paste0("ES_", pair[1])
    ct2 <- paste0("ES_", pair[2])
    plot_spatial_difference(seu, ct1, ct2, point_size)
  })

  plot_grid(plotlist = pl, ncol = 1)
}


plot_upset_expressed_cells <- function(seu, gns, expr_limit = 0) {
  m <- get_gene_expression(seu, gns, layer = "data") |> 
    filter(expression > expr_limit) |> 
    mutate(expression = 1) |> 
    pivot_wider(id_cols = cell_id, names_from = gene_symbol, values_from = expression, values_fill = 0) |> 
    ComplexHeatmap::make_comb_mat()
  
  ComplexHeatmap::UpSet(
    m,
    top_annotation = ComplexHeatmap::upset_top_annotation(m, add_numbers = TRUE, numbers_rot = 0),
    right_annotation = ComplexHeatmap::upset_right_annotation(m, add_numbers = TRUE),
    comb_col = c("black", "brown", "darkblue", "darkgreen", "darkorchid3")[ComplexHeatmap::comb_degree(m)]
  )
}

combat_upset_expressed_cells <- function(seu, gns, suffix, sample, expr_limit = 0) {
  rfile <- file.path("doc", "tables", str_glue("{suffix}_{sample}.rds"))
  m <- get_gene_expression(seu, gns, layer = "data") |> 
    filter(expression > expr_limit) |> 
    mutate(expression = 1) |> 
    pivot_wider(id_cols = cell_id, names_from = gene_symbol, values_from = expression, values_fill = 0) |> 
    ComplexHeatmap::make_comb_mat()
  readr::write_rds(m, rfile)
}


png_upset_expressed_cells <- function(seu, gns, suffix, sample, width, height, expr_limit = 0) {
  pfile <- file.path("doc", "images", str_glue("{suffix}_{sample}.png"))

  png(pfile, width = width, height = height)

  plot_upset_expressed_cells(seu, gns, expr_limit)

  dev.off()
}

plot_spatial_box <- function(seu, gene, x_cen, y_cen, size, point_size = 0.6) {
  xy <- get_spatial_coordinates(seu)
  ex <- get_gene_expression(seu, gene)

  bx <- tibble(
    x = c(x_cen - size, x_cen + size, x_cen + size, x_cen - size, x_cen - size),
    y = c(y_cen - size, y_cen - size, y_cen + size, y_cen + size, y_cen - size)
  )

  inner_join(xy, ex, by = "cell_id") |> 
    ggplot() +
    th +
    theme(legend.position = "none") +
    geom_point(aes(x = x, y = y, colour = expression), size = point_size) +
    scale_colour_distiller(palette = "Blues", direction = 1) +
    geom_path(data = bx, aes(x = x, y = y))
}

plot_resolution_comparison <- function(res_comp) {
  rg <- res_comp |> 
    filter(resolution == "8um") |> 
    get_resolution_grid()

  # Add fifth corner for completion
  rg5 <- rg |> 
    filter(corner_idx == 1) |> 
    mutate(corner_idx = 5)
  rg <- bind_rows(rg, rg5)

  rg <- rg |> 
    arrange(centre_id, corner_idx)

  res_comp |> 
    ggplot() +
    th +
    theme(legend.position = "none") +
    geom_text(aes(x = x, y = y, colour = resolution, label = expression)) +
    geom_path(data = rg, aes(x = x_corner, y = y_corner, group = centre_id)) +
    facet_wrap(~gene_symbol) +
    scale_colour_manual(values = c("grey70", "black"))
}



get_resolution_grid <- function(res_comp) {
  df <- res_comp |> 
    mutate(centre_id = cell_id |> as.factor() |> as.integer()) |> 
    distinct(cell_id, centre_id, x, y)

  # This part was written by ChatGPT o4-mini-high and improved by MG:

  knn_res <- FNN::get.knn(df[, c("x","y")], k = 2)
  nearest_idx <- knn_res$nn.index[1, 2]   # the second‐closest to row 1
  nearest_dist <- knn_res$nn.dist[1, 2]   # the distance to that neighbour

  # Extract the vector from centre[1,] to its nearest neighbour:
  v_x <- df$x[nearest_idx] - df$x[1]
  v_y <- df$y[nearest_idx] - df$y[1]

  # (3) Side‐length = s = distance between two adjacent centres:
  s <- sqrt(v_x^2 + v_y^2)

  # Form a unit‐vector along that grid direction:
  u1 <- c(v_x, v_y) / s    # length‐1
  # A perpendicular unit‐vector (rotate u1 by +90° CCW):
  u2 <- c(-u1[2], u1[1])   # also length‐1

  # (4) Corner‐offsets = ±(s/2)*u1 ± (s/2)*u2.  Define the two half‐offsets:
  half1 <- (s/2) * u1      # half along axis‐1
  half2 <- (s/2) * u2      # half along axis‐2

  # Precompute the four combinations of ( ±half1  ±half2 ):
  corner_offsets <- rbind(
    half1 + half2,
    half1 - half2,
    -half1 - half2,
    -half1 + half2
  )
  # corner_offsets is a 4×2 matrix; row 1 = +half1 +half2, etc.

  map(1:4, function(k) {
    df |> 
      mutate(
        corner_idx = k,
        x_corner = x +  corner_offsets[k, 1],
        y_corner = y + corner_offsets[k, 2]
      )
  }) |> 
    list_rbind() |> 
    group_by(centre_id) |>     
    arrange(cell_id, corner_idx) |> 
    ungroup()
}


plot_cluster_cell_type_overlap <- function(seu, ct = "Unknown") {
  ov <- cluster_cell_type_overlap(seu)

  ov |> 
    filter(cell_type == ct) |> 
    ggplot(aes(x = cluster, y = prop)) +
    th +
    geom_col() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = "Cluster", y = "Proporion of spots in cluster", title = ct)
}

plot_binned_xy_expression <- function(seu, gene_x, genes_y, bins = 6, exclude_zeroes = FALSE) {
  gx <- get_gene_expression(seu, gene_x)
  g1 <- gx |> 
    filter(expression > 0) |> 
    mutate(bin = cut_interval(expression, n = bins) |> fct_expand("0", after = 0))
  g0 <- gx |> 
    filter(expression == 0) |> 
    mutate(bin = factor("0", levels = levels(g1$bin)))
  if(exclude_zeroes) {
    gb <- g1
  } else {
    gb <- bind_rows(g0, g1)
  }

  gb <- gb |> 
    select(cell_id, bin)

  get_gene_expression(seu, genes_y) |> 
    mutate(gene_symbol = factor(gene_symbol, levels = genes_y)) |> 
    inner_join(gb, by = "cell_id")  |> 
    group_by(gene_symbol, bin) |> 
    summarise(
      m = mean(expression),
      se = sd(expression) / sqrt(n()),
      n = n()
    ) |> 
    ggplot(aes(x = bin, y = m, ymin = m - se, ymax = m + se)) +
    th +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
    geom_errorbar(width = 0, colour = "grey70") +
    geom_point() +
    facet_wrap(~gene_symbol) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    labs(x = str_glue("{gene_x} expression bin"), y = "Mean expression")
}