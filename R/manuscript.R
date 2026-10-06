oi_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#AAAAAA", "#666666")
cell_type_pal <- c(oi_palette[c(4, 1, 8, 5, 2, 3, 6, 7, 9)], "darkgreen", "lightgreen", "red", "brown", "magenta", grey(0.95))


mn_plot_meta_segmentation <- function(seu, meta_col, box = list(x1 = NA, x2 = NA, y1 = NA, y2 = NA),
                                      with_centroids = FALSE, cell_sel = NULL, palette = NULL, fill_lab = NULL) {
  poly <- get_polygons(seu)
  meta <- seu |> 
    get_cell_data() |>
    mutate(col = get(meta_col)) |> 
    select(cell_id, col)

  if(with_centroids) {
    centroids <- get_spatial_coordinates(seu, swap_xy = FALSE)
  }

  d <- inner_join(poly, meta, by = "cell_id", relationship = "many-to-many") 
  if(!is.null(cell_sel)) {
    d <- d |> 
      filter(cell_id %in% cell_sel)
  }

  pl <- d |> 
    ggplot() +
    thm +
    geom_polygon(aes(x = x, y = y, group = cell_id, fill = col), colour = NA, linewidth = 0) +
    coord_fixed(xlim = c(box$x1, box$x2), ylim = c(box$y1, box$y2)) +
    labs(x = NULL, y = NULL, fill = ifelse(is.null(fill_lab), meta_col, fill_lab))
    
  if(with_centroids) 
    pl <- pl + geom_point(data = centroids, aes(x = x, y = y), colour = "black", size = 0.5)
    
  if(is.null(palette)) {
    pl <- pl + scale_fill_brewer(palette = "Paired")
  } else {
    pl <- pl + scale_fill_manual(values = palette)
  }

  pl
}


mn_plot_binned_segmentation <- function(seu, gene_x = "INS", bins = 6, box = list(x1 = NA, x2 = NA, y1 = NA, y2 = NA)) {
  gx <- gene_binned_gradient(seu, gene_x, bins)
  
  # reorder before inserting into meta
  m <- tibble(cell_id = colnames(seu)) |> 
    left_join(gx, by = "cell_id")
  seu$bin <- m$bin

  mn_plot_meta_segmentation(seu, meta_col = "bin", box = box) +
    labs(fill = "INS intensity bin")
}


mn_plot_gene_segmentation <- function(seu, genes, box = list(x1 = NA, x2 = NA, y1 = NA, y2 = NA),
                                      seu_grid = NULL, add_grid = FALSE, with_centroids = FALSE, ncol = 2) {
  poly <- get_polygons(seu)
  ex <- get_gene_expression(seu, genes)

  if(with_centroids) {
    centroids <- get_spatial_coordinates(seu, swap_xy = FALSE)
  }

  if(!is.null(seu_grid)) {
    xy <- get_spatial_coordinates(seu_grid, swap_xy = FALSE)
    corners <- get_resolution_grid(xy)
  }

  d <- inner_join(poly, ex, by = "cell_id", relationship = "many-to-many") 
  pl <- d |> 
    ggplot() +
    thm +
    geom_polygon(aes(x = x, y = y, group = cell_id, fill = expression), colour = NA, linewidth = 0) +
    scale_fill_distiller(palette = "Purples", direction = 1) +
    #scale_fill_gradient(low = "grey90", high = "black") +
    facet_wrap(~gene_symbol, ncol = ncol) +
    coord_fixed(xlim = c(box$x1, box$x2), ylim = c(box$y1, box$y2)) +
    labs(x = NULL, y = NULL)  
  if(!is.null(seu_grid)) {
    pl <- pl +
      geom_path(data = corners, aes(x = x_corner, y = y_corner, group = cell_id), colour = "black", linewidth = 0.3)
  }
  if(with_centroids) {
    pl <- pl +
      geom_point(data = centroids, aes(x = x, y = y), colour = "black", size = 0.5)
  }
  if(add_grid) {
    pl <- pl +
      theme(
        panel.grid.major = element_line(colour = "brown", linewidth = 0.2),
        panel.background = element_rect(fill = NA),
        panel.ontop = TRUE
      )
  }
  pl +
    labs(fill = "Expression")
}





mn_plot_gene_dots <- function(seu, gns, meta_column = "cluster", with_legend = TRUE, title = NULL) {
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
      #scale_colour_gradient(low = "grey90", high = "blue") +
      scale_colour_distiller(palette = "Purples", direction = 1) +
      scale_radius(range = c(1, 6)) +
      labs(x = meta_column, y = "Gene", size = "Percent\nexpressed", colour = "Average\nexpresion", title = title)
    if(!with_legend)
      pl <- pl + theme(legend.position = "none")
    pl
  })
}


mn_plot_cell_type_dots <- function(seu, tissue_markers) {
  mn_plot_gene_dots(seu, unique(tissue_markers$gene_symbol), meta_column = "cell_type") +
    labs(x = "Cell state", y = "Marker gene")
}



mn_plot_cluster_cell_type_jaccard <- function(seu, cluster_col = "cluster", scale_limits = c(0, NA)) {
  env <- new.env(parent = globalenv())
  env$scale_limits <- scale_limits

  env$d <- seu |>
    jaccard_clusters_cell_types(cluster_col = cluster_col) |> 
    mutate(x = get(cluster_col))
  with(env, {
    d |> 
      ggplot(aes(x = x, y = cell_type, fill = j)) +
      th +
      geom_tile() +
      scale_fill_viridis_c(option = "cividis", limits = scale_limits) +
      labs(fill = "Jaccard", x = "Cluster", y = "Cell type")
  })
}

mn_plot_hormone_overlap <- function(seu, fname, width, height, genes = c("INS", "SST", "GCG", "GHRL", "PPY", "NEUROG3")) {
  up <- plot_upset_expressed_cells(seu, genes) 

  pdf(fname, width = width, height = height)
  on.exit(dev.off(), add = TRUE)  
  ComplexHeatmap::draw(up)
}


mn_plot_binned_expression <- function(seu, gene_x, genes_y, bins = 6) {
  b <- gene_binned_gradient(seu, gene_x, bins)
  seu$bin <- b$bin

  plot_gene_expression(seu, genes_y, meta_column = "bin") + ylab("Expression")
}


mn_plot_enrichment_gradient <- function(enr) {
  mx <- enr |> 
    filter(!is.infinite(odds_ratio)) |> 
    pull(odds_ratio) |> 
    max()
  d <- enr |> 
    filter(bin != 0) |>
    mutate(
      odds_ratio = if_else(is.infinite(odds_ratio), mx * 1.1, odds_ratio) |> log(10),
      ibin = as.integer(bin)
    ) |>
    group_by(term_name) |> 
    mutate(cm = sum(ibin * odds_ratio) /  sum(odds_ratio)) |> 
    ungroup() |> 
    mutate(term_name = fct_reorder(term_name, cm)) |> 
    pivot_wider(id_cols = term_name, names_from = bin, values_from = odds_ratio, values_fill = 0) |>
    arrange(term_name) |> 
    column_to_rownames("term_name")

  ggheatmap(d, order.col = FALSE, order.row = FALSE, with.y.text = TRUE, palette = "viridis", legend.name = expression(log[10]~OR))
}


mn_plot_cellchat_rank <- function(cellchat, n_top = Inf) {
  d2 <- CellChat::rankNet(cellchat, mode = "single", measure = "weight", return.data = TRUE)$signaling.contribution |>
    as_tibble() |> 
    slice_max(contribution, n = n_top)
  
  plot_one_rank(d2, "Information flow")
}

mn_plot_cellchat_chord <- function(cellchat, pathway, pfile) {

  pdf(file = file.path("manuscript", "fig", pfile), width = 8, height = 8)
  par(mfrow = c(1, 1))

  CellChat::netVisual_chord_cell(
    cellchat,
    signaling = pathway,
    lab.cex = 1.5
  ) |> 
    print()
  
  dev.off()
}

mn_plot_cellchat_chords <- function(cellchat, sample, pathways, prefix = "cellchat_chord", type = "png") {
  for(path in pathways) {
    pfile = str_glue("{prefix}_{sample}_{path}.{type}")
    mn_plot_cellchat_chord(cellchat, path, pfile)
  }
}


plot_one_simple_segmentation <- function(seu, gene, box) {
  poly <- get_polygons(seu)
  ex <- get_gene_expression(seu, gene)
  d <- inner_join(poly, ex, by = "cell_id", relationship = "many-to-many") 
  d |> 
    ggplot() +
    thm +
    theme(legend.position = "none") +
    geom_polygon(aes(x = x, y = y, group = cell_id, fill = expression), colour = NA, linewidth = 0) +
    scale_fill_distiller(palette = "Purples", direction = 1) +
    coord_fixed(xlim = c(box$x1, box$x2), ylim = c(box$y1, box$y2)) +
    labs(x = NULL, y = NULL)  
}


mn_plot_marker_segmentation <- function(seu, markers, box, label_size = 2.5, ncol = NULL) {
  plts <- map2(markers$gene_symbol, markers$cell_type, function(gene, ct) {
    pl <- plot_one_simple_segmentation(seu, gene, box)

    xmax <- box$x2
    ymax <- box$y2

    pl +
      theme(plot.margin = margin(2, 2, 2, 2)) +
      annotate("text", x = xmax, y = ymax, label = gene, hjust = 1, vjust = 0.3, size = label_size) +
      annotate("text", x = xmax, y = ymax, label = ct, hjust = 1, vjust = 2, size = label_size * 0.8, colour = "grey60")
  }) 
  plot_grid(plotlist = plts, ncol = ncol)
}



mn_plot_modules_expr_score <- function(seu, box, n_genes = "all", ncol = 2, sel_modules = NULL, ncores = 1) {
  library(UCell)

  poly <- get_polygons(seu)
  ex <-  seu |> 
    hdWGCNA::ModuleExprScore(n_genes = n_genes, method = "UCell", ncores = ncores) |>    # this is a bit slow
    hdWGCNA::GetModuleScores() |>
    as_tibble(rownames = "cell_id")

  d <- inner_join(poly, ex, by = "cell_id", relationship = "many-to-many") |> 
    pivot_longer(cols = -c(cell_id, x, y), names_to = "module", values_to = "score")
  if (!is.null(sel_modules)) {
    d <- d  |> filter(module %in% sel_modules)
  }
  d |> 
    ggplot() +
    thm +
    geom_polygon(aes(x = x, y = y, group = cell_id, fill = score), colour = NA, linewidth = 0) +
    scale_fill_distiller(palette = "Purples", direction = 1) +
    facet_wrap(~module, ncol = ncol) +
    coord_fixed(xlim = c(box$x1, box$x2), ylim = c(box$y1, box$y2)) +
    labs(x = NULL, y = NULL, fill = "Expr. score")  
}
