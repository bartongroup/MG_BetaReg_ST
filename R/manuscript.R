mn_plot_spatial_clustering <- function(seu, res = NULL, meta_column = "cluster", swap_xy = FALSE) {
  plot_one_xy(seu, resolution = res, meta_column = meta_column, what = "spatial", swap_xy = swap_xy,
             discrete_palette = modified_brewer_paired) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      legend.key.size = unit(0.5, "lines")
    ) +
    labs(fill = "Cluster")
}


mn_plot_cluster_cell_type_jaccard <- function(seu, cluster_col = "cluster", scale_limits = c(0, NA)) {
  env <- new.env(parent = globalenv())
  env$scale_limits <- scale_limits

  env$d <- seu |>
    jaccard_clusters_cell_types(cluster_col = cluster_col)
  with(env, {
    d |> 
      ggplot(aes(x = cluster, y = cell_type, fill = j)) +
      th +
      geom_tile() +
      scale_fill_viridis_c(option = "cividis", limits = scale_limits) +
      labs(fill = "Jaccard", x = "Cluster", y = "Cell type")
  })
}


mn_plot_marker_spatial <- function(seu, markers, label_size = 2.5, ncol = NULL) {
  plts <- map2(markers$gene_symbol, markers$cell_type, function(gene, ct) {
    pl <- plot_one_xy(seu, genes = gene, what = "spatial")

    # Right-aligning a label is a pain. Need to do this in original x-y coordinates.
    xmax <- max(pl$data$x)
    ymax <- max(pl$data$y)

    pl +
      theme(plot.margin = margin(2, 2, 2, 2)) +
      annotate("text", x = xmax, y = ymax, label = gene, hjust = 1, vjust = 1, size = label_size) +
      annotate("text", x = xmax, y = ymax, label = ct, hjust = 1, vjust = 3, size = label_size * 0.8, colour = "grey60")
  }) 
  plot_grid(plotlist = plts, ncol = ncol)
}



mn_plot_cell_type <- function(seu, label_size = 3) {
  pal <- c(okabe_ito_palette[c(4, 1, 8, 5, 2, 3, 6, 7, 9)], "darkgreen", "lightgreen", "red", "brown", "magenta", grey(0.95))
  plot_one_xy(seu, label_size = label_size, what = "spatial", meta_column = "cell_type", discrete_palette = pal, with_legend = TRUE) + guides(colour = guide_legend(override.aes = list(size = 3))) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      legend.key.size = unit(0.5, "lines")
    ) +
    labs(fill = "Cell type")
}


mn_plot_cell_type_dots <- function(seu_mark, tissue_markers) {
  plot_gene_dots(seu_mark, unique(tissue_markers$gene_symbol), meta_column = "cell_type") +
    labs(x = "Cell type", y = "Marker gene")
}


mn_plot_hormone_overlap <- function(seu, fname, width, height, genes = c("INS", "SST", "GCG", "GHRL", "PPY", "NEUROG3")) {
  up <- plot_upset_expressed_cells(seu, genes) 

  pdf(fname, width = width, height = height)
  on.exit(dev.off(), add = TRUE)  
  ComplexHeatmap::draw(up)
}


# Needs object calculated on the cluster and retrieved
mn_plot_hormone_overlap_2um <- function(smpl, fname, width, height) {
  m <- readr::read_rds(paste0("doc/tables/hormone_overlap_2um_", smpl, ".rds"))
  up <- ComplexHeatmap::UpSet(
    m,
    top_annotation = ComplexHeatmap::upset_top_annotation(m, add_numbers = TRUE, numbers_rot = 0),
    right_annotation = ComplexHeatmap::upset_right_annotation(m, add_numbers = TRUE),
    comb_col = c("black", "brown", "darkblue", "darkgreen", "darkorchid3")[ComplexHeatmap::comb_degree(m)]
  )

  pdf(fname, width = width, height = height)
  on.exit(dev.off(), add = TRUE)  
  ComplexHeatmap::draw(up)
}

# Need figures ready from another project. If absent, just ignore.
mn_plot_hormone_overlap_ma <- function(pcw, fname, width, height) {
  rds_file <- file.path("../BetaReg_scMark/fig", str_glue("hormone_upset_{pcw}pcw.rds"))
  if(!file.exists(rds_file))
    return(NULL)
  up <- read_rds(rds_file)
  
  pdf(fname, width = width, height = height)
  on.exit(dev.off(), add = TRUE)  
  ComplexHeatmap::draw(up)
}



mn_plot_binned_spatial <- function(seu, gene_x = "INS", bins = 6) {
  gx <- get_gene_expression(seu, gene_x)
  g1 <- gx |> 
    filter(expression > 0) |> 
    mutate(bin = cut_interval(expression, n = bins) |> fct_expand("0", after = 0))
  g0 <- gx |> 
    filter(expression == 0) |> 
    mutate(bin = factor("0", levels = levels(g1$bin)))
  gx <- bind_rows(g0, g1)
  # reorder before inserting into meta
  m <- tibble(cell_id = colnames(seu)) |> 
    left_join(gx, by = "cell_id")
  seu$bin <- m$bin

  plot_one_xy(seu, what = "spatial", meta_column = "bin") + 
    theme(
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      legend.key.size = unit(0.5, "lines")
    ) +
    labs(fill = "INS intensity bin")
}


mn_plot_binned_enrichment <- function(tbl) {
  tb <- tbl |>
    filter(bin != "0") |>
    droplevels() |> 
    mutate(
      ibin = as.integer(bin),
      infor = is.infinite(odds_ratio),
      logp = -log10(p_adjust),
      lor = log10(odds_ratio)
    )

  mx_scale <- 1.05
  mxo <- tb |> 
    filter(!infor) |> 
    pull(lor) |> 
    max()
  mxp <- tb |> 
    filter(!is.infinite(logp)) |> 
    pull(logp) |> 
    max()

  tb <- tb |> 
    mutate(
      lor = if_else(infor, mx_scale * mxo, lor),
      logp = if_else(is.infinite(logp), mxp, logp)
    )

  bin_counts <- tb |> 
    count(ibin) |> 
    pull(n)
  # Extend top and bottom to account for legend and x-scale
  bin_counts[1] <- bin_counts[1] + 5
  bin_counts[length(bin_counts)] <- bin_counts[length(bin_counts)] + 3

  plts <- map(1:max(tb$ibin), function(ibn) {
    pl <- tb |> 
      filter(ibin == ibn) |> 
      mutate(term_name = fct_reorder(term_name, enrichment)) |>
      ggplot(aes(x = lor, y = term_name)) +
      th +
      guides(colour = "none", linetype = "none") +
      geom_segment(aes(xend = 0, yend = term_name, linetype = infor), colour = "grey70") +
      geom_point(aes(size = logp, colour = infor)) +
      scale_size_continuous(limits = c(0, mxp)) +
      scale_colour_manual(values = c("black", "grey80")) +
      scale_linetype_manual(values = c("solid", "dashed")) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, mx_scale * mxo)) +
      labs(x = expression(log[10]~Odds~ratio), y = NULL, size = expression(-log[10]~FDR)) +
      facet_wrap(~ bin, strip.position = "right")
    if(ibn == 1) {
      pl <- pl + theme(legend.position = "top")
    } else {
      pl <- pl + theme(legend.position = "none")
    }
    if(ibn != max(tb$ibin)) {
      pl <- pl + theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
      )
    }
    pl

  })

  plot_grid(plotlist = plts, ncol = 1, rel_heights = bin_counts, align = "v")
}


mn_plot_spatial_selection_genes <- function(seu, genes, box = list(x1 = -Inf, x2 = Inf, y1 = -Inf, y2 = Inf),
                                            ncol = 2, n_grid = 10) {
  ex <- get_gene_expression(seu, genes) |> 
    mutate(val = expression)
  xy <- get_spatial_coordinates(seu, swap_xy = FALSE) |> 
    mutate(cell = cell_id |> str_remove("s_008um_") |> str_remove("-1$")) |> 
    separate_wider_delim(cell, delim = "_", names = c("x_grid", "y_grid")) |> 
    mutate(across(c(x_grid, y_grid), as.integer))
    
  # Select grid coordinates, not x-y
  grid_sel <-xy |> 
    filter(x >= box$x1 & x <= box$x2 & y >= box$y1 & y <= box$y2) |> 
    distinct(x_grid, y_grid)

  xy <- xy |> 
    filter(x_grid %in% grid_sel$x_grid & y_grid %in% grid_sel$y_grid)

  ay <- xy |> arrange(y_grid, x_grid)
  dx <- abs(ay[2, ]$x - ay[1, ]$x)
  ax <- xy |> arrange(x_grid, y_grid)
  dy <- abs(ax[2, ]$y - ax[1, ]$y)

  corners <- get_resolution_grid(xy)
  d <- inner_join(corners, ex, by = "cell_id", relationship = "many-to-many")
  
  # The 10x10 grid is going through the middle of pixels, need to shift by hafl width
  grid <- xy |> 
    filter(x_grid %% n_grid == 0 & y_grid %% n_grid == 0) |>     
    get_resolution_grid() |> 
    mutate(
      x_corner = x_corner + dx / 2,
      y_corner = y_corner + dy / 2
    )

  d |> 
    ggplot() +
    thm +
    geom_polygon(aes(x = x_corner, y = y_corner, group = cell_id, fill = val), colour = NA, linewidth = 0) +
    labs(x = NULL, y = NULL, fill = "Expression") +
    facet_wrap(~gene_symbol, ncol = ncol) +
    #scale_fill_viridis_c(option = "cividis") +
    scale_fill_distiller(palette = "Purples", direction = 1) +
    geom_polygon(data = grid, aes(x = x_corner, y = y_corner, group = cell_id), fill = NA, colour = "white", alpha = 0.5, linewidth = 0.3)
}


mn_plot_ins_gradient_expression <- function(seu) {
  plot_binned_expression(seu, "INS", c("NEUROG3", "GCG", "SST"), bins = 6)
}

mn_plot_wgcna_modules_expr_score <- function(seu, n_genes = "all", ncol = 2) {
  library(UCell)

  xy <- get_spatial_coordinates(seu, swap_xy = FALSE)
  corners <- get_resolution_grid(xy)
  
  ex <-  seu |> 
    hdWGCNA::ModuleExprScore(n_genes = n_genes, method = "UCell", ncores = 8) |>    # this is a bit slow
    hdWGCNA::GetModuleScores() |>
    as_tibble(rownames = "cell_id")

  d <- inner_join(corners, ex, by = "cell_id") |> 
    select(-c(centre_id, x, y, corner_idx)) |>
    pivot_longer(-c(cell_id, x_corner, y_corner), names_to = "module") 

  d |> 
    ggplot() +
    thm +
    geom_polygon(aes(x = x_corner, y = y_corner, group = cell_id, fill = value, colour = value), linewidth = 0.2) +
    facet_wrap(~ module, ncol = ncol) +
    labs(x = NULL, y = NULL, fill = "Expr score", colour = "Expr score") +
    scale_fill_distiller(palette = "Purples", direction = 1) +
    scale_colour_distiller(palette = "Purples", direction = 1)
}


mn_plot_cellchat_rank <- function(cellchat) {

  d2 <- CellChat::rankNet(cellchat, mode = "single", measure = "weight", return.data = TRUE)$signaling.contribution |> as_tibble()
  
  plot_one_rank(d2, "Information flow")
}
