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



mn_plot_cell_type_dots <- function(seu, tissue_markers) {
  plot_gene_dots(seu, unique(tissue_markers$gene_symbol), meta_column = "cell_type") +
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

