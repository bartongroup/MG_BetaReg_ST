mn_plot_spatial_clustering <- function(seu, res) {
  plot_one_xy(seu, resolution = res, what = "spatial") +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      legend.key.size = unit(0.5, "lines")
    ) +
    labs(fill = "Cluster")
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


mn_plot_hormone_overlap <- function(seu, fname, width, height) {
  up <- plot_upset_expressed_cells(seu, c("INS", "SST", "GCG", "GHRL", "PPY")) 

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

