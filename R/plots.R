okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#AAAAAA", "#000000")

th <- ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank()
  )


gs <- function(gg, exprm, name, width, height) {
  ggsave(filename = file.path("fig", paste0(exprm, "_", name, ".png")), plot = gg, device = "png",
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

plot_all_counts <- function(cell_data_all) {
  pl1 <- cell_data_all |> 
    count(sample) |> 
    ggplot(aes(x = sample, y = n)) +
    th +
    geom_col() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Sample", y = "Cell count")
  pl2 <- cell_data_all |> 
    summarise(n = sum(total_count) / 1e6, .by = "sample") |> 
    ggplot(aes(x = sample, y = n)) +
    th +
    geom_col() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Sample", y = "Read count (millions)")
  plot_grid(pl1, pl2, nrow = 1)
}

plot_seurat_clusters <- function(seu, resolutions = NULL, with_legend = FALSE, label = TRUE, ...) {
  if(is.null(resolutions)) {
    resols <- names(seu@meta.data) |> str_subset("snn_res")
  } else {
    resols <-  get_resolution_columns(resolutions)
  }
  
  P <- map(resols, function(res) {
    Idents(seu) <- res
    r <- str_remove(res, "^.*res\\.")
    dp <- DimPlot(seu, reduction = "umap", label = label, label.size = 6, ...) +
      labs(title = r)
    if(!with_legend)
      dp <- dp + theme(legend.position = "none")
    dp
  })
  plot_grid(plotlist = P)
}

