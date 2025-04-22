run_slingshot <- function(seu, cell_types = c("MPC", "EP", "Beta"), reduction = "PCA") {
  # Subset Seurat to only include cell types of interest and convert into SingleCellExperiment
  seu@meta.data <- seu@meta.data[, c("cell_type"), drop = FALSE]
  sce <- seu |> 
    subset(subset = cell_type %in% cell_types) |> 
    as.SingleCellExperiment()

  # Infer start and end from the cell_types argument
  start_clus <- head(cell_types, 1)
  end_clus <- tail(cell_types, 1)

  slingshot::slingshot(
    sce,
    clusterLabels = "cell_type",
    start.clus = start_clus,
    end.clus = end_clus,
    reducedDim = reduction
  )
}

add_pseudotime <- function(seu, sl, lineage = 1, colname = "pseudotime", ...)  {
  lng <- str_glue("slingPseudotime_{lineage}")
  
  seu@meta.data[colnames(sl), colname] <- sl[[lng]]
  seu
}

plot_slingshot_lineage <- function(sl, lineage = 1, point_size = 0.6) {
  lng <- str_glue("Lineage{lineage}")
  
  d <- slingshot::SlingshotDataSet(sl)
  
  curves <- d@curves[[lng]]$s |>
    as_tibble() |> 
    select(1, 2) |> 
    set_names(c("x", "y"))

  red <- d@reducedDim |> 
    as_tibble() |> 
    select(1, 2) |> 
    set_names(c("x", "y")) |> 
    mutate(cell_type = SummarizedExperiment::colData(sl)$cell_type)

  axes <- colnames(d@reducedDim)[1:2]

  ggplot() +
    thm +
    geom_point(data = red, aes(x = x, y = y, colour = cell_type), size = point_size) +
    scale_colour_brewer(palette = "Paired") +
    geom_path(data = curves, aes(x = x, y = y), colour = "brown", linewidth = 1.5) +
    labs(x = axes[1], y = axes[2])
}