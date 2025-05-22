run_slingshot <- function(seu, meta_col = "cell_type", clusters = c("MPC", "EP", "Beta"), reduction = "PCA", ...) {
  # Subset Seurat to only include cell types of interest and convert into SingleCellExperiment
  seu@meta.data$cl  <- seu@meta.data[[meta_col]]
  seu@meta.data <- seu@meta.data[, "cl", drop = FALSE]
  
  sce <- seu |> 
    subset(subset = cl %in% clusters) |> 
    as.SingleCellExperiment()

  # Infer start and end from the clusters argument
  start_clus <- head(clusters, 1)
  end_clus <- tail(clusters, 1)

  slingshot::slingshot(
    sce,
    clusterLabels = "cl",
    start.clus = start_clus,
    end.clus = end_clus,
    reducedDim = reduction,
    ...
  )
}

add_pseudotime <- function(seu, sl, lineage = 1, colname = "pseudotime", ...)  {
  lng <- str_glue("slingPseudotime_{lineage}")
  
  seu@meta.data[colnames(sl), colname] <- sl[[lng]]
  seu
}

plot_slingshot_lineage <- function(sl, lineage = 1, point_size = 0.6) {
  env <- new.env(parent = globalenv())
  env$point_size <- point_size
  
  lng <- str_glue("Lineage{lineage}")
  
  d <- slingshot::SlingshotDataSet(sl)
  
  env$curves <- d@curves[[lng]]$s |>
    as_tibble() |> 
    select(1, 2) |> 
    set_names(c("x", "y"))

  env$red <- d@reducedDim |> 
    as_tibble() |> 
    select(1, 2) |> 
    set_names(c("x", "y")) |> 
    mutate(cl = SummarizedExperiment::colData(sl)$cl)

  env$axes <- colnames(d@reducedDim)[1:2]

  with(env, {
    ggplot() +
      thm +
      geom_point(data = red, aes(x = x, y = y, colour = cl), size = point_size) +
      scale_colour_brewer(palette = "Paired") +
      geom_path(data = curves, aes(x = x, y = y), colour = "brown", linewidth = 1.5) +
      labs(x = axes[1], y = axes[2])
  })
}