create_cellchat <- function(seu, spot_diameter_fullres = 17.56, label_var = "cell_type", spot_size = 8,
                            sub_n_cells = 5000) {

  # seu_sub <- SketchData(seu, ncells = sub_n_cells)

  # Normalised expression data
  data_input <- GetAssayData(seu, layer = "data")
  # sub_cells <- colnames(data_input)  # sketching does not subset metadata and spatial data, need to do it manually

  # Metadata
  meta <- seu@meta.data |> 
    select(all_of(label_var)) |> 
    add_column(samples = "sample") |> 
    mutate(samples = as.factor(samples))
  # meta <- meta[sub_cells, ]

  # Spatial data
  spatial_factors <- tibble(
    ratio = spot_size / spot_diameter_fullres,
    tol = spot_size / 2
  )
  
  spatial_locs <- GetTissueCoordinates(seu) |> 
    select(x, y)
  # spatial_locs <- spatial_locs[sub_cells, ]


  CellChat::createCellChat(
    object = data_input,
    meta = meta,
    group.by = label_var,
    datatype = "spatial",
    coordinates = data.matrix(spatial_locs),
    spatial.factors = spatial_factors
  )
}

preprocess_cellchat <- function(cellchat, n_cores = 8) {
  use_db <- CellChat::CellChatDB.human
   # CellChat::subsetDB(search = "Secreted Signaling", key = "annotation")
  cellchat@DB <- use_db

  options(future.globals.maxSize = 3000 * 1024^2)
  
  # Using parallel futures makes it actually slower - uses too much memory and some overheads?
  #with(future::plan(future::multisession, workers = n_cores), {
    cellchat |> 
      CellChat::subsetData() |> 
      CellChat::identifyOverExpressedGenes() |> 
      CellChat::identifyOverExpressedInteractions(variable.both = TRUE) |> 
      CellChat::smoothData(adj = CellChat::PPI.human)
  #})
}

compute_cellchat <- function(cellchat, n_cores = 8, n_boot = 100, type = "triMean", raw_use = TRUE, population_size = FALSE,
                             tm_trim = 0.1, scale_distance = 0.1, interaction_range = 100, contact_range = 10,
                             min_cells = 10, seed = 42) {
  
  options(future.globals.maxSize = 4000 * 1024^2)

  #with(future::plan(future::multisession, workers = n_cores), {
    CellChat::computeCommunProb(
      cellchat,
      nboot = n_boot,
      type = type,
      trim = tm_trim,
      raw.use = raw_use,
      population.size = population_size,
      distance.use = TRUE,
      interaction.range = interaction_range,
      scale.distance = scale_distance,
      contact.dependent = TRUE,
      contact.range = contact_range,
      seed.use = seed
    ) |> 
      CellChat::filterCommunication(min.cells = min_cells) |> 
      CellChat::computeCommunProbPathway() |> 
      CellChat::aggregateNet()
  #})

}



plot_cellchat_circle <- function(cellchat, pfile) {

  png(filename = file.path("fig", pfile), width = 1200, height = 600)
  par(mfrow = c(1, 2), xpd = TRUE)

  cnt <- cellchat@net$count
  CellChat::netVisual_circle(
    cnt,
    vertex.weight = rowSums(cnt),
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Number of interactions"
  )
  
  wei <- cellchat@net$weight
  CellChat::netVisual_circle(
    wei,
    vertex.weight = rowSums(wei),
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Interaction weights/strength"
  )
  dev.off()
}

plot_cellchat_heatmap <- function(cellchat, pfile, measure = "weight") {

  png(filename = file.path("fig", pfile), width = 800, height = 800)
  par(mfrow = c(1, 2))

  CellChat::netVisual_heatmap(
    cellchat,
    measure = measure,
    color.heatmap = "Blues",
    font.size = 12
  ) |> 
    print()
  
  dev.off()
}

plot_cellchat_rank <- function(cellchat) {

  env <- new.env(parent =  globalenv())

  env$d1 <- CellChat::rankNet(cellchat, mode = "single", measure = "count", return.data = TRUE)$signaling.contribution |> as_tibble() 
  env$d2 <- CellChat::rankNet(cellchat, mode = "single", measure = "weight", return.data = TRUE)$signaling.contribution |> as_tibble()

  env$plot_one_ <- function(d, what) {
    d |>
      ggplot(aes(x = contribution, y = name)) +
      th +
      geom_col() +
      scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
      labs(x = what, y = NULL)
  }

  with(env, {
    pl1 <- plot_one_(d1, "Number of interactions")
    pl2 <- plot_one_(d2, "Information flow")

    plot_grid(pl1, pl2, nrow = 1)
  })
}

plot_cellchat_bubble <- function(cellchat, p_value_limit = 0.01) {
  CellChat::netVisual_bubble(
    cellchat,
    thresh = p_value_limit,
    remove.isolate = TRUE
  )
}

plot_cellchat_bubble_sel <- function(cellchat, prob_limit = 0.01, point_size = 2) {
  env <- new.env(parent = globalenv())
  env$point_size <- point_size

  env$d <- cellchat |>
    CellChat::subsetCommunication() |>
    as_tibble() |> 
    filter(prob >= prob_limit) |> 
    filter(source != "Other" & target != "Other") |> 
    unite(link, c(source, target), sep = "->")

  with(env, {
    d |> 
      ggplot(aes(x = link, y = interaction_name_2, colour = prob)) +
      th +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      ) +
      geom_point(size = point_size) +
      scale_colour_viridis_c(option = "cividis") +
      labs(x = NULL, y = NULL, colour = "Probability")
  })
  
}


plot_cellchat_chord <- function(cellchat, pathway, pfile) {

  png(filename = file.path("fig", pfile), width = 800, height = 800)
  par(mfrow = c(1, 1))

  CellChat::netVisual_chord_cell(
    cellchat,
    signaling = pathway,
    lab.cex = 1.5
  ) |> 
    print()
  
  dev.off()
}

plot_cellchat_chords <- function(cellchat, sample, paths, prefix = "cellchat_chord") {
  for(path in paths) {
    pfile = str_glue("{prefix}_{sample}_{path}.png")
    plot_cellchat_chord(cellchat, path, pfile)
  }
}


plot_cellchat_network <- function(cellchat, pathway, pfile) {

  png(filename = file.path("fig", pfile), width = 800, height = 800)
  par(mfrow = c(1, 1))

  CellChat::netAnalysis_signalingRole_network(
    cellchat,
    signaling = pathway,
    width = 12,
    height = 3,
    font.size = 12
  ) |> 
    print()
  
  dev.off()
}


plot_cellchat_networks <- function(cellchat, sample, paths, prefix = "cellchat_network") {
  cellchat <- CellChat::netAnalysis_computeCentrality(cellchat)
  for(path in paths) {
    pfile = str_glue("{prefix}_{sample}_{path}.png")
    plot_cellchat_network(cellchat, path, pfile)
  }
}