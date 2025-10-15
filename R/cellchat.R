#' Create a CellChat object from a Seurat spatial dataset
#'
#' @description
#' Wraps CellChat::createCellChat to initialize a CellChat object
#' using normalized expression data, metadata, and spatial coordinates extracted
#' from a Seurat object.
#'
#' @param seu A Seurat object containing spatial transcriptomics data.
#' @param spot_diameter_fullres Numeric. The full-resolution diameter of each capture spot
#'   (in the same units as the spatial coordinates). Default is \code{17.56}.
#' @param label_var Character. Name of the metadata column in \code{seu@meta.data}
#'   to use for grouping cells (e.g., cell type). Default is \code{"cell_type"}.
#' @param spot_size Numeric. Effective diameter to use for downsampling spatial
#'   spots when computing spatial factors. Default is \code{8}.
#' @param sub_n_cells Integer. Number of cells to randomly subsample from the full
#'   dataset (currently not active). Default is \code{5000}.
#'
#' @return A CellChat object ready for downstream analysis.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Extracts normalized expression data from the \code{data} layer.
#'   \item Builds metadata including the specified grouping variable.
#'   \item Computes spatial scaling factors based on \code{spot_diameter_fullres}
#'         and \code{spot_size}.
#'   \item Retrieves x–y coordinates via Seurat::GetTissueCoordinates.
#'   \item Calls CellChat::createCellChat with \code{datatype="spatial"}.
#' }
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

#' Preprocess a CellChat object
#'
#' @description
#' Subsets the CellChat database to human interactions, identifies over-expressed
#' genes and interactions, and smooths the data using the human PPI network.
#'
#' @param cellchat A CellChat object created by
#'   \code{create_cellchat()}.
#' @param n_cores Integer. Number of cores to use for parallel processing
#'   (currently unused). Default is \code{8}.
#'
#' @return The input \code{CellChat} object with:
#' \itemize{
#'   \item \code{@DB} set to \code{CellChat::CellChatDB.human}.
#'   \item Over-expressed genes and interactions identified.
#'   \item Expression data smoothed over the human PPI network.
#' }
#'
#' @details
#' This function increases the maximum future globals size to accommodate
#' large matrices, but does not currently run in parallel due to memory
#' overhead.
preprocess_cellchat <- function(cellchat, n_cores = 8) {
  use_db <- CellChat::CellChatDB.human
   # CellChat::subsetDB(search = "Secreted Signaling", key = "annotation")
  cellchat@DB <- use_db

  options(future.globals.maxSize = 10 * 1024^3)
  
  # Using parallel futures makes it actually slower - uses too much memory and some overheads?
  #with(future::plan(future::multisession, workers = n_cores), {
    cellchat |> 
      CellChat::subsetData() |> 
      CellChat::identifyOverExpressedGenes() |> 
      CellChat::identifyOverExpressedInteractions(variable.both = TRUE) |> 
      CellChat::smoothData(adj = CellChat::PPI.human)
  #})
}

#' Compute communication probabilities and aggregate network in CellChat
#'
#' @description
#' Runs bootstrapped communication probability calculations, filters low-support
#' interactions, computes pathway-level communications, and aggregates the network.
#'
#' @param cellchat A preprocessed \code{CellChat::CellChat} object.
#' @param n_cores Integer. Number of cores for parallel futures (unused). Default is \code{8}.
#' @param n_boot Integer. Number of bootstrap replicates for communication probability. Default is \code{100}.
#' @param type Character. Statistical measure for aggregation (e.g., \code{"triMean"}). Default is \code{"triMean"}.
#' @param raw_use Logical. Whether to use raw expression values. Default is \code{TRUE}.
#' @param population_size Logical. Whether to normalize by population size. Default is \code{FALSE}.
#' @param tm_trim Numeric. Trim fraction when computing trimmed mean. Default is \code{0.1}.
#' @param scale_distance Numeric. Scale factor for spatial distances. Default is \code{0.1}.
#' @param interaction_range Numeric. Maximum distance for signaling interactions. Default is \code{100}.
#' @param contact_range Numeric. Maximum distance for contact-dependent signaling. Default is \code{10}.
#' @param min_cells Integer. Minimum number of cells required to retain an interaction. Default is \code{10}.
#' @param seed Integer. Random seed for reproducibility. Default is \code{42}.
#'
#' @return A \code{CellChat} object with computed communication
#'   probabilities, filtered interactions, pathway probabilities, and an aggregated
#'   signaling network.
#'
#' @details
#' Adjusts \code{future.globals.maxSize} for large data objects. Runs
#' \code{CellChat::computeCommunProb} with contact-dependent smoothing,
#' then filters by \code{min_cells}, calculates pathway-level probabilities,
#' and aggregates the final network via \code{CellChat::aggregateNet}.
compute_cellchat <- function(cellchat, n_cores = 8, n_boot = 100, type = "triMean", raw_use = TRUE, population_size = FALSE,
                             tm_trim = 0.1, scale_distance = 0.1, interaction_range = 100, contact_range = 10,
                             min_cells = 10, seed = 42) {
  
  options(future.globals.maxSize = 10 * 1024^3)

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


#' Run CellChat on a selected subset of cells
#'
#' @description
#' Filters a Seurat object to a given set of cells, removes rare cell types
#' (fewer than `min_n` cells), and runs a CellChat workflow on the remaining
#' cells. Returns the fitted `CellChat` object and a tibble of
#' ligand–receptor communication probabilities ranked by `prob`.
#'
#' @param seu A **Seurat** object containing single-cell data. Must include a
#'   `cell_type` column in `seu@meta.data`.
#' @param cells A character vector of cell/barcode IDs to retain prior to
#'   filtering by cell type size.
#' @param min_n Integer (default `3`). Minimum number of cells required for a
#'   cell type to be kept. Cell types with fewer than `min_n` cells are removed
#'   before running CellChat.
#' @param scale_distance Numeric (default `0.15`). Distance scaling parameter
#'   forwarded to `compute_cellchat()` (used when spatial information is
#'   available).
#' @param n_boot Integer (default `100`). Number of bootstraps for communication
#'   probability estimation in `compute_cellchat()`.
#'
#' @return
#' A list with two elements:
#' \itemize{
#'   \item `cellchat`: A fitted `CellChat` object (or `NULL` if not enough cell
#'   types remain after filtering).
#'   \item `prob`: A tibble of inferred communications returned by
#'   `CellChat::subsetCommunication()`, sorted by descending `prob`
#'   (or `NULL` if CellChat was not run).
#' }
cellchat_on_selection <- function(seu, cells, min_n = 3, scale_distance = 0.15, n_boot = 100) {
  seu_sel <- subset(seu, cells = cells)
  
  
  meta <- seu_sel@meta.data |> 
    as_tibble(rownames = "cell_id")

  # Remove cell types with too few cells
  bad_levels <- meta |> 
    count(cell_type) |>
    filter(n < min_n) |>
    pull(cell_type)
  message(str_glue("Removing {length(bad_levels)} cell types with fewer than {min_n} cells: {str_c(bad_levels, collapse = ', ')}"))
  good_cells <- meta |> 
    filter(!(cell_type %in% bad_levels)) |>
    pull(cell_id)
  seu_sel <- subset(seu_sel, cells = good_cells)
  seu_sel$cell_type <- seu_sel$cell_type |> droplevels()
  if(length(unique(seu_sel$cell_type)) < 2) {
    message("Not enough cell types remaining after filtering, skipping CellChat")
    return(list(cellchat = NULL, prob = NULL))
  }

  cmp <- create_cellchat(seu_sel) |>
    preprocess_cellchat() |> 
    compute_cellchat(n_boot = n_boot, scale_distance = scale_distance, raw_use = TRUE, population_size = FALSE)
  prb <- cmp |> 
    CellChat::subsetCommunication() |>
    as_tibble() |>
    arrange(desc(prob))
  list(cellchat = cmp, prob = prb)
}


cellchat_insulin_gradient <- function(seu, bins = 6, min_n = 3, scale_distance = 0.15, n_boot = 100) {
  gx <- gene_binned_gradient(seu, "INS", bins, include_zeroes = FALSE) |> 
    drop_na() |> 
    droplevels()
  all_bins <- levels(gx$bin)
  map(all_bins, function(b) {
    message(str_glue("*** Processing bin {b}"))
    cells <- gx |> filter(bin == b) |> pull(cell_id)
    cellchat_on_selection(seu, cells, min_n = min_n, scale_distance = scale_distance, n_boot = n_boot) |> 
        pluck("cellchat")
  }) |> 
    set_names(all_bins)
}

plot_cellchat_insulin_gradient_ranks <- function(cc_list) {
  env <- new.env(parent = globalenv())

  cc_list <- compact(cc_list)
  bins <- names(cc_list)
  
  env$d <- map(bins, function(b) {
    CellChat::rankNet(cc_list[[b]], mode = "single", measure = "weight",
                        return.data = TRUE) |>
      pluck("signaling.contribution") |> 
      as_tibble() |> 
      add_column(bin = b)
  }) |> 
    list_rbind() |> 
    mutate(bin = factor(bin, levels = bins))
  with(env, {
    d |> ggplot(aes(x = contribution, y = name)) +
      th +
      geom_col() +
      scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
      facet_wrap(~bin, nrow = 1, scales = "free") +
      labs(x = "Information flow", y = NULL)
  })
}
 

plot_cellchat_circle <- function(cellchat, pfile) {

  png(filename = file.path("doc", "images", pfile), width = 1200, height = 600)
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

  png(filename = file.path("doc", "images", pfile), width = 800, height = 800)
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


plot_one_rank <- function(d, what) {
  d |>
    ggplot(aes(x = contribution, y = name)) +
    th +
    geom_col() +
    scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
    labs(x = what, y = NULL)
}

plot_cellchat_rank <- function(cellchat) {

  env <- new.env(parent =  globalenv())

  env$d1 <- CellChat::rankNet(cellchat, mode = "single", measure = "count", return.data = TRUE)$signaling.contribution |> as_tibble() 
  env$d2 <- CellChat::rankNet(cellchat, mode = "single", measure = "weight", return.data = TRUE)$signaling.contribution |> as_tibble()


  with(env, {
    pl1 <- plot_one_rank(d1, "Number of interactions")
    pl2 <- plot_one_rank(d2, "Information flow")

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

  png(filename = file.path("doc", "images", pfile), width = 800, height = 800)
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

  png(filename = file.path("doc", "images", pfile), width = 800, height = 800)
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