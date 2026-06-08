#########################################################################################################
# As of Seurat 5.0, slot arguments have been dropped. WGSNA still uses slot arguments in some places.
# The following functions replace original hdWGCNA functions to be compatible with Seurat 5.


CheckWGCNAName <- function(seurat_obj, wgcna_name){
  check <- wgcna_name %in% names(seurat_obj@misc) 
  if(!check){
    stop(paste0("Invalid wgcna_name supplied: ", wgcna_name))
  }  
}

CheckSeurat5 <- function(){
    startsWith(as.character(packageVersion('Seurat')), '5')
}


SetMetacellObject <- function(seurat_obj, metacell_obj, wgcna_name=NULL){
  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  # add metacell obj to Seurat obj
  seurat_obj@misc[[wgcna_name]]$wgcna_metacell_obj <- metacell_obj
  seurat_obj
}

SetWGCNAParams <- function(seurat_obj, params, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$wgcna_params <- params
  seurat_obj
}



ConstructMetacells <- function(
  seurat_obj, name='agg', ident.group='seurat_clusters', k=25,
  reduction='pca', 
  dims = NULL,
  assay='RNA',
  cells.use = NULL, # if we don't want to use all the cells to make metacells, good for train/test split
  slot='counts',  
  layer='counts',
  meta=NULL, return_metacell=FALSE,
  mode = 'average', max_shared=15,
  target_metacells=1000,
  max_iter=5000,
  verbose=FALSE,
  wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

  # check reduction
  if(!(reduction %in% names(seurat_obj@reductions))){
    stop(paste0("Invalid reduction (", reduction, "). Reductions in Seurat object: ", paste(names(seurat_obj@reductions), collapse=', ')))
  }

  # check assay
  if(!(assay %in% names(seurat_obj@assays))){
    stop(paste0("Invalid assay (", assay, "). Assays in Seurat object: ", paste(names(seurat_obj@assays), collapse=', ')))
  }

  # check slot
  if(!(slot %in% c('counts', 'data', 'scale.data'))){
    stop(paste0("Invalid slot (", slot, "). Valid options for slot: counts, data, scale.data "))
  }

  # subset seurat object by selected cells:
  if(!is.null(cells.use)){
    seurat_full <- seurat_obj
    seurat_obj <- seurat_obj[,cells.use]
  }

  # get the dim reduction
  reduced_coordinates <- as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings)
  reduc_orig <- reduced_coordinates

  # subset the dims?
  if(!is.null(dims)){
    # are these indices?
    if(is.numeric(dims)){
      if(!all(dims %in% 1:ncol(reduced_coordinates))){
        stop("Invalid selection for dims. Check that the selected dims match the columns of the selected reduction.")
      }
    } else{
      if(!all(dims %in% colnames(reduced_coordinates))){
        stop("Invalid selection for dims. Check that the selected dims match the column names of the selected reduction.")
      }
    }
    reduced_coordinates <- reduced_coordinates[,dims]
  }

  # run KNN on the chosen dim reduction
  nn_map <- FNN::knn.index(reduced_coordinates, k = (k - 1))
  row.names(nn_map) <- row.names(reduced_coordinates)
  nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
  
  # set up variables for loop
  good_choices <- seq_len(nrow(nn_map))
  choice <- sample(seq_len(length(good_choices)), size = 1,
      replace = FALSE)
  chosen <- good_choices[choice]
  good_choices <- good_choices[good_choices != good_choices[choice]]
  it <- 0
  k2 <- k * 2
  get_shared <- function(other, this_choice) {
      k2 - length(union(cell_sample[other, ], this_choice))
  }

  # loop to create metacells until convergencce
  while (length(good_choices) > 0 & length(chosen) < target_metacells & it < max_iter) {
      it <- it + 1
      choice <- sample(seq_len(length(good_choices)), size = 1,
          replace = FALSE)
      new_chosen <- c(chosen, good_choices[choice])
      good_choices <- good_choices[good_choices != good_choices[choice]]
      cell_sample <- nn_map[new_chosen, ]
      others <- seq_len(nrow(cell_sample) - 1)
      this_choice <- cell_sample[nrow(cell_sample), ]
      shared <- sapply(others, get_shared, this_choice = this_choice)

      if(max(shared) <= max_shared){
          chosen <- new_chosen
      }
  }

  shared_old <- shared
  cell_sample <- nn_map[chosen, ]

   if(length(chosen) <= 1){
      warning('Metacell failed')
      return(NULL)
    }

  # get a list of the cell barcodes that have been merged:
  cells_merged <- apply(cell_sample, 1, function(x){
    paste0(colnames(seurat_obj)[x], collapse=',')
  })

  combs <- tryCatch(
    {combn(nrow(cell_sample), 2)},
    error = function(cond){return(NA)}
  )
  if(any(is.na(combs))){
    warning('Metacell failed')
    return(NULL)
  }

  shared <- apply(combs, 2, function(x) {
      k2 - length(unique(as.vector(cell_sample[x, ])))
  })

  if(verbose){
    message(paste0("Overlap QC metrics:\nCells per bin: ",
        k, "\nMaximum shared cells bin-bin: ", max(shared),
        "\nMean shared cells bin-bin: ", mean(shared), "\nMedian shared cells bin-bin: ",
        median(shared)))
    if (mean(shared)/k > 0.1)
        warning("On average, more than 10% of cells are shared between paired bins.")
  }

  # get original expression matrix
  if(CheckSeurat5()){
    exprs_old <- SeuratObject::LayerData(seurat_obj, assay=assay, layer=layer)
  } else{
    exprs_old <- Seurat::GetAssayData(seurat_obj, assay=assay, layer=slot)
  }

  # groups of cells to combine
  mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(exprs_old)) %in%
      cell_sample[x, , drop = FALSE])
  mask <- Matrix::Matrix(mask)

  # average or sum expression?
  new_exprs <- (exprs_old %*% mask)
  if(mode == 'average'){
    new_exprs <- new_exprs / k
  }
  colnames(new_exprs) <- paste0(name, '_', 1:ncol(new_exprs))
  rownames(cell_sample) <- paste0(name, '_', 1:ncol(new_exprs))
  colnames(cell_sample) <- paste0('knn_', 1:ncol(cell_sample))

  # make seurat obj:
  metacell_obj <- CreateSeuratObject(
    counts = new_exprs,
    assay = assay
  )
  if(slot == 'scale.data'){
    metacell_obj <- SeuratObject::SetAssayData(
      metacell_obj,
      layer=slot,
      assay=assay,
      new.data=as.matrix(new_exprs)
    )
  }

  # add the cells merged info to the metacell obj
  metacell_obj$cells_merged <- as.character(cells_merged)

  # calculate stats:
  max_shared <- max(shared)
  median_shared <- median(shared)
  mean_shared <- mean(shared)

  # calculate matrix density:
  new_exprs[new_exprs > 0] <- 1
  density <- sum(Matrix::colSums(new_exprs) / (nrow(new_exprs)*ncol(new_exprs)))
  run_stats <- data.frame(
    name = name,
    max_shared = max_shared,
    mean_shared = mean_shared,
    median_shared = median_shared,
    density = density,
    n = ncol(new_exprs)
  )

  # add to metacell seurat obj
  metacell_obj@misc$run_stats <- run_stats

  # add meta-data:
  if(!is.null(meta)){
    meta_names <- names(meta)
    for(x in meta_names){
      metacell_obj@meta.data[[x]] <- meta[[x]]
    }
  } else(
    warning('meta not found')
  )

  # add seurat metacell object to the main seurat object:
  if(return_metacell){
    out <- metacell_obj
  } else{

    # revert to full seurat object if we subsetted earlier
    if(!is.null(cells.use)){
      seurat_obj <- seurat_full
    }

    # add seurat metacell object to the main seurat object:
    seurat_obj <- SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)

    # add other info
    seurat_obj <- SetWGCNAParams(
      seurat_obj, params = list(
        'metacell_k' = k,
        'metacell_reduction' = reduction,
        'metacell_slot' = slot,
        'metacell_assay' = assay
      ),
      wgcna_name
    )

    out <- seurat_obj
  }
  out

}



MetacellsByGroups_fixed <- function(
  seurat_obj, group.by=c('seurat_clusters'),
  ident.group='seurat_clusters',
  k=25, 
  reduction='pca', 
  dims=NULL,
  assay=NULL,
  slot='counts', 
  layer='counts',
  mode = 'average', 
  cells.use = NULL, 
  min_cells=100,
  max_shared=15,
  target_metacells=1000,
  max_iter=5000, 
  verbose=FALSE, 
  wgcna_name=NULL
){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  # check group.by for invalid characters:
  if(any(grepl('#', group.by))){
    stop('Invalid character # found in group.by, please re-name the group.')
  }

  # check ident.group
  if(!(ident.group %in% group.by)){
    stop('ident.group must be in group.by')
  }

  # check mode
  if(!(mode %in% c('sum', 'average'))){
    stop('Invalid choice for mode. Mode can be either sum or average.')
  }

  # check reduction
  if(!(reduction %in% names(seurat_obj@reductions))){
    stop(paste0("Invalid reduction (", reduction, "). Reductions in Seurat object: ", paste(names(seurat_obj@reductions), collapse=', ')))
  }

  # check assay:
  if(is.null(assay)){
    assay <- DefaultAssay(seurat_obj)
  } else if(!(assay %in% names(seurat_obj@assays))){
    stop(paste0('Assay ', assay, ' not found in seurat_obj. Select a valid assay: ', paste0(names(seurat_obj@assays), collapse = ', ')))
  }

  # check slot/layer:
  if(!(slot %in% c('counts', 'data', 'scale.data'))){
    stop('Invalid input for slot. Valid choices are counts, data, scale.data.')
  } else{

    # check the shape of the slot
    slot_dim <- dim(Seurat::GetAssayData(seurat_obj, assay=assay, layer=slot))
    if(any(slot_dim) == 0){
      stop(paste(c("Selected slot ", slot, " not found in this assay.")))
    }
  }

  # check that k > min_cells 
  if(min_cells < k ){
    warning("min_cells is smaller than k, this may result in downstream errors if very small groups are allowed.")
  }

  # check that max_shared is valid:
  if(max_shared < 0){
    max_shared <- 0
    warning(paste0("max_shared specified (", max_shared, ') is too low, setting max_shared <- 0'))
  }

  # subset seurat object by seleted cells:
  if(!is.null(cells.use)){
    seurat_full <- seurat_obj
    seurat_obj <- seurat_obj[,cells.use]
  }

  # setup grouping variables
  if(length(group.by) > 1){
    seurat_meta <- seurat_obj@meta.data[,group.by]
    for(col in colnames(seurat_meta)){
      seurat_meta[[col]] <- as.character(seurat_meta[[col]])
    }
    seurat_obj$metacell_grouping <- apply(seurat_meta, 1, paste, collapse='#')
  } else {
    seurat_obj$metacell_grouping <- as.character(seurat_obj@meta.data[[group.by]])
  }
  groupings <- unique(seurat_obj$metacell_grouping)
  groupings <- groupings[order(groupings)]

  # remove groups that are too small:
  group_counts <- table(seurat_obj$metacell_grouping) < min_cells
  if(any(group_counts)){
    warning(paste0("Removing the following groups that did not meet min_cells: ", paste(names(group_counts)[group_counts], collapse=', ')))
  }
  groupings <- groupings[table(seurat_obj$metacell_grouping) >= min_cells]

  if(length(groupings) == 0 ){
    stop("No groups met the min_cells requirement.")
  }

  # unique meta-data for each group
  meta_df <- as.data.frame(do.call(rbind, strsplit(groupings, '#')))
  colnames(meta_df) <- group.by

  # list of meta-data to pass to each metacell seurat object
  meta_list <- lapply(1:nrow(meta_df), function(i){
    x <- list(as.character(meta_df[i,]))[[1]]
    names(x) <- colnames(meta_df)
    x
  })

  # split seurat obj by groupings
  seurat_list <- lapply(groupings, function(x){seurat_obj[,seurat_obj$metacell_grouping == x]})
  names(seurat_list) <- groupings

  # construct metacells
  metacell_list <- mapply(
    ConstructMetacells,
    seurat_obj = seurat_list,
    name = groupings,
    meta = meta_list,
    MoreArgs = list(
      k=k, 
      reduction=reduction, 
      dims=dims,
      assay=assay, 
      layer=layer,
      return_metacell=TRUE, 
      mode=mode,
      max_shared=max_shared, 
      max_iter=max_iter, 
      target_metacells=target_metacells, 
      verbose=verbose, 
      wgcna_name=wgcna_name
    )
  )
  names(metacell_list) <- groupings

  # remove NULL
  remove <- which(sapply(metacell_list, is.null))
  if(length(remove) > 1){
    metacell_list <- metacell_list[-remove]
  }

  # get the run stats:
  run_stats <- as.data.frame(do.call(rbind, lapply(metacell_list, function(x){x@misc$run_stats})))
  rownames(run_stats) <- 1:nrow(run_stats)
  for(i in 1:length(group.by)){
    run_stats[[group.by[i]]] <- do.call(rbind, strsplit(as.character(run_stats$name), '#'))[,i]
  }

  # combine metacell objects
  if(length(metacell_list) > 1){
    metacell_obj <- merge(metacell_list[[1]], metacell_list[2:length(metacell_list)])
    
    # need to join layers if this is Seurat 5
    if(CheckSeurat5()){
      metacell_obj <- SeuratObject::JoinLayers(metacell_obj)
    }
  } else{
    metacell_obj <- metacell_list[[1]]
  }

  # set idents for metacell object:
  Idents(metacell_obj) <- metacell_obj@meta.data[[ident.group]]

  # revert to full seurat object if we subsetted earlier
  if(!is.null(cells.use)){
    seurat_obj <- seurat_full
  }

  # add seurat metacell object to the main seurat object:
  seurat_obj <- SetMetacellObject(seurat_obj, metacell_obj, wgcna_name)

  # add other info
  param_list <- list(
    'metacell_k' = k,
    'metacell_reduction' = reduction,
    'metacell_assay' = assay,
    'metacell_stats' = run_stats
  )

  if(CheckSeurat5()){
    param_list['metacell_layer'] <- layer
  } else{
    param_list['metacell_slot'] <- slot
  }

  seurat_obj <- SetWGCNAParams(
    seurat_obj, params = param_list, wgcna_name
  )

  # return the updated seurat object
  seurat_obj
}


#########################################################################################################

hd_wgcna_analysis <- function(seu, sample, wgcna_name = "wgcna", group_var = "cluster", network_type = "signed hybrid",
                              assay = "Spatial.008um", layer = "data", soft_power = NULL, k = 25, max_shared = 15,
                              deep_split = 4, detect_cut_height = 0.995, min_module_size = 25, merge_cut_height = 0.2,
                              seed = 42) {
  
  sname <- sample |> as.character()

  # Need spatial coordinates in metadata, for MetacellsByGroups() to merge
  # spatially proximal groups
  coords <- SeuratObject::GetTissueCoordinates(seu) |> 
    select(spatial_1 = x, spatial_2 = y) |> 
    as.matrix()
  seu@reductions$spatial <- SeuratObject::CreateDimReducObject(
    embeddings = coords,
    assay = assay,
    key = 'spatial_'
  )

  seu |>
    hdWGCNA::SetupForWGCNA(wgcna_name = wgcna_name, gene_select = "variable") |> 
    MetacellsByGroups_fixed(group.by = c(group_var), ident.group = group_var, reduction = "spatial", k = k, max_shared = max_shared)  |> 
    hdWGCNA::NormalizeMetacells() |> 
    hdWGCNA::SetDatExpr(assay = assay, layer = layer) |> 
    hdWGCNA::TestSoftPowers(networkType = network_type) |> 
    hdWGCNA::ConstructNetwork(
      tom_name = sname,
      overwrite_tom = TRUE,
      soft_power = soft_power,
      networkType = network_type,
      deepSplit = deep_split,
      detectCutHeight = detect_cut_height,
      minModuleSize = min_module_size,
      mergeCutHeight = merge_cut_height,
      randomSeed = seed
    ) |> 
    hdWGCNA::ModuleEigengenes() |> 
    hdWGCNA::ModuleConnectivity()
}


get_modules_expr_score <- function(seu) {
  seu |> 
    hdWGCNA::ModuleExprScore(n_genes = "all", method = "UCell", ncores = 8) |> 
    hdWGCNA::GetModuleScores() |>
    as_tibble(rownames = "barcode")
}

plot_modules_expr_score <- function(seu, n_genes = "all", point_size = 0.4, ncol = 2) {
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

plot_modules_spatial <- function(seu, features = "hMEs", point_size = 0.4, ncol = 2) {
  pls <- hdWGCNA::ModuleFeaturePlot(seu, reduction = "spatial", features = features)
  map(pls, function(pl) {
    pl$data |> 
      as_tibble() |> 
      mutate(z = (val - min(val)) / (max(val) - min(val)))
  }) |> 
    list_rbind(names_to = "module") |> 
    ggplot(aes(x = spatial_2, y = -spatial_1, colour = z)) +
    thm +
    geom_point(size = point_size, stroke = 0.1, shape = 15) +
    scale_colour_distiller(palette = "Blues", direction = 1) +
    facet_wrap(~ module, ncol = ncol) +
    labs(x = NULL, y = NULL, colour = "Rel expr")
}

plot_modules_dots <- function(seu, tbl_hub_genes) {
  modules <- tbl_hub_genes$module |> unique() |> as.character()
  pl <- map(modules, function(mdl) {
    gns <- tbl_hub_genes |> 
      filter(module == mdl) |> 
      pull(gene_name)
    plot_gene_dots(seu, gns, meta_column = "cell_type", with_legend = FALSE, title = mdl)
  })
  plot_grid(plotlist = pl, ncol = 2)
}


plot_metacells_spatial <- function(seu, group_var = "cluster", point_size = 0.9) {
  env <- new.env(parent = globalenv())
  env$point_size <- point_size
  env$group_var <- group_var

  # Spatial coordinates of individual cells
  spat <-  get_spatial_coordinates(seu, swap_xy = TRUE) |> 
    mutate(y = -y) |>
    rename(spatial_1 = x, spatial_2 = y)
  #spat <- seu@reductions$spatial@cell.embeddings |> 
  #  as_tibble(rownames = "cell_id")
  # Get metacells
  mc <- hdWGCNA::GetMetacellObject(seu)
  env$d <- mc@meta.data |>
    as_tibble(rownames = "metacell_id") |> 
    mutate(group = get(group_var)) |> 
    separate_longer_delim(cells_merged, delim = ",") |> 
    left_join(spat, by = c("cells_merged" = "cell_id")) |> 
    group_by(metacell_id, group) |> 
    summarise(x = mean(spatial_2), y = -mean(spatial_1)) |> 
    ungroup()

  with(env, {
    d |> 
      ggplot(aes(x = x, y = y, colour = group)) +
      thm +
      geom_point(size = point_size) +
      scale_colour_brewer(palette = "Paired") +
      labs(x = NULL, y = NULL, colour = group_var)
  })

}

get_modules <- function(seu) {
  seu |> 
    hdWGCNA::GetModules() |> 
    as_tibble()
}

get_hub_genes <- function(seu, n_hubs = 30) {
  seu |> 
    hdWGCNA::GetHubGenes(n_hubs = n_hubs) |> 
    as_tibble()
}

plot_hub_gene_network <- function(seu, n_hubs = 12, future_globals_max_size = 10) {
  options(future.globals.maxSize = future_globals_max_size * 1024^3)
  seu |>
    hdWGCNA::HubGeneNetworkPlot(n_hubs = n_hubs, return_graph = TRUE)
}


module_enrichment <- function(tbl_modules, all_genes, fterms, kme_limit = 0, fdr_limit = 0.05) {
  tbl_modules |> 
    select(-color) |> 
    pivot_longer(-c(gene_name, module), names_to = "module_sel", values_to = "score") |> 
    mutate(module_sel = str_remove(module_sel, "kME_")) |> 
    filter(module == module_sel & score > kme_limit) |> 
    group_split(module) |> 
    map(function(w) {
      sel_genes <- w$gene_name |> unique()
      ontologies <- names(fterms)
      map(ontologies, function(ont) {
        df <- fenr::functional_enrichment(all_genes, sel_genes, fterms[[ont]])
        if(!is.null(df)) {
          if(nrow(df) > 0)
            df <- df |> 
              filter(p_adjust < fdr_limit) |> 
              add_column(module = w$module[1], ontology = ont, .before = 1)
        }
      }) |> 
        list_rbind()
    }) |> 
      list_rbind()
}


plot_hdwgcna_dendrogram <- function(seu_hdwgcna, suffix) {
  pfile <- file.path("doc", "images", str_glue("module_dendro_{suffix}.png"))
  png(pfile, width = 2000, height = 1000, pointsize = 32)
  
  hdWGCNA::PlotDendrogram(seu_hdwgcna)
  
  dev.off()
}