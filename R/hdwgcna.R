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
    hdWGCNA::MetacellsByGroups(group.by = c(group_var), ident.group = group_var, reduction = "spatial", k = k, max_shared = max_shared)  |> 
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


plot_modules_expr_score <- function(seu, n_genes = "all", point_size = 0.4, ncol = 2) {
  library(UCell)

  xy <- get_spatial_coordinates(seu)
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
    geom_polygon(aes(x = x_corner, y = y_corner, group = cell_id, fill = value), colour = NA, linewidth = 0) +
    facet_wrap(~ module, ncol = ncol) +
    labs(x = NULL, y = NULL, fill = "Expr score") +
    scale_fill_distiller(palette = "Purples", direction = 1)
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
  spat <- seu@reductions$spatial@cell.embeddings |> 
    as_tibble(rownames = "cell_id")
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