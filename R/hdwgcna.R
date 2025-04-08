hd_wgcna_analysis <- function(seu, sample, wgcna_name = "wgcna", group_var = "cell_type",
                              group_name = c("MPC", "EP", "Beta"), assay = "Spatial.008um", layer = "data") {
  
  sname <- sample |> as.character()

  coords <- as.matrix(GetTissueCoordinates(seu)[,c('x', 'y')])
  colnames(coords) <- c('spatial_1', 'spatial_2')
  seu@reductions$spatial <- CreateDimReducObject(
    embeddings = coords,
    assay = assay,
    key = 'spatial_'
  )

  seu |>
    SetupForWGCNA(wgcna_name = wgcna_name, gene_select = "variable") |> 
    MetacellsByGroups(group.by = group_var, ident.group = group_var, reduction = "spatial", k = 50, max_shared = 15)  |> 
    NormalizeMetacells() |> 
    SetDatExpr(assay = assay, layer = layer) |> 
    TestSoftPowers(networkType = "signed") |> 
    ConstructNetwork(tom_name = sname, overwrite_tom = TRUE) |> 
    ModuleEigengenes() |> 
    ModuleConnectivity()
}