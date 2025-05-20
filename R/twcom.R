get_ligands_receptors <- function() {
  lrdb <- CellChat::CellChatDB.human
  ligands <- receptors <- vector(mode="list", length=nrow(lrdb$interaction))

  for(i in 1:length(ligands)) {
    
    if(!(lrdb$interaction$ligand[i] %in% rownames(lrdb$complex))) {
        ligands[[i]] <- lrdb$interaction$ligand[i]
    } else {
        ltemp <- c(lrdb$complex[rownames(lrdb$complex) == lrdb$interaction$ligand[i], ])
        ltemp <- ltemp[nzchar(ltemp)]
        ligands[[i]] <- ltemp
    }
    
    if(!(lrdb$interaction$receptor[i] %in% rownames(lrdb$complex))) {
        receptors[[i]] <- lrdb$interaction$receptor[i]
    } else {
        ltemp <- c(lrdb$complex[rownames(lrdb$complex) == lrdb$interaction$receptor[i], ])
        ltemp <- ltemp[nzchar(ltemp)]
        receptors[[i]] <- ltemp
    }
    
  }
  list(
    ligands = lapply(ligands, unlist),
    receptors = lapply(receptors, unlist)  
  )
}

run_fretcom <- function(seu, max_dist = 100, spot_diameter_fullres = 17.56, spot_size = 8) {
  stdat <- SeuratObject::GetAssayData(seu, layer = "data")
  cell_types <- seu@meta.data$cell_type |> unique() |> as.character()
  cell_types <- setdiff(cell_types, "Other")
  M <- seu@meta.data[, cell_types]
  coords <-  SeuratObject::GetTissueCoordinates(seu)
  lr <- get_ligands_receptors()
  xy_unit_um <- spot_diameter_fullres / spot_size

  TWCOM::FRETCOM(
    stdat = stdat,
    M = M,
    ligands = lr$ligands,
    receptors = lr$receptors,
    coordx = coords$x,
    coordy = coords$y,
    maxdist = max_dist * xy_unit_um,
    parallel = FALSE,
    ncores = 4
  ) |> 
    tibble::as_tibble()
}