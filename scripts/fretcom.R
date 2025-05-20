#!/usr/bin/env Rscript
#$ -mods l_hard m_mem_free 200G
#$ -jc long
#$ -pe smp 12

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if exactly three arguments are provided
if (length(args) != 3) {
  stop("Please provide exactly three arguments: seurat_file, ligand-receptor file and out_file", call. = FALSE)
}

# Assign arguments to variables
seurat_file <- args[1]
lr_file <- args[2]
out_file <- args[3]

# Check if files exist
if (!file.exists(seurat_file)) {
  stop(paste("File not found:", seurat_file), call. = FALSE)
}

if (!file.exists(lr_file)) {
  stop(paste("File not found:", lr_file), call. = FALSE)
}

####

run_fretcom <- function(seu, lr, n_batch = 5, max_dist = 100, spot_diameter_fullres = 17.56, spot_size = 8) {
  stdat <- SeuratObject::GetAssayData(seu, layer = "data")
  cell_types <- seu@meta.data$cell_type |> unique() |> as.character()
  cell_types <- setdiff(cell_types, "Other")
  M <- seu@meta.data[, cell_types]
  coords <-  SeuratObject::GetTissueCoordinates(seu)
  xy_unit_um <- spot_diameter_fullres / spot_size

  batch_cut <- cut(seq_along(lr$ligands), breaks = n_batch, labels = FALSE)
  ligands_batches <- split(lr$ligands, batch_cut)
  receptor_batches <- split(lr$receptors, batch_cut)

  i <- 1
  purrr::map2(ligands_batches, receptor_batches, function(lig, rec) {
    message(sprintf("====== Batch %i =======", i))
    i <<- i + 1
    TWCOM::FRETCOM(
      stdat = stdat,
      M = M,
      ligands = lig,
      receptors = rec,
      coordx = coords$x,
      coordy = coords$y,
      maxdist = max_dist * xy_unit_um,
      parallel = TRUE,
      ncores = 12
    ) |> 
      tibble::as_tibble()
  }) |> 
    purrr::list_rbind()
}

library(Seurat)
library(TWCOM)
library(tibble)
library(purrr)


seu <- readRDS(seurat_file)
lr <- readRDS(lr_file)

fret <- run_fretcom(seu, lr)

write.table(fret, file = out_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
