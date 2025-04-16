#!/usr/bin/env Rscript
#$ -mods l_hard m_mem_free 64G

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if exactly two arguments are provided
if (length(args) != 2) {
  stop("Please provide exactly two arguments: seurat_file and out_file", call. = FALSE)
}

# Assign arguments to variables
seurat_file <- args[1]
out_file <- args[2]

# Check if files exist
if (!file.exists(seurat_file)) {
  stop(paste("File not found:", seurat_file), call. = FALSE)
}

####

get_overlapping_genes <- function(seu, genes, min_overlap) {
  x <- GetAssayData(seu, layer = "counts")
  x[x > 0] <- 1   # binary matrix, expression = 1, no expression = 0
  is_sel <- x[genes, ] |> colSums()
  is_sel[is_sel > 0] <- 1  # cells where of any sel_genes is expressed
  overlap <- x %*% is_sel  # matrix multiplication, returns overlap count between all genes and sel_genes

  # Only use genes with minimum overlap (number of cells)
  names(overlap[overlap[, 1] > min_overlap, ])
}

get_top_genes <- function(seu, assay, n_top) {
  AverageExpression(seu, layer = "data", group.by = "orig_ident") |>
    pluck(assay) |> 
    as_tibble(rownames = "gene_symbol") |> 
    slice_max(order_by = V1, n = n_top) |> 
    pull(gene_symbol)
}

find_gene_correlation_cscore <- function(seu, assay = "Spatial.008um", min_overlap = 100, n_top = 3000) {
  
  all_genes <- Features(seu)

  # For testing
  all_genes <- c("NKX6-1", "PDX1", "PTF1A", "SOX9", "NEUROD1")


  # CSCORE searches for assay "RNA", so we need to rename
  seur <- seu |> 
    RenameAssays(assay, "RNA")
  cs <- CSCORE::CSCORE(seur, all_genes)

  cr <- cs$est |> 
    as_tibble(rownames = "gene_1") |> 
    pivot_longer(-gene_1, names_to = "gene_2", values_to = "correlation")
  pv <- cs$p_value |> 
    as_tibble(rownames = "gene_1") |> 
    pivot_longer(-gene_1, names_to = "gene_2", values_to = "p_value")
  full_join(cr, pv, by = c("gene_1", "gene_2")) |> 
    mutate(p_adj = p.adjust(p_value, method = "BH")) |> 
    arrange(desc(correlation)) |> 
    filter(gene_1 <= gene_2)
}

library(Seurat)
library(CSCORE)
library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(arrow)


seu <- readRDS(seurat_file)

cr <- find_gene_correlation_cscore(seu, min_overlap = 10, n_top = Inf)

print(paste("Writing to", out_file))
arrow::write_parquet(cr, out_file)