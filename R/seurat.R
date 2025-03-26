#' Download cell cycle markers
#'
#' @param species A string with species name ("human")
#' @param g2s A tibble with gene_id and gene_symbol
#'
#' @return A tibble with markers
get_cell_cycle_markers <- function(species, g2s) {
  sp <- c(
    human = "Homo_sapiens"
  )
  g2s <- g2s |> 
    select(gene_id, gene_symbol) |> 
    distinct()
  
  remote_file <- str_glue("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/{sp[species]}.csv")
  read_csv(remote_file, show_col_types = FALSE) |>
    set_names(c("phase", "gene_id", "modifies")) |> 
    left_join(g2s, by = "gene_id")
}


#' Extract cell data from a seurat object
#'
#' @param seu A seurat object
#'
#' @return A tibble with cell data.
get_cell_data <- function(seu) {
  seu@meta.data |> 
    as_tibble(rownames = "cell_id")
}

get_gene_cells <- function(seu, ids) {
  LayerData(seu, layer = "counts", features = ids) |> 
    as.matrix() |> 
    as_tibble(rownames = "feature_id") |>
    pivot_longer(-feature_id, names_to = "cell_id", values_to = "count") |> 
    filter(count > 0)
}


#' Average expression per group.
#'
#' If layer is set to 'data', this function assumes that the data has been log
#' normalized and therefore feature values are exponentiated prior to averaging
#' so that averaging is done in non-log space. Otherwise, if layer is set to
#' either 'counts' or 'scale.data', no exponentiation is performed prior to
#' averaging
#'
#' @param seu Seurat object
#' @param assay Assay
#' @param layer layer
#' @param group_var Grouping variable
#'
#' @return
#' @export
#'
#' @examples
get_average_expression <- function(seu, assay = "Spatial.008um", layer = "count", group_var = "ident") {
  AverageExpression(seu, assays = assay, layer = layer, group.by = group_var) |>
    pluck(assay) |> 
    as_tibble(rownames = "feature_id")
}

#' Load and preprocess Visium spatial transcriptomics data
#'
#' This function loads 10X Visium spatial transcriptomics data from a specified directory,
#' calculates mitochondrial gene fractions, and prepares the metadata for downstream analyses.
#'
#' @param path Character string specifying the directory containing 10X Visium data.
#' @param mito_genes Character vector of mitochondrial gene names used to compute mitochondrial fractions.
#' @param bin_size Numeric value specifying the bin size (resolution) for spatial binning. Default is 8.
#'
#' @return A Seurat object with loaded data, mitochondrial fractions calculated, and standardized metadata.
load_visium_data <- function(path, mito_genes, bin_size = 8) {
  seu <- Load10X_Spatial(data.dir = path, bin.size = bin_size)

  suffix <- sprintf("Spatial.%03dum", bin_size)
  seu_genes <- Features(seu)
  mito <- intersect(mito_genes, seu_genes)
  seu$mito_fraction <- PercentageFeatureSet(object = seu, features = mito) / 100

  seu@meta.data <- seu@meta.data |> 
    as_tibble(rownames = "cell_id") |> 
    rename(orig_ident = orig.ident, total_count = str_glue("nCount_{suffix}"), n_genes = str_glue("nFeature_{suffix}")) |> 
    mutate(cell_type = "blank") |> 
    column_to_rownames("cell_id")

  return(seu)
}


#' Extract gene data from a seurat object 
#'
#' @param seu A seurat object
#' @param group_var Grouping variable
#' @param layer Seurat's object layer to extract
#'
#' @return A tibble with gene detection counts, grouped per variable
get_gene_detection <- function(seu, group_var = "cell_type", layer = "counts") {
  mtx <- GetAssayData(seu, layer = layer)
  b2g <- set_names(seu@meta.data[[group_var]], rownames(seu@meta.data))
  groups <- b2g[colnames(mtx)] |> 
    as.character()
  # Create numeric detection matrix (0 or 1), grouped column sums
  X <- DelayedArray::colsum(as(mtx > 0, "dMatrix"), groups)
  # Remove zero-only genes
  X[rowSums(X) > 0, , drop = FALSE] |> 
    as_tibble(rownames = "gene_symbol")
}



get_resolution_columns <- function(v, suffix = "Spatial.008um") {
  str_glue("{suffix}_snn_res.{v}")
}


# Algorithm = 4 selects Leiden algorithm. Leiden improves on Louvain by providing a guarantee of connected
# clusters and generally produces higher-quality clusters by resolving issues such as disconnected or
# fragmented clusters that Louvain can sometimes create.
cluster_seurat <- function(seu, dims = 1:40, resolutions = c(0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.3, 0.4),
                           ident_res = 0.2, seed = 42, res_suffix = "Spatial.008um", min_group_count = 3,
                           algorithm = 4) {
  seu <- seu |> 
    FindVariableFeatures() |> 
    ScaleData() |> 
    RunPCA() |> 
    FindNeighbors(dims = dims) |> 
    FindClusters(resolution = resolutions, random.seed = seed, algorithm = algorithm)
    
  # Set selected resolution as Ident (will be used by default)
  res <- get_resolution_columns(ident_res, suffix = res_suffix)
  Idents(seu) <- res
  seu@meta.data$cluster <- seu[[res]][[1]]
  seu |> 
    BuildClusterTree() |> 
    RunUMAP(dims = dims, seed.use = seed)
}


find_cluster_markers <- function(seu) {
  #DefaultAssay(seu) <- "Spatial.008um"
  seu |> 
    FindAllMarkers(test.use = "wilcox") |> 
    as_tibble(rownames = "gene_id")
}


