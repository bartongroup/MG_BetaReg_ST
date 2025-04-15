
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

get_gene_expression <- function(seu, genes, norm_per_gene = FALSE) {
  d <- GetAssayData(seu, layer = "data")[genes, , drop = FALSE] |> 
    t() |> 
    as_tibble(rownames = "cell_id") |> 
    pivot_longer(-cell_id, names_to = "gene_symbol", values_to = "expression")
  if(norm_per_gene) {
    d <- d |> 
      group_by(gene_symbol) |> 
      mutate(expression = expression / max(expression)) |> 
      ungroup()
  }
  return(d)
}


get_group_gene_expression <- function(seu, df, group_var = "cell_type") {
  d <- df |> 
    mutate(group = get(group_var))

  groups <- d$group |> unique()
  cell_order <- seu@meta.data |> 
    rownames()

  map(groups, function(gr) {
    gns <- d |> 
      filter(group == gr) |> 
      pull(gene_symbol)
    get_gene_expression(seu, gns, norm_per_gene = TRUE) |> 
      summarise(expression = mean(expression), .by = cell_id) |> 
      add_column(group = gr)
  }) |> 
    list_rbind() |> 
    pivot_wider(id_cols = cell_id, names_from = group, values_from = expression) 
}

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
  v <- as.character(v)
  str_glue("{suffix}_snn_res.{v}")
}



get_meta_column <- function(seu, col) {
  seu@meta.data |> 
    as_tibble(rownames = "cell_id") |> 
    select(cell_id, value = all_of(col))
}

# Algorithm = 4 selects Leiden algorithm. Leiden improves on Louvain by providing a guarantee of connected
# clusters and generally produces higher-quality clusters by resolving issues such as disconnected or
# fragmented clusters that Louvain can sometimes create.
cluster_seurat <- function(seu, dims = 1:40, resolutions = c(0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.3, 0.4),
                           ident_res = 0.2, seed = 42, res_suffix = "Spatial.008um", nfeatures = 3000,
                           min_group_count = 3, algorithm = 4) {
  seu <- seu |> 
    FindVariableFeatures(nfeatures = nfeatures) |> 
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


find_cluster_markers <- function(seu, genes = NULL) {
  #DefaultAssay(seu) <- "Spatial.008um"
  df <- seu |> 
    FindAllMarkers(test.use = "wilcox", features = genes) 
  if(nrow(df) > 0) {
    df |> 
      as_tibble() |> 
      rename(gene_symbol = gene) |> 
      relocate(gene_symbol, .before = 1)
  } else {
    NULL
  }
    
}

get_top_cluster_markers <- function(clm, n_top = 12, fdr_limit = 0.01) {
  clm |> 
    filter(p_val_adj < fdr_limit) |> 
    group_by(cluster) |> 
    top_n(n_top, avg_log2FC) |> 
    ungroup()
}


get_cluster_cell_count <- function(seu, resolution) {
  res <- get_resolution_columns(resolution)
  seu@meta.data |>
    group_by(cluster = get(res)) |>
    tally()
}


find_gene_correlation_cor <- function(seu, sel_genes, assay = "Spatial.008um", method = "spearman", min_overlap = 100, n_top = 3000) {

  over_genes <- get_overlapping_genes(seu, sel_genes, min_overlap)
  top_genes <- get_top_genes(seu, assay, n_top)
  good_genes <- intersect(over_genes, top_genes) |> union(sel_genes)

  expr_data <- GetAssayData(seu, layer = "data")[good_genes, ]

  ex_all <- expr_data |> 
    t() |> 
    as.matrix()
  ex_sel <- ex_all[, sel_genes]

  cor(ex_sel, ex_all, method = method) |> 
    as_tibble(rownames = "gene_1") |> 
    pivot_longer(-gene_1, names_to = "gene_symbol", values_to = "correlation") |> 
    filter(gene_1 != gene_symbol)
}


find_gene_correlation_cscore <- function(seu, sel_genes, assay = "Spatial.008um", min_overlap = 100, n_top = 3000) {
  # Find overlap between sel_genes and all other genes
  over_genes <- get_overlapping_genes(seu, sel_genes, min_overlap)
  top_genes <- get_top_genes(seu, assay, n_top)
  good_genes <- intersect(over_genes, top_genes) |> union(sel_genes)

  # CSCORE searches for assay "RNA", so we need to rename
  seur <- seu |> 
    RenameAssays(assay, "RNA")
  cs <- CSCORE::CSCORE(seur, good_genes)

  cr <- cs$est |> 
    as_tibble(rownames = "gene_1") |> 
    pivot_longer(-gene_1, names_to = "gene_symbol", values_to = "correlation") |> 
    filter(gene_1 %in% sel_genes)
  pv <- cs$p_value |> 
    as_tibble(rownames = "gene_1") |> 
    pivot_longer(-gene_1, names_to = "gene_symbol", values_to = "p_value") |> 
    filter(gene_1 %in% sel_genes)
  full_join(cr, pv, by = c("gene_1", "gene_symbol")) |> 
    mutate(p_adj = p.adjust(p_value, method = "BH")) |>    # adjust on selection, other tests ignored
    arrange(desc(correlation))
}


annotate_clusters <- function(seu, key_markers, resolution = NULL) {
  if(!is.null(resolution)) {
    res <- get_resolution_columns(resolution)
    Idents(seu) <- res
  }

  all_genes <- Features(seu)
  cell_types <- unique(key_markers$cell_type)
  markers <- map(cell_types, function(ct) {
    genes <- key_markers |> 
      filter(cell_type == ct & gene_symbol %in% all_genes) |> 
      pull(gene_symbol)
    cm <- find_cluster_markers(seu, genes)
    if(!is.null(cm)) {
      cm <- cm |> add_column(cell_type = ct, .before = 1)
    }
    cm
  }) |> 
    list_rbind()

  markers |> 
    group_by(cell_type, cluster) |> 
    summarise(log_fc = mean(avg_log2FC)) |> 
    pivot_wider(id_cols = cluster, names_from = cell_type, values_from = log_fc)
}


save_for_loupe <- function(seu, fname) {
  seu@meta.data |>
    as_tibble(rownames = "barcode") |> 
    select(barcode, cell_type, starts_with("Spatial.008um_snn_res")) |> 
    rename_with(~str_replace(., "^Spatial.008um_snn_res.", "resolution_"), starts_with("Spatial.008um_snn_res")) |> 
    write_csv(fname)
}
