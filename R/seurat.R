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
    as_tibble(rownames = "cell_id") |> 
    mutate(pcw = fct_reorder(pcw, week))
}

#' Extract gene data from a seurat object 
#'
#' @param seu A seurat object
#' @param group_var Grouping variable
#'
#' @return A tibble with gene detection counts, grouped per variable
get_gene_detection <- function(seu, group_var = "pcw") {
  mtx <- GetAssayData(seu, layer = "counts")
  b2g <- set_names(seu@meta.data[[group_var]], rownames(seu@meta.data))
  groups <- b2g[colnames(mtx)] |> 
    as.character()
  # Create numeric detection matrix (0 or 1), grouped column sums
  X <- DelayedArray::colsum(as(mtx > 0, "dMatrix"), groups)
  # Remove zero-only genes
  X[rowSums(X) > 0, ] |> 
    as_tibble(rownames = "gene_id")
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
get_average_expression <- function(seu, assay = "SCT", layer = "data", group_var = "pcw") {
  AverageExpression(seu, assays = assay, layer = layer, group.by = group_var) |>
    pluck(assay) |> 
    as_tibble(rownames = "feature_id")
}

get_mean_pct <- function(seu, assay = "SCT", layer = "data", group_var = "pct") {
  
  mt <- seu@meta.data |> 
    as_tibble(rownames = "cell_id") |> 
    mutate(cluster = get(group_var)) |> 
    select(cell_id, cluster) |> 
    mutate(
      cluster = as_factor(cluster),
      iclust = as.integer(cluster)
    )
  groups <- mt |> 
    select(cluster, iclust) |> 
    distinct()
  
  GetAssayData(seu, layer = layer, assay = assay) |> 
    t() |> 
    as.data.frame() |> 
    rownames_to_column("cell_id") |> 
    as_tibble() |> 
    pivot_longer(-cell_id, names_to = "feature_id") |> 
    left_join(mt, by = "cell_id") |> 
    group_by(feature_id, iclust) |> 
    summarise(
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      n_expr = length(which(value > 0)),
      n_cells = n(),
      pct = n_expr / n_cells
    ) |> 
    left_join(groups, by = "iclust") |> 
    ungroup()
}

get_resolution_columns <- function(v, suffix = "SCT") {
  str_glue("{suffix}_snn_res.{v}")
}


#' Create initial seurat object
#'
#' @param mtx An object of class `dMatrix`, a sparse matrix containing counts:
#'   rows are genes, columns are cells
#' @param design A tibble with design of the experiment, needs columns `cell_id`
#'   and `sample`
#' @param mito_genes A vector with mitochondrial gene IDs
#' @param transpose Logical, if true mtx will be transposed before creating
#'   Surat object; necessary for Alevin-Fry data
#'
#' @return A seurat object
create_seurat <- function(mtx, design, mito_genes, transpose = TRUE) {
  # Create object, matrix is transposed by comparison with Alevin data
  if(transpose)
    mtx <- t(mtx)
  
  seu <- CreateSeuratObject(mtx)
  
  # Need to select only mito genes present in the matrix, otherwise the next step fails
  mito_genes <- intersect(mito_genes, rownames(mtx))
  
  # Calculate mitochondrial fraction
  seu$mito_fraction <- PercentageFeatureSet(object = seu, features = mito_genes) / 100
  
  # Per cell
  cell_data <- seu@meta.data |> 
    as_tibble(rownames = "cell_id") |> 
    rename(orig_ident = orig.ident, total_count = nCount_RNA, n_genes = nFeature_RNA) |> 
    left_join(design, by = "cell_id") |>
    column_to_rownames("cell_id")
  
  seu@meta.data <- cell_data 
  return(seu)
}

filter_seurat <- function(seu, params) {

  sel_cells <- seu@meta.data |> 
    filter(
      total_count > params$min_counts_per_cell &
        n_genes > params$min_detected_genes_per_cell &
        mito_fraction < params$max_mt_fraction
    ) |>
    rownames()
  
  sel_genes <- get_gene_detection(seu) |> 
    pivot_longer(-gene_id, names_to = "pcw", values_to = "detections") |> 
    filter(detections > params$min_detections_per_gene) |> 
    pull(gene_id) |> 
    unique()
  
  subset(seu, cells = sel_cells, features = sel_genes)
}


check_cell_cycle <- function(seu, cell_cycle_markers, n_variable_features = 3000) {
  seu |> 
    NormalizeData() |>
    CellCycleScoring(s.features = s_genes, g2m.features = g2m_genes) |> 
    FindVariableFeatures(selection.method = "vst", nfeatures = n_variable_features) |> 
    ScaleData() |> 
    RunPCA() |> 
    DimPlot(reduction = "pca", group.by = "Phase", split.by = "Phase")
}


normalise_seurat <- function(seu, cell_cycle_markers, params) {
  all_genes <- rownames(seu)
  g2m_genes <- cell_cycle_markers |> 
    filter(phase == "G2/M") |> 
    pull(feature_id) |> 
    intersect(all_genes)
  s_genes <- cell_cycle_markers |> 
    filter(phase == "S") |> 
    pull(feature_id) |> 
    intersect(all_genes)
  
  options(future.globals.maxSize = 16000 * 1024^2)
  
  # Perform SCT normalisation for each condition separately
  seu |>
    NormalizeData() |> 
    CellCycleScoring(s.features = s_genes, g2m.features = g2m_genes) |> 
    SCTransform(vars.to.regress = c("mito_fraction", "Phase"), variable.features.n = params$n_variable_features)
}


split_normalise_seurat <- function(seu, cell_cycle_markers, params) {
  all_genes <- rownames(seu)
  g2m_genes <- cell_cycle_markers |> 
    filter(phase == "G2/M") |> 
    pull(geneID) |> 
    intersect(all_genes)
  s_genes <- cell_cycle_markers |> 
    filter(phase == "S") |> 
    pull(geneID) |> 
    intersect(all_genes)
  
  options(future.globals.maxSize = 16000 * 1024^2)
  
  # Perform SCT normalisation for each condition separately
  seu |> 
    SplitObject(split.by = "pcw") |> 
    map(function(sr) {
      sr |> 
        NormalizeData() |>
        CellCycleScoring(s.features = s_genes, g2m.features = g2m_genes) |> 
        # SCTransform performs normalization variance stabilisation using regularised negative binomial model
        SCTransform(vars.to.regress = c("mito_fraction", "Phase"), variable.features.n = params$n_variable_features,
                    vst.flavor = "v2")
    })
}


integrate_seurat <- function(seu_split, params) {
  # Select most variable genes found by CSTransform
  integ_features <- SelectIntegrationFeatures(seu_split, nfeatures = params$n_variable_features)
  # Prepare split object for integration
  seu_split <- PrepSCTIntegration(seu_split, anchor.features = integ_features)
  # Perform canonical correlation analysis (CCA) to identify sources of variation between conditions.
  # Identify anchors (or mutual nearest neighbours or best buddies) across conditions,
  # Filter incorrect anchors
  # WARNING: this part takes a lot of time to calculate!
  seu_split <- FindIntegrationAnchors(seu_split, normalization.method = "SCT", anchor.features = integ_features)
  
  # Integrate dataset using anchors
  IntegrateData(anchorset = seu_split, normalization.method = "SCT")
}


cluster_seurat <- function(seu_integ, dims = 1:40, resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.4),
                           ident_res = 0.6, seed = 42, res_suffix = "SCT", min_group_count = 3) {
  seu <- seu_integ |> 
    RunPCA() |> 
    FindNeighbors(dims = dims) |> 
    FindClusters(resolution = resolutions, random.seed = seed)
  
  # Count cells in each cluster/pcw combination and filter out small
  barcodes <- seu@meta.data |>
    as_tibble(rownames = "barcode") |> 
    select(barcode, pcw, contains("snn_res")) |> 
    pivot_longer(-c(barcode, pcw), names_to = "resolution", values_to = "cluster") |> 
    group_by(pcw, resolution, cluster) |> 
    mutate(n = n()) |> 
    filter(n >= min_group_count) |> 
    pull(barcode)
  seu <- subset(seu, cells = barcodes)
    
  # Set selected resolution as Ident (will be used by default)
  res <- get_resolution_columns(ident_res, suffix = res_suffix)
  Idents(seu) <- res
  seu@meta.data$cluster <- seu[[res]][[1]]
  seu |> 
    BuildClusterTree() |> 
    RunUMAP(dims = dims, seed.use = seed)
}

find_all_markers <- function(seu, idents, assay = "SCT", layer = "data", test_use = "wilcox") {
  DefaultAssay(seu) <- assay
  if(idents == "pcw")
    Idents(seu) <- idents
  seu |> 
    FindAllMarkers(test.use = test_use, layer = layer, only.pos = TRUE) |> 
    as_tibble() |> 
    rename(feature_id = gene)
}

find_conserved_markers <- function(seu_clust, assay = "RNA", layer = "data", resolution = NULL,
                                   resolution_suffix = "SCT",
                                   marker_suffixes = c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")) {
  expr <- GetAssayData(seu_clust, layer = layer, assay = assay) |> 
    rowMeans()
  tb <- tibble(
    feature_id = names(expr),
    mean_expr = expr
  )
  
  if(!is.null(resolution))
    Idents(seu_clust) <- get_resolution_columns(resolution, resolution_suffix)
  clusters <- Idents(seu_clust) |> unique()
  
  # Suffixes produced by FindConservedMarkers
  suffix_pattern <- str_c(marker_suffixes, collapse = "|")
  regex_pattern <- str_glue("^(.+)_({suffix_pattern})$")

  map_dfr(clusters, function(clst) {
    seu_clust |> 
      FindConservedMarkers(assay = assay, layer = layer, test.use = "wilcox", ident.1 = clst, grouping.var = "pcw") |> 
      as_tibble(rownames = "feature_id") |> 
      add_column(cluster = clst)
  }) |>
    select(-c(max_pval, minimump_p_val)) |> 
    pivot_longer(-c(feature_id, cluster), names_pattern = regex_pattern, names_to = c("pcw", "quantity")) |>
    pivot_wider(id_cols = c("feature_id", "cluster", "pcw"), names_from = "quantity", values_from = "value") |> 
    left_join(tb, by = "feature_id")
}

top_marker_genes <- function(cons, max_n = 12) {
  cons |> 
    group_by(cluster, feature_id) |> 
    summarise(
      mfc = mean(avg_log2FC)
    ) |>
    top_n(max_n, mfc) |> 
    arrange(desc(mfc), .by_group = TRUE) |> 
    ungroup()
}



find_de_markers <- function(seu, contrasts, assay = "SCT", layer = "data", ...) {
  Idents(seu) <- "pcw"
  mt <- seu@meta.data |> 
    as_tibble(rownames = "cell_id")
  
  map_dfr(contrasts, function(ctr) {
    sam <- str_split(ctr, "-") |> 
      unlist()
    cells <- mt |> 
      filter(pcw %in% sam) |> 
      pull(cell_id)
    expr <- GetAssayData(seu, layer = layer, assay = assay)[, cells] |> 
      rowMeans()
    tb <- tibble(
      gene_id = names(expr),
      mean_expr = expr
    )
    
    seu |> 
      PrepSCTFindMarkers() |> 
      FindMarkers(ident.1 = sam[1], ident.2 = sam[2], assay = assay, layer = layer, ...) |> 
      as_tibble(rownames = "gene_id") |> 
      left_join(tb, by = "gene_id") |> 
      add_column(contrast = ctr)
  })
}


get_cross_genes <- function(conserved_markers, de_markers, fdr_limit = 0.01) {
  cm <- conserved_markers |> 
    filter(p_val_adj < fdr_limit) |> 
    mutate(score = scale(avg_log2FC))
  dem <- de_markers |> 
    filter(p_val_adj < fdr_limit) |> 
    mutate(score = scale(avg_log2FC))
  
  inner_join(cm, dem, by = "gene_id", multiple = "all")
}


correlate_gene <- function(seu, gid, min_n_overlap = 100, layer = "data", assay = "RNA") {
  tab <- GetAssayData(seu, assay = assay, layer = layer)
  gene <- tab[gid, ]
  cell_sel <- which(gene > 0)
  gene <- gene[cell_sel]
  tab <- tab[, cell_sel]
  n_over <- rowSums(tab > 0)
  gene_sel <- which(n_over > min_n_overlap)
  tab <- as.matrix(tab[gene_sel, ])
  tab[tab == 0] <- NA
  r <- cor(t(tab), gene, use = "pairwise.complete.obs", method = "pearson")
  tibble(
    gene_id = rownames(r),
    r = as.vector(r),
    n_over = n_over[gene_sel]
  ) |> 
    filter(gene_id != gid)
}


get_cluster_cell_count <- function(seu, resolution) {
  res <- get_resolution_columns(resolution)
  seu@meta.data |>
    group_by(cluster = get(res)) |>
    tally()
}
