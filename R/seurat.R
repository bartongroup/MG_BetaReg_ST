
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

  bs <- ifelse(is.numeric(bin_size), sprintf("%03dum", bin_size), str_to_title(bin_size))
  suffix <- str_glue("Spatial.{bs}")
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

#' Load a feature (gene) list from TSV
#'
#' Reads a tab-separated file with columns gene_id, gene_symbol, and a trailing dummy column and returns symbols.
#'
#' @param ffile Character path to the TSV file.
#'
#' @return Character vector of gene symbols.
load_features <- function(ffile) {
  read_tsv(ffile, col_names = c("gene_id", "gene_symbol", "dummy")) |> 
    pull(gene_symbol)
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

#' Count cells per annotated cell type
#'
#' Tallies cells by the cell_type metadata column.
#'
#' @param seu Seurat object containing a cell_type column in metadata.
#'
#' @return Tibble with columns cell_type and n.
cell_type_count <- function(seu) {
  seu@meta.data |> 
    as_tibble() |> 
    count(cell_type)
}

#' Get cells with nonzero counts for selected features
#'
#' Returns a long table of cells expressing specific feature IDs (counts layer).
#'
#' @param seu Seurat object.
#' @param ids Character vector of feature IDs (typically gene symbols).
#'
#' @return Tibble with columns feature_id, cell_id, and count (>0 only).
get_gene_cells <- function(seu, ids) {
  LayerData(seu, layer = "counts", features = ids) |> 
    as.matrix() |> 
    as_tibble(rownames = "feature_id") |>
    pivot_longer(-feature_id, names_to = "cell_id", values_to = "count") |> 
    filter(count > 0)
}

#' Get spatial coordinates in Visium layout
#'
#' Extracts tissue coordinates, reorients to a Newcastle convention
#' (swaps x/y and inverts y), and adds cell_id.
#'
#' @param seu Seurat object with Visium data.
#'
#' @return Tibble with columns cell_id, x, and y.
get_spatial_coordinates <- function(seu, swap_xy = TRUE) {
  cr <- GetTissueCoordinates(seu) |> 
    as_tibble(rownames = "cell_id") |>
    select(-cell) |> 
    relocate(cell_id, .before = 1)
  if(swap_xy) {
    cr <- cr |> 
      rename(
        t = x,
        x = y
      ) |> 
      rename(y = t)
  }
  cr |> 
    mutate(y = -y)
}

#' Get segmentation polygons coordinates in Visium layout
#'
#' Extracts tissue coordinates, reorients to a Newcastle convention
#' (swaps x/y and inverts y), and adds cell_id.
#'
#' @param seu Seurat object with Visium data.
#'
#' @return Tibble with columns cell_id, x, and y.
get_polygons <- function(seu) {
  # This is ugly, I don't know a better way to access the polygons
  seg <- seu@images$slice1.polygons$segmentations
  GetTissueCoordinates(seg) |> 
    as_tibble(rownames = "cid") |>
    select(-cid, cell_id = cell) |> 
    # Rearrange cooridnates to look like Newcastle report
    #rename(
    #  t = x,
    #  x = y
    #) |> 
    #rename(y = t) |> 
    mutate(y = -y) |> 
    relocate(cell_id, .before = 1) 
}


#' Average expression per group.
#'
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
get_average_expression <- function(seu, assay = "Spatial.008um", layer = "data", group_var = "ident") {
  AggregateExpression(seu, assays = assay, group.by = group_var, return.seurat = TRUE) |> 
    GetAssayData(layer = layer) |> 
    as_tibble(rownames = "gene_symbol") |> 
    pivot_longer(-gene_symbol, names_to = group_var, values_to = layer)
}


#' Total counts per gene
#'
#' Sums counts across all cells and ranks genes by total count.
#'
#' @param seu Seurat object.
#'
#' @return Tibble with gene_symbol, count, and rank (descending by count).
get_total_count <- function(seu) {
  cnt <- rowSums(GetAssayData(seu, layer = "counts"))
  tibble(
    gene_symbol = names(cnt),
    count = cnt,
    rank = rank(-count)
  ) |> 
    arrange(rank)
}


#' Per-cell expression for selected genes
#'
#' Fetches expression for specified genes, optionally normalizing within each gene to [0,1].
#'
#' @param seu Seurat object.
#' @param genes Character vector of gene symbols to extract.
#' @param norm_per_gene Logical; if TRUE, divides values by each gene's maximum.
#' @param layer Character layer to use (e.g., "data" or "counts").
#'
#' @return Long tibble with columns cell_id, gene_symbol, and expression.
get_gene_expression <- function(seu, genes, norm_per_gene = FALSE, layer = "data") {
  d <- GetAssayData(seu, layer = layer)[genes, , drop = FALSE] |> 
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


#' Group-wise average expression of selected gene sets
#'
#' For each group defined in df[[group_var]], averages normalized expression across
#' its gene set per cell and pivots wide.
#'
#' @param seu Seurat object.
#' @param df Data frame containing at least gene_symbol and the grouping column.
#' @param group_var Column name in df that partitions genes into groups (default "cell_type").
#'
#' @return Wide tibble of per-cell average expression, one column per group.
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

#' Genes overlapping selected set by co-detection
#'
#' Finds genes co-detected with a selected gene set in at least min_overlap
#' cells (binary overlap on counts).
#'
#' @param seu Seurat object.
#' @param genes Character vector of selected genes.
#' @param min_overlap Integer minimum number of cells required for overlap.
#'
#' @return Character vector of gene symbols passing the overlap threshold.
get_overlapping_genes <- function(seu, genes, min_overlap) {
  x <- GetAssayData(seu, layer = "counts")
  x[x > 0] <- 1   # binary matrix, expression = 1, no expression = 0
  is_sel <- x[genes, ] |> colSums()
  is_sel[is_sel > 0] <- 1  # cells where of any sel_genes is expressed
  overlap <- x %*% is_sel  # matrix multiplication, returns overlap count between all genes and sel_genes

  # Only use genes with minimum overlap (number of cells)
  names(overlap[overlap[, 1] > min_overlap, ])
}

#' Top expressed genes across samples
#'
#' Selects the top n_top genes by average expression across orig_ident groups.
#'
#' @param seu Seurat object.
#' @param assay Character assay name used in AverageExpression.
#' @param n_top Integer number of top genes to return.
#'
#' @return Character vector of top gene symbols.
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


#' Build resolution metadata column names
#'
#' Generates Seurat cluster resolution column names for a given resolution value.
#'
#' @param v Numeric/character resolution value(s).
#' @param suffix Character prefix for the assay/resolution (default "Spatial.008um").
#'
#' @return Character vector of resolution column names like "<suffix>_snn_res.<v>".

get_resolution_columns <- function(v, suffix = "Spatial.008um") {
  v <- as.character(v)
  str_glue("{suffix}_snn_res.{v}")
}


#' Extract a metadata column as a tibble
#'
#' Returns the specified metadata column with cell_id.
#'
#' @param seu Seurat object.
#' @param col Single metadata column name to extract.
#'
#' @return Tibble with columns cell_id and value.
get_meta_column <- function(seu, col) {
  seu@meta.data |> 
    as_tibble(rownames = "cell_id") |> 
    select(cell_id, value = all_of(col))
}


#' Cluster a Seurat object across multiple resolutions
#'
#' Runs HVG selection, scaling, PCA, neighbor graph, and Leiden/Louvain clustering for 
#' multiple resolutions, sets the chosen ident, builds a cluster tree, and computes UMAP.
#'
#' Algorithm = 4 selects Leiden algorithm. Leiden improves on Louvain by providing a guarantee of connected
#' clusters and generally produces higher-quality clusters by resolving issues such as disconnected or
#' fragmented clusters that Louvain can sometimes create.
#' 
#' @param seu Seurat object.
#' @param dims Integer vector of PCs to use.
#' @param resolutions Numeric vector of clustering resolutions to compute.
#' @param ident_res Single resolution to set as active identity.
#' @param seed Integer random seed.
#' @param res_suffix Assay/resolution suffix used to build resolution column names.
#' @param nfeatures Number of variable features to select.
#' @param min_group_count Unused here; reserved for downstream logic.
#' @param algorithm Integer algorithm code for FindClusters.
#'
#' @return Seurat object with clustering across resolutions and UMAP embeddings.
cluster_seurat <- function(seu, dims = 1:40, resolutions = c(0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.3, 0.4),
                           ident_res = 0.2, seed = 42, res_suffix = "Spatial.008um", nfeatures = 3000,
                           min_group_count = 3, algorithm = 4) {
  seu <- seu |> 
    FindVariableFeatures(nfeatures = nfeatures) |> 
    ScaleData() |> 
    RunPCA(approx = FALSE) |> 
    FindNeighbors(dims = dims) |> 
    FindClusters(resolution = resolutions, random.seed = seed, algorithm = algorithm)
    
  # Set selected resolution as Ident (will be used by default)
  res <- get_resolution_columns(ident_res, suffix = res_suffix)
  Idents(seu) <- res
  seu@meta.data$cluster <- seu[[res]][[1]] |> as.character() |> as.numeric() |> factor()
  seu |> 
    BuildClusterTree() |> 
    RunUMAP(dims = dims, seed.use = seed)
}


#' Minimal preprocessing (no clustering)
#'
#' Normalizes, selects HVGs, scales, runs PCA and UMAP.
#'
#' @param seu Seurat object.
#' @param dims Integer vector of PCs to use.
#' @param seed Integer random seed for UMAP.
#' @param nfeatures Number of variable features to select.
#'
#' @return Seurat object with dimensional reductions computed.
process_seurat <- function(seu, dims = 1:40,
                           seed = 42, nfeatures = 3000) {
  seu <- seu |> 
    NormalizeData() |> 
    FindVariableFeatures(nfeatures = nfeatures) |> 
    ScaleData() |> 
    RunPCA(approx = FALSE) |> 
    RunUMAP(dims = dims, seed.use = seed)
}


#' Find markers for all groups
#'
#' Wrapper around FindAllMarkers with Wilcoxon test, optionally restricting to genes.
#'
#' @param seu Seurat object.
#' @param genes Optional character vector of candidate genes to test.
#' @param group_by Metadata field to define groups (default "cluster").
#'
#' @return Tibble of marker statistics with gene_symbol first, or NULL if empty.
find_group_markers <- function(seu, genes = NULL, group_by = "cluster") {
  #DefaultAssay(seu) <- "Spatial.008um"
  df <- seu |> 
    FindAllMarkers(test.use = "wilcox", features = genes, group.by = group_by) 
  if(nrow(df) > 0) {
    df |> 
      as_tibble() |> 
      rename(gene_symbol = gene) |> 
      relocate(gene_symbol, .before = 1)
  } else {
    NULL
  }
}


#' Select top markers per cluster
#'
#' Filters by adjusted p-value and keeps the top n_top by avg_log2FC per cluster.
#'
#' @param clm Tibble/data frame of cluster markers (e.g., from find_group_markers).
#' @param n_top Number of markers to keep per cluster.
#' @param fdr_limit Adjusted p-value cutoff.
#'
#' @return Tibble of top markers per cluster.
get_top_cluster_markers <- function(clm, n_top = 12, fdr_limit = 0.01) {
  clm |> 
    filter(p_val_adj < fdr_limit) |> 
    group_by(cluster) |> 
    top_n(n_top, avg_log2FC) |> 
    ungroup()
}

#' Count cells per cluster at a resolution
#'
#' Tallies cells by cluster for a given resolution column.
#'
#' @param seu Seurat object.
#' @param resolution Numeric/character resolution used to pick the metadata column.
#'
#' @return Tibble with columns cluster and n.
get_cluster_cell_count <- function(seu, resolution) {
  res <- get_resolution_columns(resolution)
  seu@meta.data |>
    group_by(cluster = get(res)) |>
    tally()
}


#' Gene–gene correlation (classical)
#'
#' Computes correlations between a selected gene set and all genes after filtering by co-detection 
#' overlap and global expression level.
#'
#' @param seu Seurat object.
#' @param sel_genes Character vector of selected genes (rows must exist).
#' @param assay Character assay name for average/top-gene filtering.
#' @param method Correlation method for cor (e.g., "spearman").
#' @param min_overlap Integer minimum co-detection cell count.
#' @param n_top Integer number of globally top expressed genes to include.
#'
#' @return Tibble with columns gene_1, gene_symbol, and correlation (pairs excluding self-pairs).
find_gene_correlation_cor <- function(seu, sel_genes, assay = "Spatial.008um", method = "spearman",
                             min_overlap = 100, n_top = 3000) {

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

#' Gene–gene correlation using CSCORE
#'
#' Runs CSCORE on the object (temporarily renaming the assay to "RNA") and returns correlations 
#' and BH-adjusted p-values.
#'
#' @param seu Seurat object.
#' @param assay Character assay name whose data should be fed to CSCORE.
#'
#' @return Tibble with gene_1, gene_2, correlation, p_value, and p_adj (upper triangle only).
find_gene_correlation_cscore <- function(seu, assay = "Spatial.008um") {
  
  all_genes <- Features(seu)

  # For testing
  # all_genes <- c("NKX6-1", "PDX1", "PTF1A", "SOX9", "NEUROD1")


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

#' Annotate clusters by key marker enrichment
#'
#' For each provided cell type, tests its key markers as group markers across clusters and
#' summarizes mean log2FC per cell type and cluster.
#'
#' @param seu Seurat object.
#' @param key_markers Tibble with columns cell_type and gene_symbol listing key markers.
#' @param resolution Optional resolution to set as ident before testing.
#'
#' @return Wide tibble with clusters as rows and cell types as columns (mean avg_log2FC).
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
    cm <- find_group_markers(seu, genes, group_by = "cluster")
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


#' Save metadata for Loupe Browser
#'
#' Writes a CSV with barcode, cell_type, and renamed cluster resolution columns.
#'
#' @param seu Seurat object.
#' @param fname Output CSV file path.
#'
#' @return Invisibly, the path written; called for its side effect of writing a file.
save_for_loupe <- function(seu, fname, cluster_col = "cluster") {
  seu@meta.data |>
    as_tibble(rownames = "barcode") |> 
    mutate(cluster = get(cluster_col)) |> 
    select(barcode, cell_type, cluster) |> 
    write_csv(fname)
}

#' Differential expression between two groups
#'
#' Runs FindMarkers between two metadata groups with a chosen test.
#'
#' @param seu Seurat object.
#' @param group_1 Character label for group 1.
#' @param group_2 Character label for group 2.
#' @param meta_col Metadata column to use for identities.
#'
#' @return Tibble of DE results with gene_symbol and a contrast column.
seurat_de <- function(seu, group_1 = "Beta", group_2 = "MPC", meta_col = "cell_type", method = "wilcox") {
  Idents(seu) <- meta_col

  FindMarkers(seu, ident.1 = group_1, ident.2 = group_2, test.use = method) |> 
    as_tibble(rownames = "gene_symbol") |> 
    add_column(contrast = str_glue("{group_1}-{group_2}"))
}

#' Expression percentage and average per group (DotPlot helper)
#'
#' Reproduces DotPlot's underlying summary (average expression and % detected) for given genes.
#'
#' @param seu Seurat object.
#' @param gns Character vector of gene symbols.
#' @param meta_column Metadata column defining groups (default "cell_type").
#'
#' @return Tibble with gene_symbol, avg_exp, pct_exp, and group.
get_expression_percentage <- function(seu, gns, meta_column = "cell_type") {
  Idents(seu) <- meta_column
  
  dp <- DotPlot(seu, features = gns)
  dp$data |> 
    as_tibble(rownames = "gid") |> 
    select(gene_symbol = features.plot, avg_exp = avg.exp, pct_exp = pct.exp, group = id)
}

#' Expression percentage using count threshold
#'
#' Computes per-group average counts and % of cells above a count threshold for each gene.
#'
#' @param seu Seurat object.
#' @param gns Character vector of gene symbols.
#' @param meta_column Metadata column defining groups.
#' @param expr_limit Numeric threshold on counts to call expression (default 0).
#'
#' @return Tibble with gene_symbol, avg_exp, pct_exp, and group ordered by group and gene.
get_expression_percentage_2 <- function(seu, gns, meta_column = "cell_type", expr_limit = 0) {
  ex <- FetchData(seu, vars = gns, layer = "counts") |> 
    as_tibble(rownames = "cell_id") |> 
    pivot_longer(-cell_id, names_to = "gene_symbol", values_to = "count")
  grp <- FetchData(seu, vars = meta_column) |> 
    as_tibble(rownames = "cell_id") |> 
    set_names(c("cell_id", "group"))
  
  ex |> 
    left_join(grp, by = "cell_id") |> 
    group_by(group, gene_symbol) |> 
    summarise(
      avg_exp = mean(count),
      pct_exp = 100 * length(which(count > expr_limit)) / n()
    )  |> 
    ungroup() |> 
    select(gene_symbol, avg_exp, pct_exp, group) |> 
    mutate(
      gene_symbol = factor(gene_symbol, levels = gns)
    ) |> 
    arrange(group, gene_symbol)
}

#' Identify unique hormone expressors
#'
#' Flags cells expressing exactly one of the classic pancreatic hormones (GCG, INS, SST, GHRL, PPY)
#' and labels with inferred cell type.
#'
#' @param seu Seurat object.
#'
#' @return Tibble with cell_id, gene_symbol, and cell_type for uniquely expressing cells.
unique_hormone_expressors <- function(seu) {
  gnt <- tribble(
    ~gene_symbol, ~cell_type,
    "GCG", "Alpha",
    "INS", "Beta",
    "SST", "Delta",
    "GHRL", "Epsilon",
    "PPY", "PP"
  )

  get_gene_expression(seu, gnt$gene_symbol) |>
    filter(expression > 0) |>
    group_by(cell_id) |>
    summarise(n = n(), gene_symbol = str_c(gene_symbol, collapse = ",")) |>
    filter(n == 1) |> 
    select(-n) |> 
    left_join(gnt, by = "gene_symbol")
}


#' Count cells by number of hormones expressed
#'
#' Counts how many hormones are expressed per cell and summarizes distribution as counts and percentages.
#'
#' @param seu Seurat object.
#' @param gns Character vector of hormone genes to consider (default five classical hormones).
#'
#' @return Tibble with n_hormones_expressed, n, and perc of cells.
get_hormone_count <- function(seu, gns = c("INS", "SST", "GCG", "GHRL", "PPY")) {
  get_gene_expression(seu, gns, layer = "data") |> 
    filter(expression > 0) |> 
    count(cell_id, name = "n_hormones_expressed") |> 
    count(n_hormones_expressed) |> 
    mutate(perc = 100 * n / sum(n))
}

#' Compare 8 µm vs 2 µm Visium expression in a region
#'
#' Extracts spots within a square region around (x_cen, y_cen) for both resolutions and binds
#' expression tables.
#'
#' @param seu_8um Seurat object for 8 µm assay.
#' @param seu_2um Seurat object for 2 µm assay.
#' @param genes Character vector of gene symbols to extract.
#' @param x_cen Numeric x-center of region.
#' @param y_cen Numeric y-center of region.
#' @param size Half-width/height of the square region (same units as coordinates).
#' @param layer Character layer to use for expression (default "counts").
#'
#' @return Tibble with coordinates, cell_id, gene_symbol, expression, and resolution.
combine_resolution_area <- function(seu_8um, seu_2um, genes, x_cen, y_cen, size, layer = "counts") {
  sp8 <- get_spatial_coordinates(seu_8um) |> 
    filter(abs(x - x_cen) < size & abs(y - y_cen) < size)
  sp2 <- get_spatial_coordinates(seu_2um) |> 
    filter(abs(x - x_cen) < size & abs(y - y_cen) < size)
  ex8 <- get_gene_expression(seu_8um, genes, layer = layer) |> 
    filter(cell_id %in% sp8$cell_id)
  ex2 <- get_gene_expression(seu_2um, genes, layer = layer) |> 
    filter(cell_id %in% sp2$cell_id)

  d8 <- inner_join(sp8, ex8, by = "cell_id") |> 
    add_column(resolution = "8um")
  d2 <- inner_join(sp2, ex2, by = "cell_id") |> 
    add_column(resolution = "2um")
  bind_rows(d8, d2)
}

#' Distribution of nonzero counts for a gene
#'
#' Summarizes the discrete count distribution among expressing cells and cumulative percentages.
#'
#' @param seu Seurat object.
#' @param gene Single gene symbol.
#'
#' @return Tibble with expression (count), n, pct_spots, and pct_cumul.
get_count_distribution <- function(seu, gene) {
  seu |> 
    get_gene_expression(gene, layer = "counts") |>
    filter(expression > 0) |>
    count(expression) |>
    mutate(
      pct_spots = 100 * n / sum(n),
      pct_cumul = cumsum(pct_spots)
    )
}

#' Overlap between clusters and annotated cell types
#'
#' For each cell_type, counts overlapping cells with each cluster and reports proportions
#' relative to cluster size.
#'
#' @param seu Seurat object containing cell_type and cluster in metadata.
#'
#' @return Tibble with cell_type, cluster, n_overlap, n_cluster, and prop.
cluster_cell_type_overlap <- function(seu) {
  m <- seu@meta.data |> 
    as_tibble() |> 
    select(cell_type, cluster)
  tot_ct <- m |> count(cell_type, name = "n_cell_type")
  tot_cl <- m |> count(cluster, name = "n_cluster")

  cts <- m$cell_type |> unique()
  map(cts, function(ct) {
    m |> 
      filter(cell_type == ct) |> 
      count(cluster, name = "n_overlap") |> 
      left_join(tot_cl, by = "cluster") |> 
      mutate(prop = n_overlap / n_cluster) |> 
      add_column(cell_type = ct)
  }) |> 
    list_rbind()
}


relabel_clusters <- function(seu, relc = NULL, resolution = "0.3") {
  if(is.null(relc))
    return(seu)
  
  res <- get_resolution_columns(resolution)
  r <- relc[as.integer(as.character(get_meta_column(seu, res)$value))]
  f <- fct_reorder(as.character(r), r)
  
  seu[[res]] <- f
  seu$cluster <- f
  seu
}


clean_seurat <- function(seu) {
  meta <- seu@meta.data
  remove_cols <- colnames(meta) |> 
    str_subset("Spatial|cluster")
  meta <- meta[, !(names(meta) %in% remove_cols)]
  CreateSeuratObject(
    counts = GetAssayData(seu, layer = "counts"),
    meta.data = meta
  ) |> 
  subset(nCount_RNA > 0)
}

integrate_seurats <- function(seu_1, seu_2, nfeatures = 2000, seed= 42) {
  seu_1$sample <- "sample1"
  seu_2$sample <- "sample2"

  options(future.globals.maxSize = 10 * 1024^3)

  seu_list <- list(seu_1, seu_2)

  seu_list <- map(seu_list, function(seu) {
    seu |> 
      clean_seurat() |> 
      NormalizeData() |> 
      SCTransform(seed.use = seed)
  })

  # Select most variable genes found by CSTransform
  integ_features <- SelectIntegrationFeatures(seu_list, nfeatures = nfeatures)
  # Prepare split object for integration
  seu_list <- PrepSCTIntegration(seu_list, anchor.features = integ_features)
  # Perform canonical correlation analysis (CCA) to identify sources of variation between samples.
  # Identify anchors (or mutual nearest neighbours or best buddies) across samples,
  # Filter incorrect anchors
  # WARNING: this part takes a lot of time to calculate!
  seu_list <- FindIntegrationAnchors(seu_list, normalization.method = "SCT", anchor.features = integ_features)
  
  # Integrate dataset using anchors
  IntegrateData(anchorset = seu_list, normalization.method = "SCT")
}


plot_seurat_clusters <- function(seu_clust, resolutions = NULL, with_legend = FALSE, label = TRUE, ...) {
  if(is.null(resolutions)) {
    resols <- names(seu_clust@meta.data) |> str_subset("snn_res")
  } else {
    resols <-  get_resolution_columns(resolutions, suffix = "integrated")
  }
  
  P <- map(resols, function(res) {
    Idents(seu_clust) <- res
    r <- str_remove(res, "^.*res\\.")
    dp <- DimPlot(seu_clust, reduction = "umap", label = label, label.size = 6, ...) +
      labs(title = r)
    if(!with_legend)
      dp <- dp + theme(legend.position = "none")
    dp
  })
  plot_grid(plotlist = P)
}
