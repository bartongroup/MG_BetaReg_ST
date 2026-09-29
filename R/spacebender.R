#' Export Space Ranger raw HD counts, including outside-tissue bins
#'
#' Reads the raw feature-barcode HDF5 matrix and tissue_positions.parquet from
#' outs/binned_outputs/square_008um (by default). Uses hdf5r, arrow, Matrix and
#' native anndataR; neither Seurat nor Python is required. All raw barcodes,
#' including zero-count bins, are retained. Only Gene Expression features are
#' exported. X is sparse bins x genes; obsm['spatial'] contains full-resolution
#' image coordinates in pixels (x = pixel column, y = pixel row).
#'
#' @param data_path Space Ranger outs directory (EXPERIMENTS$DATA_PATH).
#' @param sample Sample name used in spacebender_data/<sample>.h5ad.
#' @param bin_size HD bin size in micrometres; defaults to 8.
#' @param overwrite Whether to replace an existing H5AD export.
#' @return Written H5AD path, invisibly (suitable for a file target).
save_data_for_spacebender <- function(data_path, sample, bin_size = 8,
                                     overwrite = TRUE) {
  if (length(bin_size) != 1L || !is.numeric(bin_size) || is.na(bin_size) ||
      !is.finite(bin_size) || bin_size < 1 || bin_size != floor(bin_size)) {
    stop("bin_size must be a positive integer.")
  }
  bin_path <- file.path(data_path, "binned_outputs", sprintf("square_%03dum", as.integer(bin_size)))
  matrix_file <- file.path(bin_path, "raw_feature_bc_matrix.h5")
  positions_file <- file.path(bin_path, "spatial", "tissue_positions.parquet")
  h5ad_file <- file.path("spacebender_data", paste0(sample, ".h5ad"))
  for (path in c(matrix_file, positions_file)) {
    if (!file.exists(path)) stop("Required Space Ranger file not found: ", path)
  }
  if (file.exists(h5ad_file) && !overwrite) {
    stop("File already exists: ", h5ad_file, ". Use overwrite = TRUE to replace it.")
  }

  message("Reading raw counts: ", matrix_file)
  h5 <- hdf5r::H5File$new(matrix_file, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  shape <- as.numeric(h5[["matrix/shape"]][])
  barcodes <- as.character(h5[["matrix/barcodes"]][])
  gene_ids <- as.character(h5[["matrix/features/id"]][])
  gene_names <- as.character(h5[["matrix/features/name"]][])
  feature_types <- as.character(h5[["matrix/features/feature_type"]][])
  if (length(shape) != 2L || any(shape < 1) ||
      !identical(shape, as.numeric(c(length(gene_names), length(barcodes)))) ||
      length(gene_ids) != shape[1L] || length(feature_types) != shape[1L] ||
      anyNA(barcodes) || any(!nzchar(barcodes)) || anyDuplicated(barcodes) ||
      anyNA(gene_names) || any(!nzchar(gene_names)) || anyNA(feature_types)) {
    stop("Invalid matrix dimensions, barcodes or feature metadata in raw HDF5 file.")
  }
  # Match Read10X_h5's unique symbols and Seurat's feature-name substitutions,
  # while retaining the original symbols and stable gene IDs in var.
  feature_names <- gsub("[_|]", "-", make.unique(gene_names))
  if (anyDuplicated(feature_names)) {
    stop("Gene names collide after Seurat-compatible name conversion.")
  }
  p <- as.numeric(h5[["matrix/indptr"]][])
  if (length(p) != shape[2L] + 1 || any(p > .Machine$integer.max) ||
      any(shape > .Machine$integer.max)) {
    stop("Matrix dimensions/pointers are invalid or exceed R Matrix's integer limits.")
  }
  # The 10x arrays store a CSC genes x bins matrix, with zero-based indices.
  counts <- Matrix::sparseMatrix(
    i = as.integer(h5[["matrix/indices"]][]), p = as.integer(p),
    x = as.numeric(h5[["matrix/data"]][]), index1 = FALSE,
    dims = as.integer(shape), dimnames = list(feature_names, barcodes)
  )
  if (any(!is.finite(counts@x)) || any(counts@x < 0) ||
      any(counts@x != floor(counts@x))) stop("Raw matrix contains invalid UMI counts.")
  keep <- feature_types == "Gene Expression"
  if (!any(keep)) stop("No Gene Expression features found.")
  counts <- counts[keep, , drop = FALSE]
  features <- data.frame(
    gene_ids = gene_ids[keep], gene_name = gene_names[keep],
    feature_types = feature_types[keep], row.names = feature_names[keep]
  )

  fields <- c("barcode", "in_tissue", "array_row", "array_col",
              "pxl_row_in_fullres", "pxl_col_in_fullres")
  positions <- as.data.frame(arrow::read_parquet(positions_file, col_select = tidyselect::all_of(fields)))
  if (anyNA(positions$barcode) || anyDuplicated(positions$barcode)) {
    stop("Tissue positions contain missing or duplicated barcodes.")
  }
  order <- match(barcodes, positions$barcode)
  if (anyNA(order)) stop("Some raw matrix barcodes are missing from tissue positions.")
  positions <- positions[order, , drop = FALSE]
  tissue <- positions$in_tissue
  if (anyNA(tissue) || !all(tissue %in% c(0, 1))) {
    stop("Tissue labels must be 0 or 1 for every exported barcode.")
  }
  obs <- data.frame(
    in_tissue = as.logical(tissue), array_row = positions$array_row,
    array_col = positions$array_col, row.names = barcodes
  )
  spatial <- as.matrix(positions[, c("pxl_col_in_fullres", "pxl_row_in_fullres")])
  storage.mode(spatial) <- "double"
  dimnames(spatial) <- NULL
  if (any(!is.finite(spatial))) stop("Spatial coordinates contain invalid values.")
  if (!any(obs$in_tissue)) stop("No in-tissue bins found in raw matrix.")
  if (all(obs$in_tissue)) warning("Raw matrix contains no outside-tissue bins.")
  zero_bins <- sum(Matrix::colSums(counts) == 0)
  if (zero_bins > 0L) {
    warning(zero_bins, " zero-count bins retained. Current SpaceBender posterior code ",
            "may fail on an entirely zero-count batch; check the runner before a full run.")
  }
  adata <- anndataR::AnnData(
    X = Matrix::t(counts), obs = obs, var = features,
    obsm = list(spatial = spatial),
    uns = list(spatial_coordinate_units = "full-resolution image pixels", bin_size_um = bin_size)
  )
  dir.create(dirname(h5ad_file), recursive = TRUE, showWarnings = FALSE)
  adata$write_h5ad(h5ad_file, compression = "gzip", mode = if (overwrite) "w" else "w-")
  message("Wrote ", ncol(counts), " bins (", sum(obs$in_tissue), " in tissue; ",
          sum(!obs$in_tissue), " outside tissue) x ", nrow(counts), " genes to ", h5ad_file)
  invisible(h5ad_file)
}


#' Export raw spatial counts for SpaceBender using native anndataR
#'
#' Exports every bin and gene in the selected assay. To export a crop, subset
#' the Seurat object spatially before calling this function. No Python is used.
#' Existing in_tissue labels are retained; if absent, all bins are assumed to
#' be in tissue (appropriate for a filtered Seurat object).
#'
#' @param seu Seurat object with a single counts layer and a spatial image.
#' @param sample Sample name used in spacebender_data/<sample>.h5ad.
#' @param assay Assay containing raw counts.
#' @param image Spatial image name; required if the assay has multiple images.
#' @param overwrite Whether to overwrite an existing export.
#' @return Path to the written H5AD file, invisibly.
save_data_for_spacebender_from_seurat <- function(seu, sample, assay = "Spatial.008um",
                                     image = NULL, overwrite = FALSE) {
  h5ad_file <- file.path("spacebender_data", paste0(sample, ".h5ad"))
  if (file.exists(h5ad_file) && !overwrite) {
    stop("File already exists: ", h5ad_file, ". Use overwrite = TRUE to replace it.")
  }
  layers <- SeuratObject::Layers(seu[[assay]], search = NA)
  if (!identical(grep("^counts($|\\.)", layers, value = TRUE), "counts")) {
    stop("Expected a single 'counts' layer; join split counts layers first.")
  }
  counts <- SeuratObject::LayerData(seu, assay = assay, layer = "counts")
  counts <- methods::as(counts, "CsparseMatrix")
  if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)) ||
      any(!is.finite(counts@x)) || any(counts@x < 0) ||
      any(counts@x != floor(counts@x))) {
    stop("Counts must be non-negative integers with unique gene and bin names.")
  }
  if (is.null(image)) {
    images <- SeuratObject::Images(seu, assay = assay)
    if (length(images) != 1L) {
      stop("Specify image explicitly: expected exactly one image for ", assay, ".")
    }
    image <- images[[1L]]
  }
  coords <- SeuratObject::GetTissueCoordinates(seu[[image]])
  # VisiumV2 supplies an explicit cell column; older images use row names.
  coord_bins <- if ("cell" %in% names(coords)) as.character(coords$cell) else rownames(coords)
  if (!all(c("x", "y") %in% names(coords)) || anyDuplicated(coord_bins) ||
      !all(colnames(counts) %in% coord_bins)) {
    stop("Spatial coordinates must have x/y columns and cover every exported bin.")
  }
  spatial <- as.matrix(coords[match(colnames(counts), coord_bins), c("x", "y"), drop = FALSE])
  storage.mode(spatial) <- "double"
  if (any(!is.finite(spatial))) stop("Spatial coordinates contain missing or infinite values.")
  dimnames(spatial) <- NULL
  meta <- seu[[]][colnames(counts), , drop = FALSE]
  if ("in_tissue" %in% names(meta)) {
    tissue <- meta$in_tissue
    if (anyNA(tissue) || !(is.logical(tissue) || is.numeric(tissue)) ||
        !all(tissue %in% c(0, 1))) {
      stop("in_tissue must contain logical or numeric 0/1 values.")
    }
    tissue <- as.logical(tissue)
  } else {
    tissue <- rep(TRUE, ncol(counts))
    message("No in_tissue labels: treating all exported bins as in tissue.")
  }
  adata <- anndataR::AnnData(
    X = Matrix::t(counts),
    obs = data.frame(in_tissue = tissue, row.names = colnames(counts)),
    var = data.frame(gene_name = rownames(counts), row.names = rownames(counts)),
    obsm = list(spatial = spatial)
  )
  dir.create(dirname(h5ad_file), recursive = TRUE, showWarnings = FALSE)
  adata$write_h5ad(h5ad_file, compression = "gzip", mode = if (overwrite) "w" else "w-")
  message("Wrote ", ncol(counts), " bins x ", nrow(counts), " genes to ", h5ad_file)
  invisible(h5ad_file)
}

#' Read SpaceBender counts into a fresh spatial Seurat object
#'
#' Reads spacebender_data/<sample>_cleaned.h5ad. Requires exactly the same
#' genes and bins as the original assay, but accepts any ordering. For a crop,
#' pass the same cropped object used for export. Images, bin metadata, project
#' name and active identities are retained. The returned object has only the
#' selected assay and a counts layer: other assays, normalized/scaled data,
#' variable features, reductions, graphs and analysis history are not copied.
#' Metadata annotations (e.g. cluster, cell_type and scores) remain annotations
#' of the ORIGINAL data until explicitly recomputed downstream.
#'
#' @param seu Original Seurat object used for export; never modified in place.
#' @param sample Sample name, without the _cleaned suffix.
#' @param assay Assay to replace with corrected counts.
#' @param mito_genes Mitochondrial feature names, preferably the same list used
#'   by load_visium_data(). If NULL, use human gene symbols starting with MT-.
#' @return A new counts-only Seurat object ready for normalization.
read_spacebender_clean_file <- function(seu, sample, assay = "Spatial.008um",
                                       mito_genes = NULL) {
  h5ad_file <- file.path("spacebender_data", paste0(sample, "_cleaned.h5ad"))
  if (!file.exists(h5ad_file)) stop("File does not exist: ", h5ad_file)
  layers <- SeuratObject::Layers(seu[[assay]], search = NA)
  if (!identical(grep("^counts($|\\.)", layers, value = TRUE), "counts")) {
    stop("Expected a single 'counts' layer; join split counts layers first.")
  }
  original <- SeuratObject::LayerData(seu, assay = assay, layer = "counts")
  genes <- rownames(original)
  bins <- colnames(original)
  if (!setequal(bins, colnames(seu))) {
    stop("The selected assay must cover every bin in the original object.")
  }
  adata <- anndataR::read_h5ad(h5ad_file, as = "InMemoryAnnData")
  if (anyDuplicated(adata$obs_names) || anyDuplicated(adata$var_names) ||
      !setequal(bins, adata$obs_names) || !setequal(genes, adata$var_names)) {
    stop("Cleaned file must have exactly the original gene and bin identifiers. ",
         "For a crop, supply the cropped Seurat object.")
  }
  # Use the storage-class conversion: direct dgRMatrix -> dgCMatrix
  # coercion is unavailable in some Matrix versions. This stays sparse.
  counts <- methods::as(Matrix::t(adata$X), "CsparseMatrix")
  dimnames(counts) <- list(adata$var_names, adata$obs_names)
  counts <- counts[genes, bins, drop = FALSE]
  if (any(!is.finite(counts@x)) || any(counts@x < 0) ||
      any(counts@x != floor(counts@x))) {
    stop("Cleaned X must contain finite, non-negative integer counts.")
  }
  # Subtraction-based SpaceBender correction must never increase raw counts.
  # Check in blocks to avoid allocating a second full difference matrix.
  for (start in seq.int(1L, length(bins), by = 1024L)) {
    j <- seq.int(start, min(start + 1023L, length(bins)))
    if (any(counts[, j, drop = FALSE] > original[, j, drop = FALSE])) {
      stop("Cleaned counts exceed original counts; check the input sample and assay.")
    }
  }
  totals <- Matrix::colSums(counts)
  detected <- Matrix::colSums(counts > 0)
  meta <- seu[[]][bins, , drop = FALSE]
  meta[[paste0("nCount_", assay)]] <- totals
  meta[[paste0("nFeature_", assay)]] <- detected
  if ("total_count" %in% names(meta)) meta$total_count <- totals
  if ("n_genes" %in% names(meta)) meta$n_genes <- detected
  if (is.null(mito_genes)) mito_genes <- grep("^MT-", genes, value = TRUE)
  mito <- intersect(mito_genes, genes)
  if (!length(mito) && "mito_fraction" %in% names(meta)) {
    warning("No mitochondrial features matched; mito_fraction will be 0 for non-empty bins. ",
            "Supply the original mito_genes list if features are not human gene symbols.")
  }
  fraction <- Matrix::colSums(counts[mito, , drop = FALSE]) / totals
  fraction[totals == 0] <- NA_real_
  meta$mito_fraction <- fraction
  if ("percent.mt" %in% names(meta)) meta$percent.mt <- 100 * fraction
  clean_assay <- if (inherits(seu[[assay]], "Assay5")) {
    SeuratObject::CreateAssay5Object(counts = counts)
  } else {
    SeuratObject::CreateAssayObject(counts = counts)
  }
  SeuratObject::Key(clean_assay) <- SeuratObject::Key(seu[[assay]])
  result <- SeuratObject::CreateSeuratObject(
    counts = clean_assay, assay = assay,
    project = SeuratObject::Project(seu), meta.data = meta
  )
  SeuratObject::Idents(result) <- SeuratObject::Idents(seu)[bins]
  for (image in SeuratObject::Images(seu, assay = assay)) {
    result[[image]] <- seu[[image]]
  }
  message("Loaded corrected counts from ", h5ad_file,
          ". Rerun normalization and downstream analyses; retained metadata labels/scores ",
          "still refer to the original analysis.")
  result
}


