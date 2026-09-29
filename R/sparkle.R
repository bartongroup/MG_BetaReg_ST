#' Import SPARKLE counts into a copy of a segmented Seurat object
#'
#' Keeps original assays, images, reductions and metadata. Existing clusters and
#' reductions describe the original data. All cells in seu must be in the export;
#' additional exported cells are ignored. Features are matched through gene IDs
#' and Space Ranger's filtered feature names. Original counts must agree exactly
#' (within numerical tolerance). Fractional corrected counts are never rounded.
#'
#' @param seu Existing segmented Seurat object, e.g. targets::tar_read(seu_seg_1).
#' @param path Completed output directory from scripts/run_sparkle.py.
#' @param source_assay Original counts assay.
#' @param assay New assay name; existing assays are never overwritten.
#' @param normalize Apply LogNormalize to the new assay.
#' @param scale_factor LogNormalize scale factor.
#' @return A copy with the new default assay, cell/feature correction metadata,
#'   and provenance in object@misc$sparkle[[assay]]. Features match source_assay;
#'   additional exported features remain available in the files.
#' @export
read_sparkle <- function(seu, path, source_assay = "Spatial.Polygons",
                         assay = "SPARKLE", normalize = TRUE,
                         scale_factor = 10000) {
  for (pkg in c("Seurat", "SeuratObject", "Matrix", "jsonlite")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package: ", pkg)
  }
  if (!inherits(seu, "Seurat")) stop("seu must be a Seurat object.")
  if (!source_assay %in% names(seu@assays)) stop("Missing source assay: ", source_assay)
  if (length(assay) != 1L || is.na(assay) || !nzchar(assay)) stop("Invalid assay name.")
  if (assay %in% names(seu@assays)) stop("Refusing to overwrite assay: ", assay)
  files <- c("diagnostics.json", "corrected_counts.mtx.gz", "original_counts.mtx.gz",
             "features.tsv.gz", "barcodes.tsv.gz", "reference_features.tsv.gz",
             "gene_diagnostics.tsv.gz")
  absent <- files[!file.exists(file.path(path, files))]
  if (length(absent)) stop("Incomplete SPARKLE export; missing: ", paste(absent, collapse = ", "))
  diagnostics <- jsonlite::fromJSON(file.path(path, "diagnostics.json"))
  if (!identical(diagnostics$status, "complete") || !identical(diagnostics$schema_version, 1L) ||
      !isTRUE(diagnostics$original_counts_verified)) stop("Unverified or unsupported SPARKLE export.")
  read_tab <- function(name, header = FALSE) {
    con <- gzfile(file.path(path, name), "rt")
    on.exit(try(close(con), silent = TRUE))
    utils::read.delim(con, header = header, stringsAsFactors = FALSE,
                     quote = "", comment.char = "", check.names = FALSE)
  }
  read_matrix <- function(name) {
    con <- gzfile(file.path(path, name), "rt")
    on.exit(try(close(con), silent = TRUE))
    methods::as(Matrix::readMM(con), "CsparseMatrix")
  }
  features <- read_tab("features.tsv.gz")
  reference <- read_tab("reference_features.tsv.gz")
  barcodes <- read_tab("barcodes.tsv.gz")[[1]]
  if (ncol(features) < 2L || ncol(reference) < 2L ||
      anyDuplicated(features[[1]]) || anyDuplicated(reference[[1]]) ||
      anyDuplicated(barcodes)) stop("Invalid or duplicate exported identifiers.")
  # Read10X makes names unique before Seurat replaces underscores and pipes.
  reference_names <- gsub("[_|]", "-", make.unique(reference[[2]]))
  if (anyDuplicated(reference_names)) stop("Ambiguous sanitized feature names.")
  layers <- SeuratObject::Layers(seu[[source_assay]], search = "^counts$")
  if (length(layers) != 1L) stop("Source assay must have one unsplit counts layer.")
  original <- SeuratObject::LayerData(seu, assay = source_assay, layer = "counts")
  if (!identical(colnames(original), colnames(seu))) {
    stop("Source assay must contain all object cells in object order.")
  }
  ref_index <- match(rownames(original), reference_names)
  row_index <- match(reference[[1]][ref_index], features[[1]])
  col_index <- match(colnames(original), barcodes)
  if (anyNA(ref_index) || anyNA(row_index) || anyNA(col_index)) {
    stop("Some Seurat genes/cells cannot be matched. Check sample and feature naming.")
  }
  expected_dim <- c(nrow(features), length(barcodes))
  baseline <- read_matrix("original_counts.mtx.gz")
  if (!identical(dim(baseline), expected_dim)) stop("Original dimensions disagree with identifiers.")
  baseline <- baseline[row_index, col_index, drop = FALSE]
  dimnames(baseline) <- dimnames(original)
  delta <- Matrix::drop0(baseline - original)
  if (any(!is.finite(baseline@x)) || any(baseline@x < 0) ||
      (length(delta@x) > 0L && any(abs(delta@x) > 1e-8))) {
    stop("Exported original counts do not match the Seurat source assay.")
  }
  rm(delta, baseline)
  corrected <- read_matrix("corrected_counts.mtx.gz")
  if (!identical(dim(corrected), expected_dim)) stop("Corrected dimensions disagree with identifiers.")
  if (any(!is.finite(corrected@x)) || any(corrected@x < 0)) stop("Invalid corrected counts.")
  corrected <- corrected[row_index, col_index, drop = FALSE]
  dimnames(corrected) <- dimnames(original)
  delta <- Matrix::drop0(corrected - original)
  if (length(delta@x) > 0L && any(delta@x > 1e-8)) stop("Corrected counts exceed original counts.")
  rm(delta)
  gene_diag <- read_tab("gene_diagnostics.tsv.gz", header = TRUE)
  if (!all(c("gene_id", "alpha", "r2", "correction_eligible") %in% names(gene_diag)) ||
      anyDuplicated(gene_diag$gene_id)) stop("Invalid gene diagnostics.")
  gene_index <- match(features[[1]][row_index], gene_diag$gene_id)
  if (anyNA(gene_index)) stop("Missing gene diagnostics.")
  feature_meta <- gene_diag[gene_index, , drop = FALSE]
  rownames(feature_meta) <- rownames(original)
  before_gene <- Matrix::rowSums(original)
  feature_meta$fraction_removed <- ifelse(before_gene > 0,
    1 - Matrix::rowSums(corrected) / before_gene, NA_real_)
  new_assay <- SeuratObject::CreateAssay5Object(counts = corrected)
  new_assay <- SeuratObject::AddMetaData(new_assay, metadata = feature_meta)
  seu[[assay]] <- new_assay
  before_cell <- Matrix::colSums(original)
  after_cell <- Matrix::colSums(corrected)
  cell_meta <- data.frame(original_count = before_cell, corrected_count = after_cell,
                         counts_removed = before_cell - after_cell,
                         fraction_removed = ifelse(before_cell > 0,
                           1 - after_cell / before_cell, NA_real_),
                         row.names = colnames(original))
  names(cell_meta) <- paste0(assay, "_", names(cell_meta))
  seu <- SeuratObject::AddMetaData(seu, metadata = cell_meta)
  SeuratObject::DefaultAssay(seu) <- assay
  if (normalize) {
    seu <- Seurat::NormalizeData(seu, assay = assay, normalization.method = "LogNormalize",
                                 scale.factor = scale_factor, verbose = FALSE)
  }
  if (is.null(seu@misc$sparkle)) seu@misc$sparkle <- list()
  seu@misc$sparkle[[assay]] <- list(
    diagnostics = diagnostics, source_assay = source_assay,
    export_path = normalizePath(path), imported_at = as.character(Sys.time()),
    normalized = normalize, scale_factor = if (normalize) scale_factor else NULL,
    note = "Existing reductions, clusters and annotations describe the original data.")
  seu
}
