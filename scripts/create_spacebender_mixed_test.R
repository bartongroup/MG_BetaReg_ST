# Run from the project root in its R environment:
# gtimeout 120s Rscript --vanilla scripts/create_spacebender_mixed_test.R
# This is an edge-case test, not a representative biological subset.
suppressPackageStartupMessages({
  library(Matrix)
  library(hdf5r)
  library(arrow)
  library(anndataR)
})

bin_path <- "data/spaceranger_output/1_18191_Results/outs/binned_outputs/square_008um"
output <- "spacebender_test/seu_1_mixed_1000.h5ad"
if (file.exists(output)) stop("Output already exists: ", output)

h5 <- H5File$new(file.path(bin_path, "raw_feature_bc_matrix.h5"), "r")
raw <- tryCatch({
  list(
    counts = sparseMatrix(
      i = as.integer(h5[["matrix/indices"]][]),
      p = as.integer(h5[["matrix/indptr"]][]),
      x = as.numeric(h5[["matrix/data"]][]), index1 = FALSE,
      dims = as.integer(h5[["matrix/shape"]][])
    ),
    barcodes = h5[["matrix/barcodes"]][],
    genes = h5[["matrix/features/name"]][],
    ids = h5[["matrix/features/id"]][],
    types = h5[["matrix/features/feature_type"]][]
  )
}, finally = h5$close_all())
gene_names <- gsub("[_|]", "-", make.unique(raw$genes))
stopifnot(!anyDuplicated(gene_names))
keep <- raw$types == "Gene Expression"
counts <- raw$counts[keep, , drop = FALSE]
totals <- colSums(counts)
pos <- as.data.frame(read_parquet(file.path(bin_path, "spatial/tissue_positions.parquet")))
idx <- match(raw$barcodes, pos$barcode)
stopifnot(!anyNA(idx), !anyDuplicated(pos$barcode))
pos <- pos[idx, , drop = FALSE]
stopifnot(all(pos$in_tissue %in% c(0, 1)))

# Choose a tissue-edge anchor, then the nearest bins in each category.
# Membership and counts come from Space Ranger, not simulated labels/zeroes.
grid_key <- paste(pos$array_row, pos$array_col, sep = ":")
boundary <- rep(FALSE, nrow(pos))
for (offset in list(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))) {
  neighbour <- match(paste(pos$array_row + offset[1], pos$array_col + offset[2], sep = ":"), grid_key)
  boundary <- boundary | (!is.na(neighbour) & pos$in_tissue[neighbour] == 0)
}
boundary[is.na(boundary)] <- FALSE
candidate <- which(boundary & pos$in_tissue == 1 & totals > 0)
if (!length(candidate)) stop("No tissue boundary found.")
inside <- pos$in_tissue == 1
centre <- c(median(pos$array_row[inside]), median(pos$array_col[inside]))
anchor <- candidate[which.min((pos$array_row[candidate] - centre[1])^2 +
                              (pos$array_col[candidate] - centre[2])^2)]
distance <- (pos$array_row - pos$array_row[anchor])^2 +
  (pos$array_col - pos$array_col[anchor])^2
nearest <- function(eligible, n) {
  indices <- which(eligible)
  if (length(indices) < n) stop("Not enough bins in one test category.")
  indices[order(distance[indices], indices)][seq_len(n)]
}
selected <- c(
  nearest(pos$in_tissue == 1 & totals > 0, 500),
  nearest(pos$in_tissue == 0 & totals > 0, 300),
  nearest(pos$in_tissue == 0 & totals == 0, 200)
)
# Empty rows deliberately come last to create entirely empty posterior batches.
# With batch_size=64, rows 833:896 and 897:960 are full empty batches.
groups <- rep(c("tissue_nonempty", "outside_nonempty", "outside_empty"), c(500, 300, 200))
spatial <- as.matrix(pos[selected, c("pxl_col_in_fullres", "pxl_row_in_fullres")])
dimnames(spatial) <- NULL
storage.mode(spatial) <- "double"
adata <- AnnData(
  X = t(counts[, selected, drop = FALSE]),
  obs = data.frame(
    in_tissue = as.logical(pos$in_tissue[selected]), test_group = groups,
    total_counts = unname(totals[selected]),
    array_row = pos$array_row[selected], array_col = pos$array_col[selected],
    row.names = raw$barcodes[selected]
  ),
  var = data.frame(gene_ids = raw$ids[keep], gene_name = raw$genes[keep],
                   feature_types = raw$types[keep], row.names = gene_names[keep]),
  obsm = list(spatial = spatial),
  uns = list(test_purpose = "Empty posterior batch stress test; not a biological benchmark",
             spatial_coordinate_units = "full-resolution image pixels", bin_size_um = 8)
)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
adata$write_h5ad(output, compression = "gzip", mode = "w-")

check <- read_h5ad(output)
stopifnot(nrow(check$obs) == 1000L, inherits(check$X, "sparseMatrix"),
          identical(check$obs_names, raw$barcodes[selected]),
          all(rowSums(check$X) == totals[selected]),
          all(rowSums(check$X)[801:1000] == 0),
          all(as.matrix(check$obsm[["spatial"]]) == spatial),
          all(check$X == t(counts[, selected, drop = FALSE])),
          identical(check$obs$in_tissue, as.logical(pos$in_tissue[selected])))
print(table(check$obs$test_group, check$obs$in_tissue))
message("Verified ", nrow(check$obs), " bins x ", nrow(check$var), " genes: ", output)
