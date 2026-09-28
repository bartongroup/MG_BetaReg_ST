# Not used. Tool designed for Visium, creates a dense distance matrix. The raw data
# contain 312506 spots (8um), so the distance matrix wouls be 312506 x 312506, which
# would take up 312506^2 * 8 bytes = 781 GB of memory. And the tool creats several 
# of these matrices, so it is not possible to run even on the larget cluster node.
#
# Relevant code in R/decontamination.R:
#
#  # Calculate Euclidean distances
#    d_mat <- dist(x, method = "euclidean", diag = TRUE)
#    d_mat <- as.matrix(d_mat)

load_spot_clean <- function(data_path) {
  matrix_path <- file.path(data_path, "binned_outputs/square_008um/raw_feature_bc_matrix")
  spatial_path <- file.path(data_path, "binned_outputs/square_008um/spatial")

  tissue_img_file <- file.path(spatial_path, "tissue_lowres_image.png")
  scale_factor_file <- file.path(spatial_path, "scalefactors_json.json")

  raw <- SpotClean::read10xRaw(matrix_path)

  # SpotClean cannot read Parquet files, so we need to convert the spatial data to a CSV file first
  positions <- arrow::read_parquet(file.path(spatial_path, "tissue_positions.parquet"))
  tissue_csv_file <- tempfile(fileext = ".csv")
  write_csv(positions, tissue_csv_file)
  
  slide_info <- SpotClean::read10xSlide(
    tissue_csv_file = tissue_csv_file,
    tissue_img_file = tissue_img_file,
    scale_factor_file = scale_factor_file
  )

  SpotClean::createSlide(count_mat = raw, slide_info = slide_info)
}