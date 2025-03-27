targets_main <- function() {

  setup_experiment <- tar_plan(
    tar_target(config_file, "config.json", format = "file"),
    config = jsonlite::read_json(config_file),
    RESOLUTIONS = c(0.05, 0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.3, 0.4)
  )

  get_annotations <- tar_plan(
    # isolate config to avoid re-runs in case of other config changes
    config_annot = list(
      species = config$species,
      ensembl = config$ensembl,
      designations = config$designations
    ), 
    mart = new_mart(config_annot),
    chromosomes = biomart_fetch_chromosomes(mart),
    gns =  biomart_fetch_genes(mart) |> filter(chr %in% chromosomes),

    cell_cycle_markers = get_cell_cycle_markers(config$species, gns),
    mito_genes =  gns |> filter(chr == "MT") |> pull(gene_symbol),

    key_markers = read_markers("info/Key markers of pancreatic development stages.xlsx")
  )

  map_experiments <- tar_map(
    values = EXPERIMENTS,
    names = NAME,
    
    tar_plan(
      
      # Wrap all seurat parameters in one list
      params = list(
        n_variable_features = N_VARIABLE_FEATURES,
        min_detected_genes_per_cell = MIN_DETECTED_GENES_PER_CELL,
        min_counts_per_cell = MIN_COUNTS_PER_CELL,
        min_detections_per_gene = MIN_DETECTIONS_PER_GENE,
        max_mt_fraction = MAX_MT_FRACTION
      ),
      
      # Create Seurat object
      seu = load_visium_data(DATA_PATH, mito_genes) |> NormalizeData(),
      seu_clust = cluster_seurat(seu, ident_res = RESOLUTION, resolutions = RESOLUTIONS),
      cluster_markers = find_cluster_markers(seu_clust),
      top_cluster_markers = get_top_cluster_markers(cluster_markers, n_top = 12),

      # Gene detection
      gene_detection = get_gene_detection(seu) |> tibble::add_column(sample = NAME),
      cell_data = get_cell_data(seu)|> tibble::add_column(sample = NAME),
      tbl_cuts = make_table_cuts(params),

      # Figures
      fig_marginals = plot_cell_marginals(cell_data, params, with_vlines = FALSE),
      fig_count_det = plot_count_det(cell_data, params, with_lines = FALSE),
      png_resolutions = plot_resolutions(seu_clust, RESOLUTIONS) |> gs("resolutions", NAME, 8, 32),
      png_top_cluster_umap = plot_top_cluster_features(seu_clust, top_cluster_markers, what = "umap", file_prefix = paste0("top_clust_umap_", NAME)),
      png_top_cluster_spatial = plot_top_cluster_features(seu_clust, top_cluster_markers, what = "spatial", file_prefix = paste0("top_clust_spatial_", NAME)),
      png_top_cluster_expr = plot_expression_clusters(seu_clust, top_cluster_markers$gene_symbol |> unique(), RESOLUTION) |> gs("top_genes", NAME, 8, 12),
      #png_key_markers = plot_key_markers_umap_spatial(seu, key_markers) |> gs("key_markers", NAME, 8, 30),

      fig_dots_cluster = plot_dots(seu_clust, top_cluster_markers, resolution = RESOLUTION)

    )
  )

  combine_samples <- list(
    tar_combine(cell_data_all, map_experiments$cell_data, command = bind_rows(!!!.x)),
    tar_combine(gene_detection_all, map_experiments$gene_detection, command = bind_rows(!!!.x))
  )

  for_report <- tar_plan(
    fig_counts = plot_all_counts(cell_data_all, gene_detection_all)
  )

  c(
    setup_experiment,
    get_annotations,
    map_experiments,
    combine_samples,
    for_report
  )

}