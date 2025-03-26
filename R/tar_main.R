targets_main <- function() {

  setup_experiment <- tar_plan(
    tar_target(config_file, "config.json", format = "file"),
    config = jsonlite::read_json(config_file)
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
    mito_genes =  gns |> filter(chr == "MT") |> pull(gene_symbol)
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
      seu_clust = cluster_seurat(seu, ident_res = RESOLUTION),
      cluster_markers = find_cluster_markers(seu_clust),

      # Gene detection
      gene_detection = get_gene_detection(seu),
      cell_data = get_cell_data(seu),
      tbl_cuts = make_table_cuts(params),


      # Figures
      fig_marginals = plot_cell_marginals(cell_data, params, with_vlines = FALSE),
      fig_count_det = plot_count_det(cell_data, params, with_lines = FALSE)
    )
  )

  merged_data <- tar_plan(
    cell_data_all = bind_rows(cell_data_1, cell_data_2, cell_data_3, cell_data_4, .id = "sample"),
    cell_count_all = cell_data_all |> count(sample),
    read_count_all = cell_data_all |> summarise(n = sum(total_count), .by = "sample")
  )

  for_report <- tar_plan(
    fig_counts = plot_all_counts(cell_data_all)
  )


  c(
    setup_experiment,
    get_annotations,
    map_experiments,
    merged_data,
    for_report
  )

}