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

    terms = get_functional_terms_sym(config_annot$designations),
    fterms = prepare_terms_fenr(terms, SeuratObject::Features(seu_1), feature_name = "gene_symbol"),

    cell_cycle_markers = get_cell_cycle_markers(config$species, gns),
    mito_genes =  gns |> filter(chr == "MT") |> pull(gene_symbol),

    key_markers = read_markers("info/key_markers.xlsx"),

    mpc_markers = c("PTF1A", "SOX9", "NKX6-1", "PDX1"),
    ep_markers = c("NEUROD1", "NEUROG3", "NKX2-2", "PAX4", "PAX6"),
    beta_markers = "INS",
    nervous_markers = c("DBH", "PRPH", "GAL", "GAP43"),
    vsm_markers = c("MYH11", "TAGLN", "ADIRF"),
    immune_markers = c("LYZ", "S100A9"),
    matrix_markers = c("COL16A1", "LUM", "POSTN"),

    select_markers = get_unique_markers(mpc_markers, ep_markers, beta_markers, nervous_markers, vsm_markers, immune_markers, matrix_markers)
  )

  map_experiments <- tar_map(
    values = EXPERIMENTS,
    names = NAME,
    
    tar_plan(
      
      # Create Seurat object
      seu = load_visium_data(DATA_PATH, mito_genes) |> NormalizeData() |> cluster_seurat(ident_res = RESOLUTION, resolutions = RESOLUTIONS),
      seu_mark = seu |> mark_cell_types(select_markers),
      cluster_markers = find_cluster_markers(seu),
      top_cluster_markers = get_top_cluster_markers(cluster_markers, n_top = 12),
      key_marker_genes = key_markers$gene_symbol |> unique() |> intersect(SeuratObject::Features(seu)),
      top_cluster_genes = top_cluster_markers$gene_symbol |> unique(),

      # Gene detection
      gene_detection = get_gene_detection(seu) |> tibble::add_column(sample = NAME),
      cell_data = get_cell_data(seu)|> tibble::add_column(sample = NAME),

      # Genes
      tbl_key_corr = find_gene_correlation(seu, key_marker_genes),
      top_key_corr = tbl_key_corr |> dplyr::filter(correlation > 0.3),
      tbl_cluster_enrichment = group_enrichment(top_cluster_markers, fterms, group_var = "cluster"),
      tbl_keycorr_enrichment = group_enrichment(top_key_corr, fterms, group_var = "sel_gene"),

      # Figures
      fig_marginals = plot_cell_marginals(cell_data, params, with_vlines = FALSE),
      fig_count_det = plot_count_det(cell_data, params, with_lines = FALSE),

      # Fig clustering
      png_resolutions = plot_resolutions(seu, RESOLUTIONS) |> gs("resolutions", NAME, 8, 32),
      png_top_cluster = plot_top_cluster_umap_spatial(seu, top_cluster_markers, file_prefix = paste0("top_clust_", NAME)),
      png_top_cluster_expr = plot_gene_expression(seu, top_cluster_genes) |> gs("top_genes", NAME, 8, 12),
      fig_dots_top_cluster = plot_gene_dots(seu_1, top_cluster_genes),

      # Fig key markers
      png_key_markers = plot_key_markers_umap_spatial(seu, key_markers, file_prefix = paste0("key_markers_", NAME)),
      png_key_markers_expr = plot_gene_expression(seu, key_marker_genes) |> gs("key_marker_expr", NAME, 8, 12),
      fig_dots_key_markers = plot_gene_dots(seu, key_marker_genes),
      fig_cor_heat_key_markers = plot_gene_cor(tbl_key_corr, min_cor = 0.3),
      png_xy_marker_mpc_ep = plot_gene_vs_gene_cluster(seu, mpc_markers, ep_markers, clust_sel = c(2, 6)) |> gs("xy_marker_mpc_ep", NAME, 6, 5),
      png_xy_marker_mpc_beta = plot_gene_vs_gene_cluster(seu, mpc_markers, beta_markers, clust_sel = c(2, 6)) |> gs("xy_marker_mpc_beta", NAME, 6, 5),

      # Cell types
      fig_cell_type = plot_cell_type(seu_mark) |> gs("cell_type", NAME, 8, 4)
    )
  )

  combine_samples <- list(
    tar_combine(cell_data_all, map_experiments$cell_data, command = bind_rows(!!!.x)),
    tar_combine(gene_detection_all, map_experiments$gene_detection, command = bind_rows(!!!.x))
  )

  for_report <- tar_plan(
    fig_counts = plot_all_counts(cell_data_all, gene_detection_all),
    png_marginals = plot_grid(fig_marginals_1, fig_marginals_2, ncol = 1, labels = 1:2)  |> gs("marginals", "all", 8, 6),
    png_detections = plot_grid(fig_count_det_1, fig_count_det_2, labels = 1:2) |> gs("detections", "all", 8, 3)
  )

  c(
    setup_experiment,
    get_annotations,
    map_experiments,
    combine_samples,
    for_report
  )

}