targets_2um <- function() {
  map_experiments <- tar_map(
    values = EXPERIMENTS,
    names = NAME2,
    
    tar_plan(
      
      # Create Seurat object
      seur = load_visium_data(DATA_PATH, mito_genes, bin_size = 2),
      seu = process_seurat(seur, nfeatures = 2000),
      seu_mark = annotate_cells_page(seu, tissue_markers, n_boot = 100, p_val_limit = 0.05, n_cores = 1, verbose = TRUE),
      #seu_clust = seur |> NormalizeData() |> cluster_seurat(ident_res = RESOLUTION, resolutions = RESOLUTIONS, nfeatures = 2000),
  
      # General figures
      fig_marginals = plot_cell_marginals(cell_data |> dplyr::select(total_count, n_genes, mito_fraction), params, with_vlines = FALSE),
      fig_count_det = plot_count_det(cell_data |> dplyr::select(total_count, n_genes, mito_fraction), params, with_lines = FALSE),
  
      # Gene detection
      gene_detection = get_gene_detection(seu) |> tibble::add_column(sample = NAME2),
      cell_data = get_cell_data(seu)|> tibble::add_column(sample = NAME2),
      tbl_hormone_overlap = combat_upset_expressed_cells(seu, c("INS", "SST", "GCG", "GHRL", "PPY"), "hormone_overlap", NAME2),
      tbl_hormone_count = get_hormone_count(seu) |> tb("hormone_count", NAME2),
      tbl_mpc_overlap = combat_upset_expressed_cells(seu, c("SOX9", "PTF1A", "NKX6-1", "PDX1"), "mpc_overlap", NAME2),
      tbl_ep_overlap = combat_upset_expressed_cells(seu, c("PAX4", "PAX6", "NEUROD1", "NEUROG3"), "ep_overlap", NAME2),

      key_marker_genes = key_markers$gene_symbol |> unique() |> intersect(SeuratObject::Features(seu)),
      key_marker_genes_endoc = key_markers_endoc$gene_symbol |> levels() |> intersect(SeuratObject::Features(seu)),
  
      # Clustering
      #cluster_markers = find_group_markers(seu, group_by = "cluster"),
      #top_cluster_markers = get_top_cluster_markers(cluster_markers, n_top = 12),
      #top_cluster_genes = top_cluster_markers$gene_symbol |> unique(),
      #tbl_cluster_enrichment = group_enrichment(top_cluster_markers, fterms, group_var = "cluster"),
      #png_resolutions = plot_resolutions(seu, RESOLUTIONS) |> gs("resolutions", NAME2, 8, 32),
      #png_top_cluster = plot_top_cluster_umap_spatial(seu, top_cluster_markers, file_prefix = paste0("top_clust_", NAME2)),
      #png_top_cluster_expr = plot_gene_expression(seu, top_cluster_genes) |> gs("top_genes", NAME2, 8, 12),
      #fig_dots_top_cluster = plot_gene_dots(seu, top_cluster_genes),
  
      # Key markers
      png_key_markers = plot_key_markers_umap_spatial(seu, key_markers, file_prefix = paste0("key_markers_", NAME2), width = 32),
      #png_key_markers_expr = plot_gene_expression(seu, key_marker_genes) |> gs("key_marker_expr", NAME2, 8, 12),
      png_dots_key_markers = plot_gene_dots(seu, key_marker_genes)  |> gs("dots_key_markers", NAME2, 8, 5),
      #png_xy_marker_mpc_ep = plot_gene_vs_gene(seu, goi$tissue_markers_positive$MPC, goi$tissue_markers_positive$EP, group_sel = c(2, 6)) |> gs("xy_marker_mpc_ep", NAME2, 6, 5),
      #png_xy_marker_mpc_beta = plot_gene_vs_gene(seu, goi$tissue_markers_positive$MPC, goi$tissue_markers_positive$Beta, group_sel = c(2, 6)) |> gs("xy_marker_mpc_beta", NAME2, 6, 5),
  
      # Cell types
      png_cell_type_scores = plot_cell_type_scores(seu_mark, tissue_markers) |> gs("cell_type_scores", NAME2, 8, 22),
      png_cell_type = plot_cell_type(seu_mark) |> gs("cell_type", NAME2, 32, 16),
      #fig_cluster_cell_type = plot_cluster_cell_type_jaccard(seu_mark),
      png_dots_cell_type = plot_gene_dots(seu_mark, unique(tissue_markers$gene_symbol), meta_column = "cell_type") |> gs("dots_cell_type", NAME2, 6, 6),
      tbl_cell_type_count = cell_type_count(seu_mark) |> tb("cell_type_count", NAME2),
      tbl_expr_pct = get_expression_percentage(seu_mark, unique(tissue_markers$gene_symbol), meta_column = "cell_type") |> tb("expr_pct", NAME2),
      png_dif_spatial = plot_spatial_differences(seu_mark, pairs = list(c("EP", "MPC"), c("Beta", "MPC"), c("Beta", "Alpha"), c("Beta", "Delta"), c("Beta", "Epsilon")), point_size = 0.85) |> gs("dif_spatial", NAME2, 8, 30),
      cell_type_markers = find_group_markers(seu_mark, group_by = "cell_type"),
      png_cell_type_marker_dots = plot_cell_type_marker_dots(seu_mark, cell_type_markers) |> gs("cell_type_marker_dots", NAME2, 8, 35),
      cell_type_markers_sel = cell_type_markers |> dplyr::filter(p_val_adj < 0.01 & avg_log2FC > 1),
      #cell_type_enrichment = group_enrichment(cell_type_markers_sel, fterms, group_var = "cluster"),
      png_gene_gene_ins_sst = plot_gene_vs_gene(seu_mark, "INS", "SST", group_var = "cell_type", group_sel = c("Beta", "Delta")) |> gs("gene_gene_ins_sst", NAME2, width = 6, height = 5),
      png_gene_gene_ins_gcg = plot_gene_vs_gene(seu_mark, "INS", "GCG", group_var = "cell_type", group_sel = c("Beta", "Alpha")) |> gs("gene_gene_ins_gcg", NAME2, width = 6, height = 5),
      png_gene_gene_ins_ghrl = plot_gene_vs_gene(seu_mark, "INS", "GHRL", group_var = "cell_type", group_sel = c("Beta", "Epsilon")) |> gs("gene_gene_ins_ghrl", NAME2, width = 6, height = 5),
  
      # Insulin pathway
      png_ins_path = plot_expression(seu, goi$insulin_signalling_pathway, umap_point_size = 0.5, spatial_point_size = 0.5) |> gs("insulin_path", NAME2, 7, 48),
      png_ins_path_expr = plot_gene_expression(seu, goi$insulin_signalling_pathway, meta_column = "cell_type") |> gs("insulin_path_expr", NAME2, 8, 12),
      fig_dots_ins_path = plot_gene_dots(seu, goi$insulin_signalling_pathway, meta_column = "cell_type"),
  
    )
  )

  combine_samples <- list(
    tar_combine(cell_data_all_2um, map_experiments$cell_data, command = bind_rows(!!!.x)),
    tar_combine(gene_detection_all_2um, map_experiments$gene_detection, command = bind_rows(!!!.x))
    #tar_combine(gse_all, map_experiments$gse_mpc, command = bind_rows(!!!.x))
  )

  for_report <- tar_plan(
    png_counts_2um = plot_all_counts(cell_data_all_2um, gene_detection_all_2um) |> gs("counts_2um", "all", 8, 3),
    png_marginals_2um = plot_grid(fig_marginals_2um_1, fig_marginals_2um_2, ncol = 1, labels = 1:2)  |> gs("marginals_2um", "all", 8, 6),
    png_detections_2um = plot_grid(fig_count_det_2um_1, fig_count_det_2um_2, labels = 1:2) |> gs("detections_2um", "all", 8, 3),

    sct_2um_1 = tibble::tibble(x = c(5120, 5230), y = c(-3400, -5400)),
    png_section_ins_sst_2um = plot_spatial_section(seu_2um_1, "INS", "SST", sct_2um_1, slice = "slice1.002um") |> gs("spatial_section_ins_sst_2um", width = 32, height = 32)
    #png_section_ins_gcg = plot_spatial_section(seu_1, "INS", "GCG", sct_1) |> gs("spatial_section_ins_gcg", width = 8, height = 8),
    #png_section_ins_ghrl = plot_spatial_section(seu_1, "INS", "GHRL", sct_1) |> gs("spatial_section_ins_ghrl", width = 8, height = 8)
    
  )


  c(
    map_experiments,
    combine_samples,
    for_report
  )
}