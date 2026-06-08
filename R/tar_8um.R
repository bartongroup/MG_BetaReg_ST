targets_8um <- function() {

  map_experiments <- tar_map(
    values = EXPERIMENTS,
    names = NAME,
    
    tar_plan(
      
      # Create Seurat object
      seu = load_visium_data(DATA_PATH, mito_genes) |>
        NormalizeData() |>
        cluster_seurat(ident_res = RESOLUTION, resolutions = RESOLUTIONS, nfeatures = 2000) |> 
        relabel_clusters(RELABEL_CLUSTERS, RESOLUTION),
      seu_mark = annotate_cells_page(seu, tissue_markers, n_boot = 1000, p_val_limit = 0.05),

      # Save seurat for shiny
      sav_seu = save_for_shiny(seu_mark, NAME),

      # General figures
      fig_marginals = plot_cell_marginals(cell_data |> dplyr::select(total_count, n_genes, mito_fraction), params, with_vlines = FALSE),
      fig_count_det = plot_count_det(cell_data |> dplyr::select(total_count, n_genes, mito_fraction), params, with_lines = FALSE),

      # Gene detection
      gene_count = dim(seu)[1],
      gene_detection = get_gene_detection(seu) |> tibble::add_column(sample = NAME),
      cell_data = get_cell_data(seu)|> tibble::add_column(sample = NAME),
      fig_hormone_overlap = plot_upset_expressed_cells(seu, c("INS", "SST", "GCG", "GHRL", "PPY")),
      hormone_count = get_hormone_count(seu),
      fig_mpc_overlap = plot_upset_expressed_cells(seu, c("SOX9", "PTF1A", "NKX6-1", "PDX1")),
      fig_ep_overlap = plot_upset_expressed_cells(seu, c("PAX4", "PAX6", "NEUROD1", "NEUROG3")),
      sry_expression_count = get_gene_expression(seu, "SRY") |> summarise(n = length(which(expression > 0))),

      # Clustering
      cluster_markers = find_group_markers(seu, group_by = "cluster"),
      max_best_clusters = 30,
      tbl_best_cluster_markers = find_best_markers(cluster_markers, max_per_cluster = max_best_clusters),
      tbl_best_cluster_negative_markers = find_best_markers(cluster_markers, max_per_cluster = max_best_clusters, negative = TRUE),
      top_cluster_markers = get_top_cluster_markers(cluster_markers, n_top = 12),
      top_cluster_genes = top_cluster_markers$gene_symbol |> unique(),
      tbl_cluster_enrichment = group_enrichment(top_cluster_markers, fterms, group_var = "cluster"),
      png_resolutions = plot_resolutions(seu, RESOLUTIONS) |> gs("resolutions", NAME, 8, 32),
      png_top_cluster = plot_top_cluster_umap_spatial(seu, top_cluster_markers, file_prefix = paste0("top_clust_", NAME)),
      png_top_cluster_expr = plot_gene_expression(seu, top_cluster_genes) |> gs("top_genes", NAME, 8, 12),
      fig_dots_top_cluster = plot_gene_dots(seu, top_cluster_genes),

      # Key markers
      png_key_markers = plot_key_markers_umap_spatial(seu, key_markers, file_prefix = paste0("key_markers_", NAME)),
      png_key_markers_expr = plot_gene_expression(seu, key_marker_genes) |> gs("key_marker_expr", NAME, 8, 12),
      fig_dots_key_markers = plot_gene_dots(seu, key_marker_genes),
      png_xy_marker_mpc_ep = plot_gene_vs_gene(seu, goi$tissue_markers_positive$MPC, goi$tissue_markers_positive$EP, group_sel = c(2, 6)) |> gs("xy_marker_mpc_ep", NAME, 6, 5),
      png_xy_marker_mpc_beta = plot_gene_vs_gene(seu, goi$tissue_markers_positive$MPC, goi$tissue_markers_positive$Beta, group_sel = c(2, 6)) |> gs("xy_marker_mpc_beta", NAME, 6, 5),

      # Insulin gradient
      n_insulin_bins = 6,
      max_best_insulin = 100,
      tbls_binned_enrichment = enrichment_gradient(seu, fterms, "INS", bins = n_insulin_bins),
      tbl_best_insulin_gradient_markers = find_best_markers(tbls_binned_enrichment$markers, max_per_cluster = max_best_insulin),
      tbl_best_insulin_gradient_negative_markers = find_best_markers(tbls_binned_enrichment$markers, max_per_cluster = max_best_insulin, negative = TRUE),
      fig_enrichment_gradient = plot_enrichment_gradient(tbls_binned_enrichment$enrichment),
      genes_for_binned = c("PTF1A", "PDX1", "NKX6-1", "SOX9", "PAX4", "PAX6", "NKX2-2", "NEUROD1", "NEUROG3", "GCG", "SST", "GHRL"),
      fig_binned_expression = plot_binned_expression(seu, "INS", genes_for_binned, bins = n_insulin_bins),
      fig_binned_ins_prop = plot_binned_xy_expression(seu, "INS", genes_for_binned, bins = n_insulin_bins, what = "proportion", include_x_zeroes = TRUE),
      fig_binned_ins_spatial = plot_binned_spatial(seu, "INS", bins = n_insulin_bins) |> gs("binned_ins_spatial", NAME, 8, 6),

      # Cell types
      png_cell_type_scores = plot_cell_type_scores(seu_mark, tissue_markers) |> gs("cell_type_scores", NAME, 8, 22),
      png_cell_type = plot_cell_type(seu_mark) |> gs("cell_type", NAME, 8, 4),
      fig_cluster_cell_type = plot_cluster_cell_type_jaccard(seu_mark),
      fig_dots_cell_type = plot_gene_dots(seu_mark, unique(tissue_markers$gene_symbol), meta_column = "cell_type"),
      tbl_cell_type_count = cell_type_count(seu_mark),
      tbl_expr_pct = get_expression_percentage(seu_mark, unique(tissue_markers$gene_symbol), meta_column = "cell_type"),
      png_dif_spatial = plot_spatial_differences(seu_mark, pairs = list(c("EP", "MPC"), c("Beta", "MPC"), c("Beta", "EP"), c("Beta", "Alpha"), c("Beta", "Delta"), c("Beta", "Epsilon"))) |> gs("dif_spatial", NAME, 8, 36),
      cell_type_markers = find_group_markers(seu_mark, group_by = "cell_type"),
      png_cell_type_marker_dots = plot_cell_type_marker_dots(seu_mark, cell_type_markers) |> gs("cell_type_marker_dots", NAME, 8, 35),
      cell_type_markers_sel = cell_type_markers |> dplyr::filter(p_val_adj < 0.01 & avg_log2FC > 1),
      cell_type_enrichment = group_enrichment(cell_type_markers_sel, fterms, group_var = "cluster"),
      png_gene_gene_ins_sst = plot_gene_vs_gene(seu_mark, "INS", "SST", group_var = "cell_type", group_sel = c("Beta", "Delta")) |> gs("gene_gene_ins_sst", NAME, width = 6, height = 5),
      png_gene_gene_ins_gcg = plot_gene_vs_gene(seu_mark, "INS", "GCG", group_var = "cell_type", group_sel = c("Beta", "Alpha")) |> gs("gene_gene_ins_gcg", NAME, width = 6, height = 5),
      png_gene_gene_ins_ghrl = plot_gene_vs_gene(seu_mark, "INS", "GHRL", group_var = "cell_type", group_sel = c("Beta", "Epsilon")) |> gs("gene_gene_ins_ghrl", NAME, width = 6, height = 5),
      fig_sankey_cell_type = plot_sankey_cell_type(seu_mark),

      # Unknown
      fig_cell_type_cluster_overlap = plot_cluster_cell_type_overlap(seu_mark),
      seu_unknown = seu_mark |> subset(cell_type == "Unknown") |>
        NormalizeData() |>
        cluster_seurat(ident_res = RESOLUTION_UNKNOWN, resolutions = RESOLUTIONS_UNKNOWN, nfeatures = 2000),
      cluster_markers_unknown = find_group_markers(seu_unknown, group_by = "cluster"),
      top_cluster_markers_unknown = get_top_cluster_markers(cluster_markers_unknown, n_top = 12),
      top_cluster_genes_unknown = top_cluster_markers_unknown$gene_symbol |> unique(),
      tbl_cluster_enrichment_unknown = group_enrichment(top_cluster_markers_unknown, fterms, group_var = "cluster"),
      png_resolutions_unknown = plot_resolutions(seu_unknown, RESOLUTIONS_UNKNOWN) |> gs("resolutions_unknown", NAME, 8, 32),
      png_top_cluster_unknown = plot_top_cluster_umap_spatial(seu_unknown, top_cluster_markers_unknown, file_prefix = paste0("top_clust_unknown_", NAME)),
      png_top_cluster_expr_unknown = plot_gene_expression(seu_unknown, top_cluster_genes_unknown) |> gs("top_genes_unknown", NAME, 8, 12),
      fig_dots_top_cluster_unknown = plot_gene_dots(seu_unknown, top_cluster_genes_unknown),


      # Insulin pathway
      png_ins_path = plot_expression(seu, goi$insulin_signalling_pathway, umap_point_size = 0.5) |> gs("insulin_path", NAME, 7, 48),
      png_ins_path_expr = plot_gene_expression(seu, goi$insulin_signalling_pathway, meta_column = "cell_type") |> gs("insulin_path_expr", NAME, 8, 12),
      fig_dots_ins_path = plot_gene_dots(seu, goi$insulin_signalling_pathway, meta_column = "cell_type"),

      # Risk genes
      fig_dots_t1d = plot_gene_dots(seu_mark, risk_genes$t1d, meta_column = "cell_type"),
      fig_dots_t2d = plot_gene_dots(seu_mark, risk_genes$t2d, meta_column = "cell_type"),

      # Correlation
      key_marker_genes = key_markers$gene_symbol |> unique() |> intersect(SeuratObject::Features(seu)),
      key_marker_genes_endoc = key_markers_endoc$gene_symbol |> levels() |> intersect(SeuratObject::Features(seu)),
      tbl_key_corr_spear = find_gene_correlation_cor(seu, key_marker_genes_endoc, method = "spearman", n_top = 3000),
      top_key_corr_spear = tbl_key_corr_spear |> dplyr::filter(correlation > 0.3),
      
      # The next command requires > 36 GB memory, it crashes my laptop; need to run scripts/cscore_correlation_all.R
      # on the cluster and copy back results. Input data for the script are generated by save_for_external_correlation().
      # tbl_key_corr_cscore = find_gene_correlation_cscore(seu),
      
      # Next two lines are for running correlation on a large-memory cluster node.
      sav_external_cor = save_for_cluster_run(seu, NAME),  # prepare date for external script
      cor_file = paste0("ext/cor_", NAME, ".parquet"),              # load results after running the external script
      tbl_key_corr_cscore = read_correlation_markers(cor_file, key_marker_genes_endoc),

      tbl_key_corr_cscore_endoc = tbl_key_corr_cscore |> dplyr::filter(gene_2 %in% key_marker_genes_endoc),
      tbl_key_corr_cscore_others = tbl_key_corr_cscore |> dplyr::filter(!(gene_2 %in% key_marker_genes_endoc)),
      tbl_key_corr_cscore_ins_path = tbl_key_corr_cscore |> dplyr::filter(gene_2 %in% goi$insulin_signalling_pathway),
      tbl_top_cor = tabulate_correlated_genes(tbl_key_corr_cscore_others, n_top = 30, fdr_limit = 0.05),
      fig_cor_heat_key_markers = plot_cor_heat(tbl_key_corr_cscore_others, n_top = 5, fdr_limit = 0.05, cluster_rows = FALSE, cluster_columns = FALSE),
      fig_cor_heat_only_markers = plot_cor_heat(tbl_key_corr_cscore_endoc, fdr_limit = 0.05),
      fig_cor_heat_ins_path = plot_cor_heat(tbl_key_corr_cscore_ins_path, cluster_rows = FALSE, cluster_columns = FALSE),

      # Correlation enrichment
      top_key_corr_cscore = tbl_key_corr_cscore |> dplyr::filter(correlation > 0.6 & p_adj < 0.01),
      tbl_keycorr_enrichment = group_enrichment(top_key_corr_cscore, fterms, group_var = "gene_1", gene_var = "gene_2"),

      # Correlation WGCNA
      # wg = wgcna_analysis(cor_file),

      # hdWGCNA
      seu_hdwgcna = hd_wgcna_analysis(seu, NAME, soft_power = NULL, merge_cut_height = WGCNA_MERGE_CUT_HEIGHT, min_module_size = WGCNA_MIN_MODULE_SIZE),
      png_module_dendro = plot_hdwgcna_dendrogram(seu_hdwgcna, NAME),
      fig_metacells_spatial = plot_metacells_spatial(seu_hdwgcna),
      png_modules_spatial = plot_modules_spatial(seu_hdwgcna) |> gs("modules_spatial", NAME, 8, 13),
      png_modules_expr_score = plot_modules_expr_score(seu_hdwgcna, n_genes = "all") |> gs("modules_score", NAME, 8, 18),
      igr_hub_gene_network = plot_hub_gene_network(seu_hdwgcna, n_hubs = 12),
      tbl_module_score = get_modules_expr_score(seu_hdwgcna),
      tbl_modules = get_modules(seu_hdwgcna),
      tbl_hub_genes = get_hub_genes(seu_hdwgcna),
      tbl_module_enrichment = module_enrichment(tbl_modules, gene_detection$gene_symbol, fterms, kme_limit = 0, fdr_limit = 0.05),
      png_modules_dots = plot_modules_dots(seu_mark, tbl_hub_genes) |> gs("modules_dots", NAME, 8, 25),

      # CellChat
      cellchat = create_cellchat(seu_mark) |> preprocess_cellchat(),
      cellchat_comp = compute_cellchat(cellchat, n_boot = 100, scale_distance = CELLCHAT_SCALE_DISTANCE, raw_use = TRUE, population_size = FALSE),
      cellchat_prob = CellChat::subsetCommunication(cellchat_comp) |> tibble::as_tibble() |> dplyr::arrange(desc(prob)),
      # This gives huge amount of interactions with very low probabilities:
      # cellchat_comp_smth = compute_cellchat(cellchat, n_boot = 100, scale_distance = CELLCHAT_SCALE_DISTANCE, raw_use = FALSE),
      # This gives exactly the same interactions as the main one, but probabilities are very low
      # cellchat_comp_pop = compute_cellchat(cellchat, n_boot = 100, scale_distance = CELLCHAT_SCALE_DISTANCE, population_size = TRUE),
      png_cellchat_circle = plot_cellchat_circle(cellchat_comp, paste0("cellchat_circle_", NAME, ".png")),
      png_cellchat_heatmap_count = plot_cellchat_heatmap(cellchat_comp, paste0("cellchat_heatmap_count", NAME, ".png"), measure = "count"),
      png_cellchat_heatmap_weight = plot_cellchat_heatmap(cellchat_comp, paste0("cellchat_heatmap_weight", NAME, ".png"), measure = "weight"),
      fig_cellchat_rank = plot_cellchat_rank(cellchat_comp),
      fig_cellchat_bubble = plot_cellchat_bubble_sel(cellchat_comp, prob_limit = 0.01),
      png_cellchat_chords = plot_cellchat_chords(cellchat_comp, NAME, CELLCHAT_PATHWAYS),
      png_cellchat_networks = plot_cellchat_networks(cellchat_comp, NAME, CELLCHAT_PATHWAYS),

      # TWCOM
      # ABANDONED. Uses more than 500 GB of memory.
      # sav_external_twcom = save_for_cluster_run(seu_mark, NAME),

      # Pseudo time
      pseudo_time_cl26 = run_slingshot(seu_mark, meta_col = "cluster", clusters = c("2", "6"), reduction = "UMAP"),
      fig_lineage_cl26 = plot_slingshot_lineage(pseudo_time_cl26),
      pseudo_time_beta = run_slingshot(seu_mark, meta_col = "cell_type", clusters = c("MPC", "Beta")),
      fig_lineage_beta = plot_slingshot_lineage(pseudo_time_beta),
      pseudo_time_delta = run_slingshot(seu_mark, meta_col = "cell_type", clusters = c("MPC", "Delta")),
      fig_lineage_delta = plot_slingshot_lineage(pseudo_time_delta),
      seu_mark_pt = seu_mark |>
        add_pseudotime(pseudo_time_beta, colname = "pseudotime_beta") |> 
        add_pseudotime(pseudo_time_delta, colname = "pseudotime_delta") |> 
        add_pseudotime(pseudo_time_cl26, colname = "pseudotime_cl26"),
      png_spatial_pt = cowplot::plot_grid(
        plot_spatial_meta(seu_mark_pt, meta_column = "pseudotime_beta"),
        plot_spatial_meta(seu_mark_pt, meta_column = "pseudotime_delta"),
        ncol = 1, labels = c("Beta", "Delta")
      ) |> gs("spatial_pt", NAME, 8, 12),

      # Differential expression late/MPC
      de_beta_mpc = seurat_de(seu_mark, group_1 = "Beta"),
      de_delta_mpc = seurat_de(seu_mark, group_1 = "Delta"),
      de_mpc = bind_rows(de_beta_mpc, de_delta_mpc),
      fig_volcano_mpc = plot_volcano(de_mpc, fdr_limit = 0.01),

      gse_mpc = fgsea_all_terms(de_mpc, fterms, feature_var = "gene_symbol", rank_expr = "avg_log2FC") |> tibble::add_column(experiment = NAME),

      top_genes_de_beta_mpc = de_beta_mpc |> dplyr::filter(p_val_adj < 0.01 & avg_log2FC > 0) |> dplyr::arrange(p_val) |> head(12) |> pull(gene_symbol),
      png_expression_pt_beta_mpc = plot_expression_pseudotime(seu_mark_pt, top_genes_de_beta_mpc, pt_column = "pseudotime_beta") |> gs("expr_pt_beta", NAME, 8, 7),

      # Save data for external use
      sav_loupe = save_for_loupe(seu_mark_pt, paste0("tab/clusters_", NAME, ".csv")),
      sav_cors = save_correlation(tbl_key_corr_cscore, paste0("tab/correlation_", NAME, ".csv")),
      sav_enr = save_enrichment(tbl_keycorr_enrichment, paste0("tab/keycorr_enrichment_", NAME, ".csv"), group_var = "group"),
      sav_mods = save_spatial_modules(tbl_modules, paste0("tab/spatial_modules_", NAME, ".csv")),
      sav_mod_enr = save_enrichment(tbl_module_enrichment, paste0("tab/module_enrichment_", NAME, ".csv"), group_var = "module"),
      sav_cellchat = save_cellchat(cellchat_comp, NAME),
      sav_cellchat_prob = readr::write_csv(cellchat_prob, paste0("tab/cellchat_prob_", NAME, ".csv")),

      # Prepare for shiny
      dex = make_dexdash(seu_mark, de_mpc, NAME),

      # Integration
      seu_integ = add_clusters_integ(seu_mark, integ_clust, NAME),
      cluster_markers_integ = find_group_markers(seu_integ, group_by = "cluster_integ"),
      png_clust_integ = mn_plot_spatial_clustering(seu_integ, meta_column = "cluster_integ", swap_xy = FALSE) |> gs("clustering_integ", sample = NAME, width = 7, height = 6),
      fig_clust_integ_jaccard = mn_plot_cluster_cell_type_jaccard(seu_integ, cluster_col = "cluster_integ"),
      top_cluster_markers_integ = get_top_cluster_markers(cluster_markers_integ, n_top = 30),
      top_cluster_genes_integ = top_cluster_markers_integ$gene_symbol |> unique(),
      tbl_cluster_enrichment_integ = group_enrichment(top_cluster_markers_integ, fterms, group_var = "cluster") |> dplyr::rename(cluster = group),
      png_dots_top_cluster_integ = plot_gene_dots(seu_integ, top_cluster_genes_integ, meta_column = "cluster_integ") |> gs("dot_clusters_integ", sample = NAME, width = 3.5, height = 29),
      sav_loupe_integ = save_for_loupe(seu_integ, paste0("tab/clusters_integ_", NAME, ".csv"), cluster_col = "cluster_integ"),
  

      # Segmentation data
      seu_seg = load_visium_data(SEG_DATA_PATH, mito_genes, bin_size = "polygons") |>
        NormalizeData() |> 
        cluster_seurat(ident_res = RESOLUTION, resolutions = RESOLUTIONS, nfeatures = 2000, res_suffix = "Spatial.Polygons"),
      fig_seg_hormone_overlap = plot_upset_expressed_cells(seu_seg, c("INS", "SST", "GCG", "GHRL", "PPY")),
      fig_seg_hormneu_overlap = plot_upset_expressed_cells(seu_seg, c("INS", "SST", "GCG", "GHRL", "NEUROG3")),
      seu_seg_mark = annotate_cells_page(seu_seg, tissue_markers, n_boot = 1000, p_val_limit = 0.05),

      # For the manuscript
      mn_clust_integ = mn_plot_spatial_clustering(seu_integ, meta_column = "cluster_integ") |> gp("clustering_integ", sample = NAME, width = 3.3, height = 3, path = "manuscript/fig"),
      mn_clust_jaccard = mn_plot_cluster_cell_type_jaccard(seu_integ, cluster_col = "cluster_integ") |> gp("jaccard_integ", sample = NAME, width = 3.5, height = 2.7, path = "manuscript/fig"),
      mn_clust = mn_plot_spatial_clustering(seu_mark, RESOLUTION) |> gp("clustering", sample = NAME, width = 3.3, height = 3, path = "manuscript/fig"),
      mn_all_markers = mn_plot_marker_spatial(seu_mark, tissue_markers, label_size = 2.5) |> gp("all_markers", sample = NAME, width = 8, height = 8, path = "manuscript/fig"),
      mn_endoc_markers = mn_plot_marker_spatial(seu_mark, tissue_markers_endoc, label_size = 2.5) |> gp("endoc_markers", sample = NAME, width = 8, height = 8, path = "manuscript/fig"),
      mn_cell_types = mn_plot_cell_type(seu_mark) |> gp("cell_type", sample = NAME, width = 3.55, height = 3, path = "manuscript/fig"),
      mn_cell_type_dots = mn_plot_cell_type_dots(seu_mark, tissue_markers) |> gp("cell_type_dots", sample = NAME, width = 4.5, height = 7.5, path = "manuscript/fig"),
      mn_hormone_overlap_8um = mn_plot_hormone_overlap(seu_mark, file.path("manuscript/fig", paste0("hormone_overlap_", NAME, "_8um.pdf")), 7, 3),
      mn_insulin_bins_spatial = mn_plot_binned_spatial(seu_mark) |> gp("insulin_bins", sample = NAME, width = 3.8, height = 3, path = "manuscript/fig"),
      mn_insulin_gradient_expression = plot_binned_expression(seu, "INS", c("NEUROG3", "GCG", "SST"), bins = 6) |> gp("insulin_expression", sample = NAME, width = 6, height = 3, path = "manuscript/fig"),
      mn_metadata = write_meta_for_archive(seu_integ, paste0("manuscript/data_submit/metadata_", NAME, ".csv")),
      mn_wgcna_modules_expr_score = mn_plot_wgcna_modules_expr_score(seu_hdwgcna, n_genes = "all") |> gp("wgcna_modules_score", sample = NAME, width = 8, height = 24, path = "manuscript/fig")
    )
  )

  nmf_analysis <- tar_plan(
    seu_list = list("1" = seu_1, "2" = seu_2),
    nmf = gene_nmf_analysis(seu_list, k = 6:12, nmp = 12),
    seu_nmf_1 = add_module_score_ucell(seu_mark_1, nmf$metaprograms.genes),
    seu_nmf_2 = add_module_score_ucell(seu_mark_2, nmf$metaprograms.genes),

    png_nmf_metaprograms_1 = plot_nmf_metaprograms(seu_nmf_1) |> gs("nmf_metaprograms_1", width = 8, height = 18),
    png_nmf_metaprograms_2 = plot_nmf_metaprograms(seu_nmf_2) |> gs("nmf_metaprograms_2", width = 8, height = 18),

    png_nmf_dots_1 = plot_metaprogram_dots(seu_mark_1, nmf) |> gs("nmf_metaprograms_dots_1", width = 8, height = 30),
    png_nmf_dots_2 = plot_metaprogram_dots(seu_mark_2, nmf) |> gs("nmf_metaprograms_dots_2", width = 8, height = 30)
  )

  combine_samples <- list(
    tar_combine(cell_data_all, map_experiments$cell_data, command = bind_rows(!!!.x)),
    tar_combine(gene_detection_all, map_experiments$gene_detection, command = bind_rows(!!!.x)),
    tar_combine(gse_all, map_experiments$gse_mpc, command = bind_rows(!!!.x))
  )

  cluster_integrated_samples <- tar_plan(
    integ_res = c(0.15, 0.2, 0.25, 0.3),
    seu_integ = integrate_seurats(seu_1, seu_2, nfeatures = 3000),
    seu_integ_clust = cluster_seurat(seu_integ, dims = 1:30, ident_res = 0.25, resolutions = integ_res, nfeatures = 2000, pca_approx = TRUE, res_suffix = "integrated", seed = 123) |> 
                      JoinLayers(assay = "RNA"),  # after integrative analysis we still have two layers
    seu_integ_ord = reorder_clusters(seu_integ_clust, clust_order = c(3, 5, 6, 7, 8, 1, 9, 4, 2, 10)),
                      
    integ_clust = extract_clusters_integ(seu_integ_ord),

    cluster_markers_integ = find_group_markers(seu_integ_ord, group_by = "cluster"),
    max_best_clusters_integ = 30,
    tbl_best_cluster_markers_integ = find_best_markers(cluster_markers_integ, max_per_cluster = max_best_clusters_integ),
    tbl_best_cluster_negative_markers_integ = find_best_markers(cluster_markers_integ, max_per_cluster = max_best_clusters_integ, negative = TRUE),
    top_cluster_markers_integ = get_top_cluster_markers(cluster_markers_integ, n_top = 30),
    top_cluster_genes_integ = top_cluster_markers_integ$gene_symbol |> unique(),

    tbl_cluster_enrichment_integ = group_enrichment(top_cluster_markers_integ, fterms, group_var = "cluster") |> dplyr::rename(cluster = group),
    sav_clust_integ_top_markers = top_cluster_markers_integ |> dplyr::mutate(dplyr::across(where(is.numeric), ~signif(., 4))) |> readr::write_csv("tab/cluster_integ_top_markers.csv"),
    sav_clust_integ_enr = save_enrichment(tbl_cluster_enrichment_integ, paste0("tab/cluster_integ_enrichment.csv"), group_var = "cluster"),

    png_clusters_integ = plot_seurat_clusters(seu_integ_ord, resolutions = integ_res) |> gs("integ_clust", width = 8, height = 6),
    png_top_clust_gene_integ = plot_top_cluster_features(seu_integ_ord, cluster_markers_integ, gns)
  )

  for_report <- tar_plan(
    fig_counts = plot_all_counts(cell_data_all, gene_detection_all),
    png_marginals = plot_grid(fig_marginals_1, fig_marginals_2, ncol = 1, labels = 1:2)  |> gs("marginals", "all", 8, 6),
    png_detections = plot_grid(fig_count_det_1, fig_count_det_2, labels = 1:2) |> gs("detections", "all", 8, 3),

    sct_1 = tibble::tibble(x = c(5120, 5230), y = c(-3400, -5400)),
    png_section_ins_sst = plot_spatial_section(seu_1, "INS", "SST", sct_1) |> gs("spatial_section_ins_sst", width = 8, height = 8),
    png_section_ins_gcg = plot_spatial_section(seu_1, "INS", "GCG", sct_1) |> gs("spatial_section_ins_gcg", width = 8, height = 8),
    png_section_ins_ghrl = plot_spatial_section(seu_1, "INS", "GHRL", sct_1) |> gs("spatial_section_ins_ghrl", width = 8, height = 8),

    box_1 = list(x1 = 4300, x2 = 5600, y1 = -5300, y2 = -3600),
    ins_spots = read_sx_sy("info/insulin_spots.csv", seu_1),
    png_spatial_selection_genes = plot_spatial_selection_genes(seu_1, c("SOX9", "NEUROG3", "INS", "SST", "GCG", "GHRL"), box_1, spots = ins_spots) |> gs("spatial_selection_hormones", width = 8, height = 12),
    png_spatial_selection_meta = plot_spatial_selection_meta(seu_mark_1, "cell_type", box_1) |> gs("spatial_selection_cell_type", width = 5, height = 4),

    box_2 = list(x1 = 4900, x2 = 5200, y1 = -4500, y2 = -4150),
    png_segmentation_grid = plot_gene_segmentation(seu_seg_1, c("INS", "SST"), box_2, seu_grid = seu_1) |> gs("segmentation_grid", width = 8, height = 4.5),

    gif_spatial_ins_gcg = animate_spatial_genes(seu_1, c("INS", "GCG"), box_1),
    gif_spatial_ins_sst = animate_spatial_genes(seu_1, c("INS", "SST"), box_1),
    gif_spatial_ins_ghrl = animate_spatial_genes(seu_1, c("INS", "GHRL"), box_1),
    gif_spatial_hormones = animate_spatial_genes(seu_1, c("INS", "GCG", "SST", "GHRL", "PPY"), box_1),
    gif_spatial_neurog3_ins = animate_spatial_genes(seu_1, c("NEUROG3", "INS"), box_1),
    gif_spatial_neurog3_dbh = animate_spatial_genes(seu_1, c("NEUROG3", "DBH"), box_1),
    gif_spatial_full = animate_spatial_genes(seu_1, c("PTF1A", "PDX1", "SOX9", "PAX4", "NEUROD1", "NEUROG3", "INS"), box_1),

    seu_2um = load_visium_data(EXPERIMENTS$DATA_PATH[1], mito_genes, bin_size = 2),
    box_coords = list(x = 5154, y = -4196, size = 28),
    fig_box_ins = plot_spatial_box(seu_1, "INS", x_cen = box_coords$x, y_cen = box_coords$y, size = box_coords$size),
    res_comp = combine_resolution_area(seu_1, seu_2um, c("INS", "SST", "GCG", "GHRL"), x_cen = box_coords$x, y_cen = box_coords$y, size = box_coords$size),
    fig_res_comp = plot_resolution_comparison(res_comp),

    tbl_key_markers_count = get_total_count(seu_1) |> dplyr::filter(gene_symbol %in% key_markers$gene_symbol),

    tbl_cluster_count = integ_clust |> dplyr::count(cluster, sample),
    tbl_cell_type_count = dplyr::bind_rows(tbl_cell_type_count_1 |> tibble::add_column(sample = "1"), tbl_cell_type_count_2 |> tibble::add_column(sample = "2")),
    fig_cluster_integ_cell_count = plot_counted_cell_perc(tbl_cluster_count, "cluster"),
    fig_cell_type_cell_count = plot_counted_cell_perc(tbl_cell_type_count, "cell_type"),

    tbl_best_cluster_markers = merge_best_markers(tbl_best_cluster_markers_1, tbl_best_cluster_markers_2),
    tbl_best_insulin_gradient_markers = merge_best_markers(tbl_best_insulin_gradient_markers_1, tbl_best_insulin_gradient_markers_2),
    tbl_all_markers = merge_all_markers(tbl_best_cluster_markers, tbl_best_insulin_gradient_markers, key_markers, gns),
    sav_all_markers = readr::write_csv(tbl_all_markers, "tab/markers_for_xenium.csv"),
    tbl_best_cluster_negative_markers = merge_best_markers(tbl_best_cluster_negative_markers_1, tbl_best_cluster_negative_markers_2),
    tbl_best_insulin_gradient_negative_markers = merge_best_markers(tbl_best_insulin_gradient_negative_markers_1, tbl_best_insulin_gradient_negative_markers_2),
    tbl_all_negative_markers = merge_all_markers(tbl_best_cluster_negative_markers, tbl_best_insulin_gradient_negative_markers, key_markers, gns),
    sav_all_negative_markers = readr::write_csv(tbl_all_negative_markers, "tab/negative_markers_for_xenium.csv"),
    top_genes = get_highest_expressors(seu_1, n_top = 19),
    bad_xenium_genes = c(top_genes, "LT", "FOS", "TTR", "HBA2", "MGP"),
    final_markers = prepare_top_markers(tbl_all_markers, n_top = 480, wrong_genes = bad_xenium_genes),

    cellchat_ins = cellchat_insulin_gradient(seu_mark_1, bins = 6),
    fig_cellchat_ins = plot_cellchat_insulin_gradient_ranks(cellchat_ins),
    png_cellchat_ins_insulin = plot_cellchat_chord(cellchat_ins$`(2.94,3.95]`, "INSULIN", "insulin_chord_2.94_3.95.png"),

    tbl_clust_clust_jaccard = jaccard_cluster_cluster(cluster_markers_1, cluster_markers_2, n_top = 50),
    fig_clust_clust_jaccard = plot_cluster_match_samples(cluster_markers_1, cluster_markers_2, n_top = 50)
  )

  for_shiny <- tar_plan(
    sav_de = save_for_shiny_de(dex_1, dex_2, gns, fterms, gse_all)
  )

  for_external_run <- tar_plan(
    # lr = get_ligands_receptors(),
    # sav_lr = saveRDS(lr, file = file.path("ext", "cellchat_ligand_receptor.rds"))
  )

  c(
    map_experiments,
    nmf_analysis,
    combine_samples,
    for_report,
    for_shiny,
    for_external_run,
    cluster_integrated_samples
  )

}