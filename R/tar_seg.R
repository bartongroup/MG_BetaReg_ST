targets_seg <- function() {

  map_experiments <- tar_map(
    values = EXPERIMENTS,
    names = NAME,
    
    tar_plan(
      # Read segmented data
      seus = load_visium_data(SEG_DATA_PATH, mito_genes, bin_size = "polygons") |>
        NormalizeData() |> 
        cluster_seurat(ident_res = RESOLUTION, resolutions = RESOLUTIONS, nfeatures = 2000, res_suffix = "Spatial.Polygons") |> 
        annotate_cells_page(tissue_markers, n_boot = 1000, p_val_limit = 0.05),

      # Integration
      seus_integ = add_clusters_integ(seus, integ_seg_clust, NAME),
      seg_cluster_markers_integ = find_group_markers(seus_integ, group_by = "cluster_integ"),

      # Insulin gradient
      ins_binned_enrichment = enrichment_gradient(seus, fterms, "INS", bins = 6),
      top_50_ins_markers = top_n_cluster_markers(ins_binned_enrichment$markers, n_top = 50),

      # CellChat
      cellchat_seg = create_cellchat(seus) |> preprocess_cellchat(),
      cellchat_seg_comp = compute_cellchat(cellchat_seg, n_boot = 100, scale_distance = CELLCHAT_SCALE_DISTANCE, raw_use = TRUE, population_size = FALSE),

      # WGCNA
      seus_hdwgcna = hd_wgcna_analysis(seus, NAME, assay = "Spatial.Polygons", soft_power = NULL, merge_cut_height = WGCNA_MERGE_CUT_HEIGHT, min_module_size = WGCNA_MIN_MODULE_SIZE),

      # SpaceBender
      # This bit needs to be done on a cluster node due to computational requirements (~256 GB memory).
      # First, save data for spacebender (create targets sav_seus_sb), copy the files from spacebender_data
      # to the cluster node, run spacebender on the cluster (script "scripts/run_spacebender.py"),
      # and copy the results back to the local machine. Then, you can create "seus_sb".
      sav_seus_sb = save_data_for_spacebender_from_seurat(seus, NAME, assay = "Spatial.Polygons"),
      #seus_sb = read_spacebender_clean_file(seus, NAME) |> NormalizeData(),

      # Create figures
      mn_cell_type = mn_plot_meta_segmentation(seus, meta_col = "cell_type", box = BOX, palette = cell_type_pal, fill_lab = "Cell state") |> 
        gp("cell_type", NAME, width = 4.5, height = 5, path = "manuscript/fig"),
      mn_cell_type_dots = mn_plot_cell_type_dots(seus, tissue_markers) |>
        gp("cell_type_dots", sample = NAME, width = 4.5, height = 7.5, path = "manuscript/fig"),
      mn_clust_integ = mn_plot_meta_segmentation(seus_integ, meta_col = "cluster_integ", box = BOX, fill_lab = "Cluster") |>
        gp("clustering_integ", sample = NAME, width = 7, height = 6, path = "manuscript/fig"),
      mn_clust_integ_jaccard = mn_plot_cluster_cell_type_jaccard(seus_integ, cluster_col = "cluster_integ") |> 
        gp("clustering_integ_jaccard", sample = NAME, width = 4, height = 4, path = "manuscript/fig"),
      mn_hormone_overlap = mn_plot_hormone_overlap(seus, file.path("manuscript/fig", paste0("hormone_overlap_", NAME, ".pdf")), 8, 3),

      mn_insulin_bins_spatial = mn_plot_binned_segmentation(seus, box = BOX) |> gp("insulin_bins", sample = NAME, width = 3.8, height = 3, path = "manuscript/fig"),
      mn_insulin_gradient_expression = mn_plot_binned_expression(seus, "INS", c("NEUROG3", "GCG", "SST"), bins = 6) |> gp("insulin_expression", sample = NAME, width = 6, height = 3, path = "manuscript/fig"),
      mn_enrichment_gradient = mn_plot_enrichment_gradient(ins_binned_enrichment$enrichment) |> 
        gp("enrichment_gradient", sample = NAME, width = 6.7, height = 10, path = "manuscript/fig"),
      mn_ins_neurog3 = mn_plot_gene_segmentation(seus, genes = c("INS", "NEUROG3"), box = ZOOM, add_grid = TRUE) |> 
        gp("ins_neurog3", sample = NAME, width = 7, height = 4, path = "manuscript/fig"),
      mn_ins_neurog3_dlk = mn_plot_gene_segmentation(seus, genes = c("INS", "DLK1", "NEUROG3"), box = BOX, ncol = 3) |> 
        gp("ins_neurog3_dlk", sample = NAME, width = 8, height = 4, path = "manuscript/fig"),
      mn_gcg_pcsk = mn_plot_gene_segmentation(seus, genes = c("GCG", "PCSK1", "PCSK1N", "PCSK2"), box = ZOOM, ncol = 2, add_grid = TRUE) |> 
        gp("gcg_pcsk", sample = NAME, width = 7, height = 8, path = "manuscript/fig"),


      mn_cellchat_rank = mn_plot_cellchat_rank(cellchat_seg_comp, n_top = 20) |>
        gp("cellchat_rank_top20", sample = NAME, width = 4, height = 4.5, path = "manuscript/fig"),
      mn_cellchat_chords = mn_plot_cellchat_chords(cellchat_seg_comp, NAME, pathways =  c("FLRT", "NOTCH")),

    )
  )

  separate_figures <- tar_plan(
    box1 = EXPERIMENTS$BOX[1][[1]],
    box2 = EXPERIMENTS$BOX[2][[1]],

    mn_cell_markers_1 = mn_plot_marker_segmentation(seus_1, tissue_markers_endoc, box1) |> gp("cell_markers", "1", 5.5, 8, path = "manuscript/fig"),
    mn_cell_markers_2 = mn_plot_marker_segmentation(seus_2, tissue_markers_endoc, box2) |> gp("cell_markers", "2", 6.6, 8, path = "manuscript/fig"),

    mn_mpc_endoc_modules_1 = mn_plot_modules_expr_score(seus_hdwgcna_1, box1, sel_modules = c("turquoise", "pink"), ncores = 6) |> gp("mpc_endoc_modules", "1", 7, 5, path = "manuscript/fig"),
    mn_mpc_endoc_modules_2 = mn_plot_modules_expr_score(seus_hdwgcna_2, box2, sel_modules = c("turquoise", "brown"), ncores = 6) |> gp("mpc_endoc_modules", "2", 7, 5, path = "manuscript/fig"),
  )

  cluster_integrated_samples <- tar_plan(
    seus_integ = integrate_seurats(seus_1, seus_2, nfeatures = 3000),
    seus_integ_clust = cluster_seurat(seus_integ, dims = 1:30, ident_res = 0.3, resolutions = c(0.2, 0.25, 0.3, 0.35, 0.4), nfeatures = 2000, pca_approx = TRUE, res_suffix = "integrated", seed = 123) |> 
                       JoinLayers(assay = "RNA"),  # after integrative analysis we still have two layers
    seus_integ_ord = reorder_clusters(seus_integ_clust, clust_order = c(6, 2, 4, 7, 3, 1, 5, 8, 9)),
    integ_seg_clust = extract_clusters_integ(seus_integ_ord)
   
    #png_clusters_integ = mn_plot_meta_segmentation(seu_integ_ord, resolutions = integ_res) |> gs("integ_clust", width = 8, height = 6),

  )


  c(
    map_experiments,
    cluster_integrated_samples,
    separate_figures
  )

}