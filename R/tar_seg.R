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

      # Create figures
      mn_cell_type = mn_plot_meta_segmentation(seus, meta_col = "cell_type", box = BOX, palette = cell_type_pal, fill_lab = "Cell state") |> 
        gp("seg_cell_type", NAME, width = 4.5, height = 5, path = "manuscript/fig"),
      mn_cell_type_dots = mn_plot_cell_type_dots(seus, tissue_markers) |>
        gp("cell_type_dots", sample = NAME, width = 4.5, height = 7.5, path = "manuscript/fig"),
      mn_clust_integ = mn_plot_meta_segmentation(seus_integ, meta_column = "cluster_integ") |>
        gp("clustering_integ", sample = NAME, width = 7, height = 6),
      mn_clust_integ_jaccard = mn_plot_cluster_cell_type_jaccard(seus_integ, cluster_col = "cluster_integ") |> 
        gp("clustering_integ_jaccard", sample = NAME, width = 6, height = 4.5),
    )
  )

  cluster_integrated_samples <- tar_plan(
    seus_integ = integrate_seurats(seus_1, seus_2, nfeatures = 3000),
    seus_integ_clust = cluster_seurat(seus_integ, dims = 1:30, ident_res = 0.3, resolutions = c(0.2, 0.25, 0.3, 0.35, 0.4), nfeatures = 2000, pca_approx = TRUE, res_suffix = "integrated", seed = 123) |> 
                       JoinLayers(assay = "RNA"),  # after integrative analysis we still have two layers
    seus_integ_ord = reorder_clusters(seus_integ_clust, clust_order = c(6, 2, 4, 7, 3, 1, 5, 8, 9)),
    integ_seg_clust = extract_clusters_integ(seus_integ_ord)
   
    #png_clusters_integ = mn_plot_meta_segmentation(seu_integ_ord, resolutions = integ_res) |> gs("integ_clust", width = 8, height = 6),

  )

  make_figures <- tar_plan(

  )


  c(
    map_experiments,
    cluster_integrated_samples,
    make_figures
  )

}