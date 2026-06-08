tar_read(seu_integ) -> seu

dims = 1:30
resolutions = c(0.18, 0.2, 0.22, 0.25)
nfeatures = 2000


seu_integ_clust = cluster_seurat(
  seu,
  dims = dims,
  ident_res = 0.2,
  resolutions = resolutions,
  nfeatures = nfeatures,
  pca_approx = TRUE,
  res_suffix = "integrated",
  seed = 123
) |> 
  JoinLayers(assay = "RNA")

res <- "0.25"
res_col <- str_glue("integrated_snn_res.{res}")
levels(seu_integ_clust@meta.data[[res_col]])

seu_integ_clust$cluster = seu_integ_clust@meta.data[[res_col]]
seu_integ_ord = reorder_clusters(seu_integ_clust, clust_order = c(3, 5, 6, 7, 8, 1, 9, 4, 2, 10))
integ_clust = extract_clusters_integ(seu_integ_ord)

tar_read(seu_mark_1) -> seu_mark

seu_integ = add_clusters_integ(seu_mark, integ_clust, "1")
mn_plot_cluster_cell_type_jaccard(seu_integ, cluster_col = "cluster_integ")
