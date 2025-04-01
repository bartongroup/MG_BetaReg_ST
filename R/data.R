
top_data_path <- "data/spaceranger_output"

EXPERIMENTS <- tibble::tibble(
  SAMPLE = c("1_18191", "2_17961", "3_18171", "4_17966"),
  NAME = c("1", "2", "3", "4"),
  RESOLUTION = c(0.3, 0.3, 0.05,  0.05),
  N_VARIABLE_FEATURES = c(3000, 3000, 3000, 3000),
  MIN_DETECTED_GENES_PER_CELL = c(10, 10, 10, 10),
  MIN_COUNTS_PER_CELL = c(10, 10, 10, 10),
  MIN_DETECTIONS_PER_GENE = c(5, 5, 5, 5),
  MAX_MT_FRACTION = c(0.1, 0.1, 0.1, 0.1)
) |> 
  dplyr::mutate(
    DATA_PATH = file.path(top_data_path, stringr::str_glue("{SAMPLE}_Results"), "outs"),
  )


read_markers <- function(fname) {
  read_excel(fname) |> 
    select(-starts_with("...")) |> 
    pivot_longer(-TF, names_to = "cell_type", values_to = "is") |> 
    drop_na() |> 
    select(-is) |> 
    separate_longer_delim(TF, delim = " / ") |> 
    rename(gene_symbol = TF) |> 
    filter(gene_symbol != "NKX6.1") |>   # redundant gene
    mutate(
      cell_type = factor(cell_type, levels = c("MPC", "EP", "Beta", "Alpha", "Delta", "Epsilon", "PP"))
    )
}


cluster_enrichment <- function(top_cluster_markers, fterms, fdr_limit = 0.05) {
  onts <- names(fterms)
  clusts <- top_cluster_markers$cluster |> unique()
  map(onts, function(ont) {
    trm <- fterms[[ont]]
    map(clusts, function(clust) {
      genes_all <- trm$feature2term |> names()
      genes_sel <- top_cluster_markers |> 
        filter(cluster == clust) |> 
        pull(gene_symbol)
      en <- fenr::functional_enrichment(genes_all, genes_sel, trm)
      if(!is.null(en)) {
        en <- en |> 
          add_column(
            ontology = ont,
            cluster = clust
          )
      }
      en
    }) |> 
      list_rbind()
  }) |> 
    list_rbind() |> 
    filter(p_adjust < fdr_limit)
}


get_unique_markers <- function(mpc_markers, ep_markers, beta_markers) {
  bind_rows(
    tibble(cell_type = "MPC", gene_symbol = mpc_markers),
    tibble(cell_type = "EP", gene_symbol = ep_markers),
    tibble(cell_type = "Beta", gene_symbol = beta_markers)
  )
}


mark_cell_types <- function(seu, select_markers) {
  cell_types <- select_markers$cell_type |> unique()
  cell_order <- seu@meta.data |> 
    rownames()

  d <- map(cell_types, function(ct) {
    gns <- select_markers |> 
      filter(cell_type == ct) |> 
      pull(gene_symbol)
    get_gene_expression(seu, gns, norm_per_gene = TRUE) |> 
      summarise(expression = mean(expression), .by = cell_id) |> 
      add_column(cell_type = ct)
  }) |> 
    list_rbind() |> 
    pivot_wider(id_cols = cell_id, names_from = cell_type, values_from = expression) |> 
    mutate(
      cell_type = case_when(
        Beta > 0.5 ~ "Beta",
        EP > 0.15 ~ "EP",
        EP <= 0.15 & MPC > 0 ~ "MPC",
        .default = "Other"
      )
    ) |> 
    arrange(factor(cell_id, levels = cell_order))

  seu@meta.data$cell_type <- d$cell_type

  return(seu)
}