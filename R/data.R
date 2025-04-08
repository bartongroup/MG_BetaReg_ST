
top_data_path <- "data/spaceranger_output"

EXPERIMENTS_ALL <- tibble::tibble(
  SAMPLE = c("1_18191", "2_17961", "3_18171", "4_17966"),
  NAME = c("1", "2", "3", "4"),
  RESOLUTIONS = list(
    c(0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5),
    c(0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5),
    c(0.01, 0.025, 0.05),
    c(0.01, 0.025, 0.05)
  ),
  RESOLUTION = c(0.4, 0.3, 0.05,  0.05),
  N_VARIABLE_FEATURES = c(3000, 3000, 3000, 3000),
  MIN_DETECTED_GENES_PER_CELL = c(10, 10, 10, 10),
  MIN_COUNTS_PER_CELL = c(10, 10, 10, 10),
  MIN_DETECTIONS_PER_GENE = c(5, 5, 5, 5),
  MAX_MT_FRACTION = c(0.1, 0.1, 0.1, 0.1)
) |> 
  dplyr::mutate(
    DATA_PATH = file.path(top_data_path, stringr::str_glue("{SAMPLE}_Results"), "outs"),
  )

# Only use samples 1 and 2. Samples 3 and 4 are failed.
EXPERIMENTS <- EXPERIMENTS_ALL[c(1,2), ]

read_markers <- function(fname) {
  read_excel(fname) |> 
    pivot_longer(-gene_symbol, names_to = "cell_type", values_to = "is") |> 
    mutate(cell_type = as_factor(cell_type)) |> 
    drop_na() |> 
    select(-is) |> 
    arrange(cell_type) |> 
    mutate(gene_symbol = as_factor(gene_symbol))
}


group_enrichment <- function(df, fterms, group_var = "cluster", gene_var = "gene_symbol", fdr_limit = 0.05) {
  onts <- names(fterms)
  d <- df |> 
    mutate(
      group = get(group_var),
      gene_symbol = get(gene_var)
    )
  groups <- d$group |> as.character() |> unique()
  res <- map(onts, function(ont) {
    trm <- fterms[[ont]]
    map(groups, function(grp) {
      genes_all <- trm$feature2term |> names()
      genes_sel <- d |> 
        filter(group == grp) |> 
        pull(gene_symbol)
      en <- NULL
      if(length(genes_sel) > 1)
        en <- fenr::functional_enrichment(genes_all, genes_sel, trm)
      if(!is.null(en)) {
        en <- en |> 
          add_column(
            ontology = ont,
            group = grp
          )
      }
      en
    }) |> 
      list_rbind()
  }) |> 
    list_rbind()
  if(!is.null(res) & nrow(res) > 0)
    res <- res |> filter(p_adjust < fdr_limit)

  return(res)
}


get_unique_markers <- function(mpc_markers, ep_markers, beta_markers, nervous_markers,
                               vsm_markers, immune_markers, matrix_markers) {
  bind_rows(
    tibble(cell_type = "MPC", gene_symbol = mpc_markers),
    tibble(cell_type = "EP", gene_symbol = ep_markers),
    tibble(cell_type = "Beta", gene_symbol = beta_markers),
    tibble(cell_type = "Nervous", gene_symbol = nervous_markers),
    tibble(cell_type = "VSM", gene_symbol = vsm_markers),
    tibble(cell_type = "Immune", gene_symbol = immune_markers),
    tibble(cell_type = "Matrix", gene_symbol = matrix_markers)
  )
}


mark_cell_types <- function(seu, select_markers) {
  cell_types <- select_markers$cell_type |> unique()
  cell_order <- seu@meta.data |> 
    rownames()

  mex <- get_group_gene_expression(seu, select_markers, group_var = "cell_type")

  d <- mex |> 
    mutate(
      cell_type = case_when(
        Beta > 0.4 ~ "Beta",
        EP > 0.05 ~ "EP",
        EP <= 0.05 & MPC > 0 ~ "MPC",
        Matrix > 0.2 ~ "Matrix",
        Immune > 0.2 ~ "Immune",
        Nervous > 0.2 ~ "Nervous",
        .default = "Other"
      )
    ) |> 
    arrange(factor(cell_id, levels = cell_order))

  seu@meta.data$cell_type <- d$cell_type

  return(seu)
}


find_top_cor_genes <- function(tbl_cor, fdr_limit = 0.01, min_cor = NULL, n_top = NULL) {
  tb <- tbl_cor |> 
    filter(p_adj < fdr_limit)

  if(!is.null(min_cor)) {
    good_genes <- tb |> 
      filter(correlation >= min_cor) |> 
      pull(gene_2) |> 
      unique()      
  } else if(!is.null(n_top)) {
    good_genes <- tb |> 
      group_by(gene_1) |> 
      slice_max(order_by = correlation, n = n_top) |> 
      pull(gene_2) |> 
      unique()
  } else {
    good_genes <- unique(tb$gene_2)
  }
  good_genes
}

tabulate_correlated_genes <- function(tbl_cor, fdr_limit = 0.01, n_top = 16) {
  tbl_cor |>
    filter(p_adj < fdr_limit) |> 
    group_by(gene_1) |> 
    arrange(desc(correlation)) |> 
    mutate(r = row_number()) |> 
    filter(r <= n_top) |> 
    summarise(genes = gene_2 |> str_c(collapse = ", ")) |> 
    rename(key_gene = gene_1)
}

# CS-SCORE correlation does not fit in memory of my laptop. Hence,
# I need to save the needed data and run a standalone script on the cluster.
save_for_external_correlation <- function(seu, key_marker_genes_endoc, suffix) {
  saveRDS(seu, file = file.path("tab", str_glue("seu_{suffix}.rds")))
  saveRDS(key_marker_genes_endoc, file = file.path("tab", str_glue("key_marker_genes_endoc_{suffix}.rds")))
}

read_correlation_table <- function(fname, key_marker_genes_endoc) {
  stopifnot(file.exists(fname))
  tb <- qs_read(fname)

  # 'tb' is a long pivot of a triangualar matrix. Here, we want to filter 
  # 
  # 
  # We need all marker genes as gene_1,
  # but they can be both gene_1 and gene_2. Hence, we split and join.

  tb1 <- tb |> 
    filter(gene_1 %in% key_marker_genes_endoc & !(gene_2 %in% key_marker_genes_endoc))
  tb2 <- tb |> 
    filter(gene_2 %in% key_marker_genes_endoc & !(gene_1 %in% key_marker_genes_endoc)) |> 
    rename(g1 = gene_1, g2 = gene_2) |> 
    rename(gene_1 = g2, gene_2 = g1) |> 
    relocate(gene_1, .before = gene_2)
  tb3 <- tb |> 
    filter(gene_1 %in% key_marker_genes_endoc & gene_2 %in% key_marker_genes_endoc)
  bind_rows(tb1, tb2, tb3) |> 
    mutate(gene_1 = factor(gene_1, levels = key_marker_genes_endoc))
}


symmetrise_table <- function(cors) {
  cors |> 
    filter(gene_1 != gene_2) |> 
    mutate(g1 = gene_2, g2 = gene_1) |> 
    select(-c(gene_1, gene_2)) |> 
    rename(gene_1 = g1, gene_2 = g2) |> 
    bind_rows(cors)
}