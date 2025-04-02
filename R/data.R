
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
    select(-is)
}


group_enrichment <- function(df, fterms, group_var = "cluster", fdr_limit = 0.05) {
  onts <- names(fterms)
  d <- df |> 
    mutate(group = get(group_var))
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
        EP > 0.15 ~ "EP",
        EP <= 0.15 & MPC > 0 ~ "MPC",
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