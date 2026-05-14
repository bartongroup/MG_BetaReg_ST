
top_data_path <- "data/spaceranger_output"
seg_data_path <- "data/spaceranger_new"

#RELABEL_CLUSTERS_2 <- c(6, 9, 3, 1, 2, 7, 8, 10, 5, 4)

EXPERIMENTS_ALL <- tibble::tibble(
  SAMPLE = c("1_18191", "2_17961", "3_18171", "4_17966"),
  NAME = c("1", "2", "3", "4"),
  NAME2 = c("2um_1", "2um_2", "2um_3", "2um_4"),
  RESOLUTIONS = list(
    c(0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5),
    c(0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5),
    c(0.01, 0.025, 0.05),
    c(0.01, 0.025, 0.05)
  ),
  RESOLUTION = c(0.4, 0.3, 0.05,  0.05),
  RELABEL_CLUSTERS = list(NULL, c(6, 9, 3, 1, 2, 7, 8, 10, 5, 4), NULL, NULL),
  RESOLUTIONS_UNKNOWN = list(
    c(0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4),
    c(0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4),
    c(0.01, 0.025, 0.05),
    c(0.01, 0.025, 0.05)
  ),
  RESOLUTION_UNKNOWN = c(0.4, 0.2, 0.05,  0.05),
  N_VARIABLE_FEATURES = c(3000, 3000, 3000, 3000),
  MIN_DETECTED_GENES_PER_CELL = c(10, 10, 10, 10),
  MIN_COUNTS_PER_CELL = c(10, 10, 10, 10),
  MIN_DETECTIONS_PER_GENE = c(5, 5, 5, 5),
  MAX_MT_FRACTION = c(0.1, 0.1, 0.1, 0.1),
  CELLCHAT_SCALE_DISTANCE = c(0.15, 0.15, 0.1, 0.1),
  CELLCHAT_PATHWAYS = list(
    c("NOTCH", "FLRT", "MK", "THBS", "COLLAGEN", "LAMININ"),
    c("NOTCH", "FLRT", "CADM"),
    c("NOTCH"),
    c("NOTCH")
  )
) |> 
  dplyr::mutate(
    DATA_PATH = file.path(top_data_path, stringr::str_glue("{SAMPLE}_Results"), "outs"),
    SEG_DATA_PATH = file.path(seg_data_path, stringr::str_glue("{SAMPLE}_Results"), "outs")
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

read_risk_alleles <- function(t1d_file, t2d_file) {
  t1d <- read_csv(t1d_file)
  t2d <- read_csv(t2d_file)
  list(
    t1d = t1d$Name,
    t2d = t2d$Name
  )
}


tb <- function(df, name, sample = NULL) {
  if(is.null(sample)) {
    fname <- file.path("doc", "tables", paste0(name, ".txt"))
  } else {
    fname <- file.path("doc", "tables", paste0(name, "_", sample, ".txt"))
  }

  readr::write_tsv(df, fname)
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


# At the moment, negative markers are not used
get_unique_markers <- function(tissue_markers_positive, tissue_markers_negative) {
  pos <- tissue_markers_positive |> 
    enframe(name = "cell_type", value = "gene_symbol") |> 
    unnest(gene_symbol) |> 
    add_column(effect = "positive")
  neg <- tissue_markers_negative |> 
    enframe(name = "cell_type", value = "gene_symbol") |> 
    unnest(gene_symbol) |> 
    add_column(effect = "negative")
  # bind_rows(pos, neg)
  pos
}


mark_cell_types <- function(seu, tissue_markers, beta_limit = 0.5) {
  cell_types <- tissue_markers$cell_type |> unique()
  cell_order <- seu@meta.data |> 
    rownames()

  mex <- get_group_gene_expression(seu, tissue_markers, group_var = "cell_type")

  d <- mex |> 
    mutate(
      cell_type = case_when(
        EP <= 0.05 & MPC > 0 & Beta < beta_limit ~ "MPC",
        EP > 0.05 & Beta < beta_limit ~ "EP",
        Alpha > 0.05 & Beta < beta_limit ~ "Alpha",
        Delta > 0.2 & Beta < beta_limit ~ "Delta",
        Epsilon > 0.2 & Beta < beta_limit ~ "Epsilon",
        Beta >= beta_limit ~ "Beta",
        Matrix > 0.2 ~ "Matrix",
        Immune > 0.2 ~ "Immune",
        Nervous > 0.2 ~ "Nervous",
        VSM > 0.2 ~ "VSM",
        .default = "Other"
      )
    ) |> 
    arrange(factor(cell_id, levels = cell_order))

  seu@meta.data$cell_type <- d$cell_type |> 
    factor(levels = c(cell_types, "Other"))

  return(seu)
}


mark_cell_types_by_max <- function(seu, tissue_markers) {
  cell_types <- tissue_markers$cell_type |> unique()
  cell_order <- seu@meta.data |> 
    rownames()

  mex <- get_group_gene_expression(seu, tissue_markers, group_var = "cell_type")

  d <- mex |> 
    pivot_longer(-cell_id) |> 
    group_by(cell_id) |> 
    summarise(cell_type = if_else(sum(value) == 0, "Other", name[which.max(value)])) |> 
    arrange(factor(cell_id, levels = cell_order))

  seu@meta.data$cell_type <- d$cell_type |> 
    factor(levels = c(cell_types, "Other"))

  return(seu)
}


annotate_cell_types <- function(seu, tissue_markers, seed = 42) {
  set.seed(seed)
  cell_types <- tissue_markers$cell_type |> unique()
  tm <- tissue_markers |> 
    group_by(cell_type) |> 
    summarise(genes = list(gene_symbol)) |> 
    deframe()
  s <- AddModuleScore(
    seu,
    features = tm,
    name = "celltype"
  )
  meta_scores <- s@meta.data |> 
    as_tibble(rownames = "cell_id") |> 
    select(cell_id, starts_with("celltype")) |> 
    set_names(c("cell_id"), names(tm))
  ct <- meta_scores  |> 
    pivot_longer(-cell_id) |> 
    mutate(value = if_else(value < 0, 0, value)) |> 
    group_by(cell_id) |> 
    summarise(cell_type = if_else(sum(value) == 0, "Other", name[which.max(value)]))

  meta_extra <- left_join(meta_scores, ct, by = "cell_id") |> 
    select(-cell_id) |> 
    mutate(cell_type = factor(cell_type, levels = c(cell_types, "Other")))

  seu@meta.data$cell_type <- NULL
  seu@meta.data <- bind_cols(seu@meta.data, meta_extra)
  seu
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

# Some tools do not fit in memory of my laptop. Hence,
# I need to save the needed data and run a standalone script on the cluster.
save_for_cluster_run <- function(seu, suffix) {
  saveRDS(seu, file = file.path("ext", str_glue("seu_{suffix}.rds")))
}


# Read correlation data from a parquet file (contains all correlations)
# and select only marker genes in the first column
read_correlation_markers <- function(fname, key_marker_genes_endoc) {
  stopifnot(file.exists(fname))
  tb <- arrow::read_parquet(fname)

  # 'tb' is a long pivot of a triangualar matrix.
  # We need all marker genes as gene_1,
  # but they can be both gene_1 and gene_2. Hence, we split and join.

  tb1 <- tb |> 
    filter(gene_1 %in% key_marker_genes_endoc & !(gene_2 %in% key_marker_genes_endoc)) |> 
    collect()
  tb2 <- tb |> 
    filter(gene_2 %in% key_marker_genes_endoc & !(gene_1 %in% key_marker_genes_endoc)) |> 
    rename(g1 = gene_1, g2 = gene_2) |> 
    rename(gene_1 = g2, gene_2 = g1) |> 
    relocate(gene_1, .before = gene_2) |> 
    collect()
  tb3 <- tb |> 
    filter(gene_1 %in% key_marker_genes_endoc & gene_2 %in% key_marker_genes_endoc) |> 
    collect()
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

save_correlation <- function(cors, fname, fdr_limit = 0.05) {
  cors |> 
    filter(p_adj < fdr_limit) |> 
    mutate(across(where(is.numeric), ~signif(.x, 4))) |> 
    write_csv(fname)
}

save_enrichment <- function(enr, fname, fdr_limit = 0.05, group_var = "group") {
  enr |> 
    select({{ group_var }}, ontology, term_id, term_name, N_with, n_with_sel, odds_ratio, p_value, p_adjust, genes = ids) |> 
    arrange(.data[[group_var]], desc(odds_ratio)) |> 
    mutate(across(where(is.numeric), ~signif(.x, 4))) |> 
    write_csv(fname)
}

save_spatial_modules <- function(mods, fname, fdr_limit = 0.05) {
  mods |> 
    select(gene_symbol = gene_name, module, starts_with("kME")) |> 
    mutate(across(where(is.numeric), ~signif(.x, 4))) |> 
    write_csv(fname)
}

save_cellchat <- function(cellchat, suffix) {
  fname <- file.path("shiny", "data", str_glue("cellchat_{suffix}.qs"))
  qs2::qs_save(cellchat, fname)
}

save_for_shiny <- function(seu, suffix) {
  fname <- file.path("shiny", "data", str_glue("seu_{suffix}.qs"))
  qs2::qs_save(seu, fname)
}


jaccard_clusters_cell_types <- function(seu, cluster_col = "cluster") {
  meta <- seu@meta.data |>
    as_tibble(rownames = "cell_id") |> 
    mutate(cluster = get(cluster_col)) |> 
    select(cell_id, cluster, cell_type)
  m_cl <- meta |>
    group_by(cluster) |> 
    summarize(ids = list(cell_id), .groups = "drop") |> 
    deframe()    
  m_ct <- meta |>
    group_by(cell_type) |> 
    summarize(ids = list(cell_id), .groups = "drop") |> 
    deframe()    

  expand_grid(cluster = levels(meta$cluster), cell_type = levels(meta$cell_type)) |> 
    mutate(
      cluster = factor(cluster, levels = levels(meta$cluster)),
      cell_type = factor(cell_type, levels = levels(meta$cell_type))
    ) |> 
    rowwise() |> 
    mutate(
      n_un = union(m_cl[[cluster]], m_ct[[cell_type]]) |> length(),
      n_int = intersect(m_cl[[cluster]], m_ct[[cell_type]]) |> length(),
      j = n_int / n_un
    )
}

jaccard_cluster_cluster <- function(markers_1, markers_2, n_top = 30) {
  top_1 <- get_top_cluster_markers(markers_1, n_top = n_top)
  top_2 <- get_top_cluster_markers(markers_2, n_top = n_top)
  t1 <- top_1 |>
    group_by(cluster) |> 
    summarize(ids = list(gene_symbol), .groups = "drop") |> 
    deframe()
  t2 <- top_2 |>
    group_by(cluster) |> 
    summarize(ids = list(gene_symbol), .groups = "drop") |> 
    deframe()    

  expand_grid(cluster_1 = levels(top_1$cluster), cluster_2 = levels(top_2$cluster)) |> 
    rowwise() |> 
    mutate(
      n_un = union(t1[[cluster_1]], t2[[cluster_2]]) |> length(),
      n_int = intersect(t1[[cluster_1]], t2[[cluster_2]]) |> length(),
      j = n_int / n_un
    )
}


# Enrichment score from fold change and a vector of gene names
ES <- function(fc, genes, mu, delta) {
  n <- length(genes)
  if(n == 0)
    return(0)
  Sm <- colMeans(fc[genes, , drop = FALSE], na.rm = TRUE)  # mean fold change for the cell type
  ES <- sqrt(n) * (Sm - mu) / delta    # Enrichment score
}

annotate_cells_page <- function(seu, tissue_markers, n_boot = 1000, p_val_limit = 0.05, n_cores = 10, verbose = FALSE) {
  # Normalised data matrix
  m <- GetAssayData(seu, layer = "data")
  # Mean per gene across spots
  gm <- rowMeans(m, na.rm = TRUE)
  keep <- gm > 0
  m <- m[keep, , drop = FALSE]
  gm <- gm[keep]
  good_genes <- names(gm)
  # Fold changes
  fc <- m - gm
  # Mean and std per spot across genes
  mu <- colMeans(fc, na.rm = TRUE)
  delta <- MatrixGenerics::colSds(as.matrix(fc), na.rm = TRUE)

  options(future.globals.maxSize = 100 * 1024^3)
  future::plan(future::multicore, workers = n_cores)

  annot <- tissue_markers |> 
    group_split(cell_type) |> 
    furrr::future_map(function(w) {
      cell_type <- w$cell_type[1]
      if(verbose)
        cat(cell_type, "\n")
      # Enrichment score for positive genes
      genes_pos <- w |> filter(effect == "positive") |> pull(gene_symbol)
      stopifnot(all(genes_pos %in% good_genes))
      ES_obs <- ES(fc, genes_pos, mu, delta)

      # Null distribution
      null_mat <- replicate(n_boot, {
        boot_genes_pos <- sample(good_genes, length(genes_pos))
        ES(fc, boot_genes_pos, mu, delta)
      })

      p_val <- map_dbl(1:length(ES_obs), function(i) {
        (1 + sum(null_mat[i, ] >= ES_obs[i])) / (1 + n_boot)
      })

      tibble(
        cell_type = cell_type,
        cell_id = names(ES_obs),
        ES = ES_obs,
        p_val = p_val
      )
    },
    .options = furrr::furrr_options(seed = TRUE)
  ) |> 
      list_rbind()
  
  future::plan(future::sequential)
  
  cell_types <- tissue_markers$cell_type |> unique()

  score_tab <- annot |> 
    pivot_wider(id_cols = cell_id, names_from = cell_type, values_from = c(ES, p_val))

  cell_annots <- annot |> 
    filter(p_val < p_val_limit) |> 
    group_by(cell_id) |> 
    summarise(cell_type = cell_type[which.max(ES)])

  meta_extra <- tibble(
    cell_id = rownames(seu@meta.data )
  ) |> 
    left_join(score_tab, by = "cell_id") |> 
    left_join(cell_annots, by = "cell_id") |> 
    mutate(cell_type = replace_na(cell_type, "Unknown")) |> 
    select(-cell_id) |> 
    mutate(cell_type = factor(cell_type, levels = c(cell_types, "Unknown")))

  seu@meta.data$cell_type <- NULL
  seu@meta.data <- bind_cols(seu@meta.data, meta_extra)
  seu
}

read_sx_sy <- function(fname, seu) {
  xy <- get_spatial_coordinates(seu) |> 
    mutate(s = str_remove(cell_id, "-1")) |>
    separate_wider_delim(s, delim = "_", names = c("ss", "um", "ix", "iy")) |>
    mutate(sx = as.numeric(ix), sy = as.numeric(iy)) |> 
    select(sx, sy, x, y)
 read_csv(fname, show_col_types = FALSE) |> 
  left_join(xy, by = c("sx", "sy"))
}


gene_binned_gradient <- function(seu, gene, bins, include_zeroes = TRUE) {
  gx <- get_gene_expression(seu, gene)
  g1 <- gx |> 
    filter(expression > 0) |> 
    mutate(bin = cut_interval(expression, n = bins) |> fct_expand("0", after = 0))
  g0 <- gx |> 
    filter(expression == 0) |> 
    mutate(bin = factor("0", levels = levels(g1$bin)))
  
  if(include_zeroes) {
    gx <- bind_rows(g0, g1)
  } else {
    gx <- g1
  }

  # reorder to match metadata
  tibble(cell_id = colnames(seu)) |> 
    left_join(gx, by = "cell_id") |> 
    select(cell_id, bin)
}


enrichment_gradient <- function(seu, fterms, gene = "INS", bins = 6, fdr_limit = 0.01, logfc_limit = 1, pct_limit = 0.1) {
  all_genes <- rownames(seu)
  
  gx <- gene_binned_gradient(seu, gene, bins, include_zeroes = TRUE)
  seu$bin <- gx$bin

  mrk <- find_group_markers(seu, group_by = "bin")

  enr <- mrk |> 
    filter(p_val_adj < fdr_limit & avg_log2FC > logfc_limit & pct.1 > pct_limit) |> 
    group_split(cluster) |> 
    map(function(w) {
      onts <- names(fterms)
      map(onts, function(ont) {
        enr <- fenr::functional_enrichment(all_genes, w$gene_symbol, fterms[[ont]])
        if(!is.null(enr))
          enr <- enr |> filter(p_adjust < 0.01)
        enr
      }) |> 
        bind_rows() |> 
        add_column(bin = w$cluster[1])
    }) |> 
      list_rbind()
  list(
    markers = mrk,
    enrichment = enr
  )
}


# Get top max_per_cluster markers (by logFC rank) from each cluster,
# merge results. Ignore cluster 0. If negative = TRUE, get negative markers.
# For each gene find
#   - best_rank - smallest rank across all clusters
#   - best_cluster - where the smallest rank is
#   - best_pcs - percentage expressed in best cluster
find_best_markers <- function(markers, pct_limit = 0.1, fdr_limit = 0.01, max_per_cluster = 200, negative = FALSE) {
  if(negative) {
    pct_limit <- 0
    markers <- markers |> 
      mutate(avg_log2FC = -avg_log2FC)
  }
  markers |>
    filter(cluster != "0" & p_val_adj < 0.01 & pct.1 > pct_limit) |> 
    group_by(cluster) |> 
    mutate(rank = min_rank(-avg_log2FC)) |> 
    ungroup() |> 
    filter(rank <= max_per_cluster) |>
    group_by(gene_symbol) |> 
    summarise(
      best_rank = min(rank),
      best_cluster = cluster[which.min(rank)],
      best_pct = pct.1[which.min(rank)],
      best_pct_out = pct.2[which.min(rank)]
    ) |>
    arrange(best_rank)
}

merge_best_markers <- function(best_markers_1, best_markers_2) {
  m <- best_markers_1 |>
    full_join(best_markers_2, by = "gene_symbol", suffix = c(".1", ".2"))
  m |> 
    select(gene_symbol, starts_with("best_rank")) |> 
    pivot_longer(-gene_symbol) |>
    group_by(gene_symbol) |>
    mutate(mean_rank = mean(value, na.rm = TRUE)) |> 
    ungroup() |>
    pivot_wider(id_cols = c(gene_symbol, mean_rank)) |> 
    arrange(mean_rank) |> 
    left_join(m |> select(-starts_with("best_rank")), by = join_by(gene_symbol))
}


merge_all_markers <- function(tbl_best_cluster_markers, tbl_best_insulin_gradient_markers, key_markers, gns) {
  km <- key_markers |> 
    distinct(gene_symbol) |> 
    add_column(mean_rank = 0, source = "Key marker")

  mrk <- bind_rows(
    tbl_best_cluster_markers |> add_column(source = "Clustering"),
    tbl_best_insulin_gradient_markers |> add_column(source = "Insulin gradient")
  ) |> 
    filter(!(gene_symbol %in% km$gene_symbol)) |> # don't duplicated key markers
    # deduplicate cluster and insulin, remove lower ranks
    group_by(gene_symbol) |> 
    slice_min(mean_rank, n = 1, with_ties = FALSE) |>
    ungroup()

  g <- gns |> 
    select(gene_symbol, description) |> 
    distinct()
  bind_rows(km, mrk) |> 
    arrange(mean_rank) |> 
    left_join(g, by = "gene_symbol") |> 
    relocate(description, .after = gene_symbol)
}

prepare_top_markers <- function(tbl_all_markers, n_top = 500, wrong_genes = character(0)) {
  tb <- tbl_all_markers |> 
    filter(!(gene_symbol %in% wrong_genes))
  
  tb1 <- tb |> 
    filter(source == "Key marker") |> 
    add_column(final_rank = 0)

  n_part <- ceiling((n_top - nrow(tb1)) / 2)
  tb2 <- tb |>
    filter(source == "Clustering") |>
    slice_min(order_by = mean_rank, n = n_part, with_ties = FALSE) |> 
    mutate(final_rank = seq_along(mean_rank))
  tb3 <- tb |>
    filter(source == "Insulin gradient") |>
    slice_min(order_by = mean_rank, n = n_part, with_ties = FALSE) |> 
    mutate(final_rank = seq_along(mean_rank))

  bind_rows(tb1, tb2, tb3) |> 
    arrange(final_rank) |>
    head(n_top) |> 
    pull(gene_symbol)
}

# The top 19 is a manual selecttion based on the elbow in the total counts plot.
get_highest_expressors <- function(seu, n_top = 19) {
  get_total_count(seu) |> 
    head(n_top) |> 
    pull(gene_symbol)
}


# Written by ChatGPT
fetch_uniprot_fasta <- function(genes, reviewed = TRUE, canonical_only = TRUE, chunk_size = 100) {
  q_parts <- sprintf("gene_exact:%s", genes)
  filters <- c("organism_id:9606")
  if (reviewed)      filters <- c(filters, "reviewed:true")
  if (canonical_only) filters <- c(filters, "is_isoform:false")  # canonical entry per gene

  # Split into chunks to keep URLs reasonable
  q_parts <- sprintf("gene_exact:%s", genes)
  chunks <- split(q_parts, ceiling(seq_along(q_parts) / chunk_size))

  fas <- character()
  for (ch in chunks) {
    # THIS is the important bit: collapse with OR, don't pass the vector itself
    core <- paste(ch, collapse = " OR ")
    q <- sprintf("(%s) AND %s", core, paste(filters, collapse = " AND "))

    url <- paste0("https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=",
                  URLencode(q, reserved = TRUE))

    resp <- request(url) |> req_perform()
    fas <- c(fas, resp_body_string(resp))
  }

  fas
}


reorder_clusters <- function(seu, clust_order, clust_col_name = "cluster") {
  cl2cl <- tibble(
    old_clust = clust_order |> as.character(),
    new_clust = seq_along(old_clust) |> as.character()
  )
  m <- get_cell_data(seu) |> 
    rename(old_clust = !!clust_col_name) |> 
    left_join(cl2cl, by = join_by(old_clust)) |> 
    select(-old_clust) |> 
    mutate(new_clust = factor(new_clust, levels = cl2cl$new_clust)) |> 
    rename(!!clust_col_name := new_clust) |>
    column_to_rownames("cell_id")
  seu@meta.data <- m
  seu
}

# Extract join clustering from integrated object and split it to samples 1 and 2
extract_clusters_integ <- function(seu, clust_order = NULL) {
  m <- seu@meta.data |> 
    as_tibble(rownames = "cell_id") |> 
    select(cell_id, cluster) |> 
    separate_wider_regex(cell_id, patterns = c(cell_id = ".+", "_", sample = "\\d+$")) 
    
  if(!is.null(clust_order)) {
    cl2cl <- tibble(
      old_clust = clust_order |> as.character(),
      new_clust = seq_along(old_clust) |> as.character()
    )
    m <- m |> 
      left_join(cl2cl, by = join_by(cluster == old_clust)) |> 
      select(-cluster, cluster = new_clust)
  }
  
  m |> 
      mutate(cluster = fct_reorder(cluster, cluster |> as.character() |> as.numeric()))
}


add_clusters_integ <- function(seu, clust_info, smpl) {
  ci <- clust_info |> 
    filter(sample == smpl)
  meta <- seu@meta.data[ci$cell_id, ]
  meta$cluster_integ <- ci$cluster
  seu@meta.data <- meta
  seu
}

# The last point of (x, y) is the first point repeated.
polygon_area <- function(x, y, coord2um = 2.195) {
  n <- length(x)
  dn <- sum(x[1:(n - 1)] * y[2:n]) / coord2um^2
  up <- sum(x[2:n] * y[1:(n - 1)]) / coord2um^2
  abs(dn - up) / 2
}

write_meta_for_archive <- function(seu, fname) {
  seu@meta.data |> 
    as_tibble(rownames = "spot_id") |> 
    select(
      spot_id,
      total_count,
      n_genes,
      cluster = cluster_integ,
      cell_type,
      starts_with("ES_"),
      starts_with("p_val_")
    ) |> 
    write_csv(fname)
}