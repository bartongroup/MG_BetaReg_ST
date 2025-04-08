
wgcna_power <- function(tab, R2_cut = 0.85) {
  # very hard to stop the bloody function from printing stuff
  invisible({capture.output({
    sft <- WGCNA::pickSoftThreshold.fromSimilarity(abs(tab), verbose = 0, RsquaredCut = R2_cut)
  })})
  sft
}

adjacency_to_pairs <- function(adj) {
  adj |> 
    as_tibble(rownames = "node_1") |>
    pivot_longer(-node_1, names_to = "node_2", values_to = "adjacency") |> 
    filter(node_1 < node_2)
}

merge_modules <- function(tab, mods, clust_method, me_diss_cut) {
  merged <- WGCNA::mergeCloseModules(tab, as.integer(mods), cutHeight = me_diss_cut)
  mods_merged <- merged$colors
  
  tibble(
    gene_symbol := colnames(tab),
    module = as.integer(mods),
    merged_module = mods_merged,
    colour = WGCNA::labels2colors(module),
    merged_colour = WGCNA::labels2colors(merged_module)
  ) |> 
    mutate(across(c(colour, merged_colour), \(x) x |> as_factor() |> fct_relevel("grey")))
}

wgcna_analysis <- function(cors_file, 
                           gene_sel = NULL, ncores = 4, R2_cut = 0.85, power = NULL,
                           adjacency_type = "unsigned", clust_method = "complete",
                           min_module_size = 30, me_diss_cut = 0.25, deep_split = 2) {
  
  cors <- qs_read(cors_file)
  
  if(!is.null(gene_sel)) {
    cors <- cors |> 
      filter(gene_1 %in% gene_sel & gene_2 %in% gene_sel)
  }

  # Make correlation matrix
  cat("Symmetrising table\n")
  cors <- symmetrise_table(cors)

  cat("Creating correlation matrix\n")
  tab <- cors |> 
    pivot_wider(id_cols = gene_1, names_from = gene_2, values_from = correlation) |> 
    column_to_rownames("gene_1") |> 
    as.matrix()
  cat(dim(tab))
  cat("\n")
  genes <- colnames(tab)

  # Make it symmetric
  tab <- tab[genes, ]

  # Threads break power calculation in WGCN (at least on my laptop)
  # WGCNA::disableWGCNAThreads()

  WGCNA::enableWGCNAThreads(nThreads = ncores)

  # thresholding power
  cat("Power\n")
  sft <- wgcna_power(tab, R2_cut)

  if(is.null(power)) {
    pwr <- sft$powerEstimate
    if(is.na(pwr))
      stop(str_glue("Unable to find power with R2 cutoff {R2_cut}. Try lowering R2_cut."))
    cat(paste("Best power = ", pwr, "\n"))
  } else {
    pwr <- power
    cat(str_glue("Using provided power of {pwr}.\n"))
  }
  sftc <- sft$fitIndices |>
    as_tibble() |>
    janitor::clean_names()

  #adjacency matrix
  cat("Adjacency matrix\n")
  adj <- WGCNA::adjacency.fromSimilarity(abs(tab), power = pwr, type = adjacency_type)
  cat("\n")
  
  # Topological Overlap Matrix (TOM)
  cat("TOM\n")
  tom <- WGCNA::TOMsimilarity(adj)
  
  # clustering
  cat("Clustering\n")
  tom_diss <- 1 - tom
  tree <- hclust(as.dist(tom_diss), method = clust_method)
  
  # module identification - dynamic tree cut
  cat("Modules\n")
  mods <- dynamicTreeCut::cutreeDynamic(dendro = tree, distM = tom_diss, deepSplit = deep_split,
                                        pamRespectsDendro = FALSE, minClusterSize = min_module_size)
  
  cat("Eigengenes\n")
  me <- WGCNA::moduleEigengenes(tab, as.integer(mods))
  me_diss <- 1 - cor(me$eigengenes)
  me_tree <-  hclust(as.dist(me_diss), method = clust_method)

  tb_mods <- tibble(
    gene_symbol = genes,
    module = as.integer(mods)
  )

  list(
    power = pwr,
    sft = sftc,
    min_module_size = min_module_size,
    clust_method = clust_method,
    me_diss_cut = me_diss_cut,
    deep_split = deep_split,
    genes = genes,
    tom = tom,
    tree = tree,
    me_tree = me_tree,
    mods = tb_mods
  )
}


wgcna_create_network <- function(wg, mod, tom_cut = 0.1) {
  sel <- wg$mods |> 
    mutate(i = seq_along(gene_symbol)) |> 
    filter(merged_module == mod) 

  # TOM subset for genes in a given module
  tom <- wg$tom[sel$i, sel$i]
  # create binary adjacency matrix based on TOM cutoff
  tom[tom > tom_cut] <- 1
  tom[tom != 1] <- 0
  
  ntw <- tom |> 
    igraph::graph_from_adjacency_matrix() |> 
    igraph::simplify()
  
  igraph::V(ntw)$name <- sel$gene_symbol
  ntw <- igraph::delete_vertices(ntw, igraph::degree(ntw) == 0)
  
  ntw
}


plot_wgcna_profiles <- function(set, wg) {
  d <- set$dat |>
    filter(sample %in% wg$samples) |> 
    group_by(gene_id) |>
    mutate(
      m = mean(rlog),
      r = rlog - m
    ) |>
    ungroup() |>
    right_join(wg$mods, by = "gene_id") |> 
    filter(merged_module > 0)
  dm <- d |> 
    group_by(merged_module, sample) |> 
    summarise(prof = mean(r))

  d |> 
    ggplot() +
    geom_line(aes(x = sample, y = r, group = gene_id), colour = "grey50", alpha = 0.3) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    geom_line(data = dm, aes(x = sample, y = prof, group = 1), colour = "black") +
    facet_wrap(~ merged_module)
}

plot_wgcna_scale <- function(wg, threshold = 0.9) {
  g1 <- wg$sft |> 
    ggplot(aes(x = power, y = -sign(slope) * truncated_r_sq)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_point() +
    geom_hline(yintercept = threshold, colour = "red") +
    labs(x = "Power", y = "signed R2")
  g2 <- wg$sft |> 
    ggplot(aes(x = power, y = mean_k)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_point() +
    labs(x = "Power", y = "Mean connectivity")
  cowplot::plot_grid(g1, g2, nrow = 1, align= "h")
}

plot_wgcna_heatmap <- function(set, wg, mod) {
  gns <- wg$mods |>
    filter(module == mod) |>
    pull(gene_id)
  plot_de_heatmap(set, gns)
}

plot_wgcna_dendro <- function(wg, pdf_file) {
  pdf(pdf_file)
  WGCNA::plotDendroAndColors(wg$tree, cbind(wg$mods$colour, wg$mods$merged_colour), dendroLabels = FALSE)
  dev.off()
}


tabulate_wgcna <- function(wg, genes) {
  genes <- genes |> 
    select(gene_id, gene_symbol = gene_name, description) |> 
    distinct()
  wg$mods |> 
    filter(merged_module > 0) |> 
    left_join(genes, by = "gene_id") |> 
    select(gene_id, module = merged_module, gene_symbol, description) |> 
    arrange(module)
}


wgcna_enrichment <- function(wg, all_genes, fterms) {
  mods <- wg$mods |> 
    filter(merged_module > 0)
  modules <- unique(mods$merged_module) |> 
    sort()
  ontologies <- names(fterms)
  map_dfr(ontologies, function(ont) {
    map_dfr(modules, function(mod) {
      sel_genes <- mods |> 
        filter(merged_module == mod) |> 
        pull(gene_symbol) |> 
        unique()
      fe <- fenr::functional_enrichment(all_genes, sel_genes, fterms[[ont]])
      if(!is.null(fe))
        fe <- fe |> add_column(module = mod)
    })
  })
}