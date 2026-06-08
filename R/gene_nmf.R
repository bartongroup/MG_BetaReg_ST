gene_nmf_analysis <- function(seu_list, assay = "Spatial.008um", k = 8:12, nmp = 10,
                              weight_explained = 0.7, max_genes = 300, seed = 42) {
  
  seu_list |> 
    GeneNMF::multiNMF(assay = assay, k = k, seed = seed) |> 
    GeneNMF::getMetaPrograms(
      nMP = nmp,
      weight.explained = weight_explained,
      max.genes = max_genes
    )
}

add_module_score_ucell <- function(seu, mp_genes, assay = "Spatial.008um") {
  UCell::AddModuleScore_UCell(seu, features = mp_genes, assay = assay)
}

plot_nmf_metaprograms <- function(seu, ncol = 2) {
  xy <- get_spatial_coordinates(seu, swap_xy = FALSE)
  corners <- get_resolution_grid(xy)
  
  mp_scores <- get_cell_data(seu) |>
    select(cell_id, matches("^MP\\d+_UCell$")) |> 
    pivot_longer(-cell_id, names_to = "program") |> 
    mutate(
      program = str_remove(program, "_UCell"),
      i = str_remove(program, "MP") |> as.integer(),
      program = fct_reorder(program, i)
    ) |> 
    select(-i)

  d <- inner_join(mp_scores, corners, by = "cell_id", relationship = "many-to-many")
  
  d |> 
    ggplot(aes(x = x_corner, y = y_corner, fill = value)) +
    thm +
    geom_polygon(aes(x = x_corner, y = y_corner, group = cell_id, fill = value), colour = NA, linewidth = 0) +
    facet_wrap(~program, ncol = ncol) +
    scale_fill_distiller(palette = "Purples", direction = 1) +
    labs(x = NULL, y = NULL, fill = "Score")
}



plot_metaprogram_dots <- function(seu, nmf) {
  mps <- nmf$metaprograms.genes
  pl <- map(names(mps), function(mp) {
    gns <- mps[[mp]]
    plot_gene_dots(seu, gns, meta_column = "cell_type", with_legend = FALSE, title = mp)
  })
  plot_grid(plotlist = pl, ncol = 2)
}
