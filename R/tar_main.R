targets_main <- function() {

  setup_experiment <- tar_plan(
    tar_target(config_file, "config.json", format = "file"),
    config = jsonlite::read_json(config_file)
  )

  get_annotations <- tar_plan(
    # isolate config to avoid re-runs in case of other config changes
    config_annot = list(
      species = config$species,
      ensembl = config$ensembl,
      designations = config$designations
    ), 
    mart = new_mart(config_annot),
    chromosomes = biomart_fetch_chromosomes(mart),
    gns =  biomart_fetch_genes(mart) |> filter(chr %in% chromosomes),
    all_features = load_features(config$feature_file),

    terms = get_functional_terms_sym(config_annot$designations),
    fterms = prepare_terms_fenr(terms, all_features, feature_name = "gene_symbol"),

    cell_cycle_markers = get_cell_cycle_markers(config$species, gns),
    mito_genes =  gns |> filter(chr == "MT") |> pull(gene_symbol),

    tar_target(key_marker_file, "info/key_markers.xlsx", format = "file"),
    key_markers = read_markers(key_marker_file),
    key_markers_endoc = key_markers |> filter(cell_type %in% c("MPC", "EP", "Beta", "Alpha", "Delta", "Epsilon", "PP")) |> droplevels(),

    tar_target(t1d_file, "info/T1D Risk Alleles.csv", format = "file"),
    tar_target(t2d_file, "info/T2D Risk Alleles.csv", format = "file"),
    risk_genes = read_risk_alleles(t1d_file, t2d_file),

    tar_target(genes_of_interest_file, "info/genes.json", format = "file"),
    goi = jsonlite::read_json(genes_of_interest_file, simplifyVector = TRUE),

    tissue_markers = get_unique_markers(goi$tissue_markers_positive, goi$tissue_markers_negative),
    tissue_markers_endoc = tissue_markers |> filter(cell_type %in% c("MPC", "EP", "Beta", "Alpha", "Delta", "Epsilon", "PP")) |> droplevels()
  )

  for_manuscript <- tar_plan(
    sav_tissue_markers = tissue_markers |> readr::write_tsv("manuscript/tab/key_markers.txt"),
    mn_hormone_overlap_2um_1 = mn_plot_hormone_overlap_2um("1", "manuscript/fig/hormone_overlap_1_2um.pdf", 7, 3),
    mn_hormone_overlap_2um_2 = mn_plot_hormone_overlap_2um("2", "manuscript/fig/hormone_overlap_2_2um.pdf", 7, 3),
    mn_hormone_overlap_ma_8pwc = mn_plot_hormone_overlap_ma("8", "manuscript/fig/hormone_overlap_ma_8pcw.pdf", 3.7, 3),
    mn_hormone_overlap_ma_11pwc = mn_plot_hormone_overlap_ma("11", "manuscript/fig/hormone_overlap_ma_11pcw.pdf", 7, 3)
  )




  c(
    setup_experiment,
    get_annotations,
    for_manuscript
  )

}