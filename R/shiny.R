merge_de <- function(de_mpc, de_delta) {
  bind_rows(de_mpc, de_delta)
}

make_dexdash <- function(seu, de, suffix) {
  dex_de <- de |> 
    select(id = gene_symbol, log_fc = avg_log2FC, expr = pct.1, p_value = p_val, fdr = p_val_adj, contrast) |> 
    mutate(p_value = if_else(p_value == 0, .Machine$double.xmin, p_value)) # need to be > 0  for volcano plot
  dex_data <- get_average_expression(seu, group_var = "cell_type") |> 
    select(id = gene_symbol, sample = cell_type, value = data)
  dex_meta <- dex_data |> 
    select(sample) |> 
    distinct() |> 
    mutate(cell_type = sample) |> 
    mutate(
      group = case_match(
        cell_type,
        c("MPC", "EP") ~ "Progenitor",
        c("Alpha", "Beta", "Delta", "Epsilon", "PP") ~ "Endocrine",
        c("Nervous", "VSM", "Immune", "Matrix", "Blood") ~ "Other",
        "Other" ~ "Unmarked"
      )
    ) |> 
    mutate(
      cell_type = factor(cell_type, levels = c("MPC", "EP", "Alpha", "Beta", "Delta", "Epsilon", "PP", "Nervous", "VSM", "Immune", "Matrix", "Blood", "Other"))
    )
  dexdash::dexdash_set(dex_de, dex_data, dex_meta, name = suffix)
}



write_shiny_file <- function(obj, method = c("qs", "rds")) {
  method <- match.arg(method)
  path <- file.path("shiny", "data")
  if (!dir.exists(path)) dir.create(path)
  obj_name <- deparse(substitute(obj))
  cat(paste("  Writing", obj_name))
  if(method == "qs") {
    file_name <- file.path(path, str_glue("{obj_name}.qs"))
    qs2::qs_save(obj, file_name)
  } else if(method == "rds") {
    file_name <- file.path(path, str_glue("{obj_name}.rds"))
    readr::write_rds(obj, file_name, compress = "xz")
  }
  cat("\n")
}


save_for_shiny_de <- function(dex_1, dex_2, gns, fterms, gse) {
  dexset <- dexdash::dexdash_list(dex_1, dex_2)
  
  features <- gns |>
    select(id = gene_symbol, description) |>
    mutate(name = id) |>
    distinct()

  write_shiny_file(dexset)
  write_shiny_file(features)
  write_shiny_file(fterms)
  write_shiny_file(gse)
}