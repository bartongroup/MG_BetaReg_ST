require(shiny)
require(bslib)
require(bsicons)
require(htmltools)
require(SeuratObject)
require(tibble)
require(dplyr)
require(tidyr)
require(purrr)
require(readxl)
require(rlang)
require(qs2)
require(arrow)

#' Create a Download Link
#'
#' This function creates a Shiny download link with a specified ID.
#'
#' @param id A character string specifying the ID for the download link.
#' @return A Shiny download link UI component.
download_link <- function(id) {
  shiny::downloadLink(
    outputId = shiny::NS(id, "handle_download"),
    label = bsicons::bs_icon("download")
  )
}

#' Create an Info Icon with Popover
#'
#' This function creates an info icon that, when clicked, displays a popover
#' with content from a markdown file.
#'
#' @param helper A character string specifying the name of the markdown file
#'   (without the extension) located in the "inst/helpers" directory.
#' @return A Shiny UI component with an info icon and a popover.
info_icon <- function(helper) {
  bslib::popover(
    bsicons::bs_icon("info-circle"),
    htmltools::includeMarkdown(system.file(stringr::str_glue("helpers/{helper}.md"), package = "dexdash")),
    options = list(customClass = "info-pop")
  )
}

#' Create a Gear Icon with Popover
#'
#' This function creates a gear icon that, when clicked, displays a popover with the specified content.
#'
#' @param ... Additional arguments to be passed to the popover content.
#' @return A Shiny UI component with a gear icon and a popover.
gear_icon <- function(...) {
  bslib::popover(
    bsicons::bs_icon("gear"),
    ...
  )
}

#' Create reactive container for app state
#'
#' Creates a container for app state to exchange information between modules. The
#' initialisation of fields is not necessary, as by default all undefined
#' fields are NULL. However, this seems like a good spot to gather all the field
#' names and comments together.
#'
#' @return A reactive values list with pre-initialised fields.
#' @noRd
new_app_state <- function() {
  shiny::reactiveValues(
    sample = NULL,           # sample name, from the drop-down menu
    gene_set = NULL,         # gene set selection from the drop-down menu
    gene_input = NULL,       # gene selection from the drop-down menu
    gene_cor = NULL,         # gene selected manually from the correlation table
    gene_plot = NULL,        # Gene to plot expression; calculated in communication mod
    sel_enrichment = NULL,   # all genes from the correlation table
    x_var = NULL,            # variable to put on x-axis of dot and expression plots
  )
}

#' Calculate percentage of cells expressing each gene
#'
#' This function calculates the percentage of cells in which each gene is expressed 
#' (count > 0) in a Seurat object.
#'
#' @param seu A Seurat object containing gene expression data.
#'
#' @return A tibble with two columns:
#' \itemize{
#'   \item gene_symbol: Character vector of gene names
#'   \item perc: Numeric vector of percentages (0-100) of cells expressing each gene
#' }
percentage_expressed <- function(seu) {
  ex <- SeuratObject::GetAssayData(seu, layer = "count")
  
  nonzero_cell_counts <- rowSums(ex != 0)
  total_cells <- ncol(ex)

  tibble::tibble(
    gene_symbol = names(nonzero_cell_counts),
    perc = (nonzero_cell_counts / total_cells) * 100
  )
}

#' Read gene sets from Excel file and match with available genes
#'
#' This function reads gene sets organized in Excel sheets and matches them against
#' a list of available genes. Each sheet in the Excel file represents a different
#' gene set or pathway.
#'
#' @param excel_file Path to the Excel file containing gene sets (one per sheet)
#' @param all_genes Character vector of all available genes to match against
#'
#' @return A named list where each element contains the genes found in both the 
#' Excel sheet and the all_genes vector. Names correspond to Excel sheet names.
read_genes_for_pathways <- function(excel_file, all_genes) {
  groups <- readxl::excel_sheets(excel_file)
  purrr::map(groups, function(group) {
    tb <- readxl::read_excel(excel_file, sheet = group, na = "#N/A", range = cell_cols("A")) |> 
      dplyr::rename(gene_symbol = 1) |> 
      tidyr::drop_na() |> 
      tidyr::separate_longer_delim(gene_symbol, delim = ",") |> 
      dplyr::mutate(gene_symbol = str_trim(gene_symbol)) |> 
      dplyr::distinct() |> 
      dplyr::mutate(found = gene_symbol %in% all_genes)
    missing_genes <- tb |> 
      dplyr::filter(!found) |> 
      dplyr::pull(gene_symbol)
    if(length(missing_genes) > 0) {
    #  print(paste(group, ": ", stringr::str_c(missing_genes, collapse = ", ")))
    }
    tb |> 
      dplyr::filter(found) |> 
      dplyr::pull(gene_symbol)
  }) |> 
    rlang::set_names(groups)
}



#' Read all data files for analysis
#'
#' This function reads multiple data files required for analysis, including:
#' - Seurat objects (gene expression data)
#' - Correlation tables
#' - Genes of interest from Excel file
#' - Gene descriptions
#' - Functional terms
#'
#' @return A list containing:
#' \itemize{
#'   \item seu: Named list of Seurat objects
#'   \item cor: Named list of correlation tables
#'   \item genes: Gene sets from Excel file
#'   \item gene_info: Gene descriptions
#'   \item fterms: Functional terms
#' }
read_all_data <- function() {
  
  sets <- c("1", "2")

  # Read Seurat objects
  seu <- purrr::map(sets, ~qs2::qs_read(stringr::str_glue("../data/seu_{.x}.qs"))) |> 
    rlang::set_names(sets)

  # Read correlation tables
  cor <- purrr::map(sets, ~arrow::open_dataset(stringr::str_glue("../data/cor_{.x}.parquet"))) |> 
    rlang::set_names(sets)
    
  # Read genes of interest
  all_genes <- SeuratObject::Features(data$seu$`1`)
  genes <- read_genes_for_pathways("genes_of_interest.xlsx", all_genes)

  # Read gene descriptions
  gene_info <- qs2::qs_read("../data/gene_info.qs")

  # Read functional terms
  fterms <- qs2::qs_read("../data/fterms.qs")

  list(
    seu = seu,
    cor = cor,
    genes = genes,
    gene_info = gene_info,
    fterms = fterms
  )
}