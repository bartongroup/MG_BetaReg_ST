#' Experimental GO Evidence Codes
#'
#' A vector of GO evidence codes that are considered experimental.
#' These codes are used to filter GO terms based on experimental evidence.
#' See https://geneontology.org/docs/guide-go-evidence-codes/
#' 
#' @format A character vector containing the following GO evidence codes:
#' \describe{
#'   \item{EXP}{Inferred from Experiment}
#'   \item{IDA}{Inferred from Direct Assay}
#'   \item{IPI}{Inferred from Physical Interaction}
#'   \item{IMP}{Inferred from Mutant Phenotype}
#'   \item{IGI}{Inferred from Genetic Interaction}
#'   \item{IEP}{Inferred from Expression Pattern}
#'   \item{HTP}{Inferred from High Throughput Experiment}
#'   \item{HDA}{Inferred from High Throughput Direct Assay}
#'   \item{HMP}{Inferred from High Throughput Mutant Phenotype}
#'   \item{HGI}{Inferred from High Throughput Genetic Interaction}
#'   \item{HEP}{Inferred from High Throughput Expression Pattern}
#' }
EXPERIMENTAL_GO_CODES <- c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP") 

#' Get Functional Terms using Gene Symbols
#'
#' This function retrieves functional terms (GO, Reactome, and KEGG) for a given set of gene symbols and IDs.
#' It filters GO terms to only include those with experimental evidence codes.
#'
#' @param designations A list containing species designations for GO, Reactome, and KEGG.
#' @param codes A character vector of GO evidence codes to filter by (default is EXPERIMENTAL_GO_CODES).
#'
#' @return A list containing GO, Reactome, and KEGG terms and mappings.
get_functional_terms_sym <- function(designations, codes = EXPERIMENTAL_GO_CODES) {
  
  cat("Loading GO term data\n")
  go <- fenr::fetch_go(species = designations$go, use_cache = TRUE)
  cat("Loading Reactome data\n")
  re <- fenr::fetch_reactome(species = designations$re)
  cat("Loading KEGG data\n")
  kg <- fenr::fetch_kegg(species = designations$kg)

  go_cc <- go
  go_bp <- go
  go_mf <- go

  go_cc$terms <- go_cc$terms |> 
    filter(namespace == "cellular_component")
  go_bp$terms <- go_bp$terms |> 
    filter(namespace == "biological_process")
  go_mf$terms <- go_mf$terms |> 
    filter(namespace == "molecular_function")

  # Filter non-existent terms and only experimental codes
  go_cc$mapping <- go_cc$mapping |> 
    drop_na() |>
    filter(term_id %in% go_cc$terms$term_id & evidence %in% codes)
  go_bp$mapping <- go_bp$mapping |> 
    drop_na() |>
    filter(term_id %in% go_bp$terms$term_id & evidence %in% codes)
  go_mf$mapping <- go_mf$mapping |> 
    drop_na() |>
    filter(term_id %in% go_mf$terms$term_id & evidence %in% codes)

  re$mapping <- re$mapping

  kg$mapping <- kg$mapping

  list(
    go_cc = go_cc,
    go_bp = go_bp,
    go_mf = go_mf,
    re = re,
    kg = kg
  )
}


#' Prepare Terms for Functional Enrichment Analysis
#'
#' This function prepares functional terms for enrichment analysis using the `fenr` package.
#'
#' @param terms A list of functional terms (GO, Reactome, KEGG) to prepare.
#' @param all_features A vector of all features (e.g., gene symbols) to consider in the enrichment analysis.
#' @param feature_name The name of the feature column to use (default is "gene_symbol").
#'
#' @return A list of prepared terms for each ontology.
prepare_terms_fenr <- function(terms, all_features, feature_name = "gene_symbol") {
  ontologies <- names(terms)

  map(ontologies, function(ont) {
    trm <- terms[[ont]]
    fenr::prepare_for_enrichment(trm$terms, trm$mapping, all_features, feature_name)
  }) |>
    set_names(ontologies)
}

