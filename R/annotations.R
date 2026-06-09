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
  cat("Loading WIKI data\n")
  wi <- fenr::fetch_wiki(species = designations$wi)

  # Older version of fenr (<1.4.2) does not return evidence column. In such case do no filtering.
  if(is.null(go$mapping[["evidence"]])) {
    go$mapping$evidence <- "none"
    codes <- "none"
  }

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

  wi$mapping <- wi$mapping |> 
    filter(type == "GeneProduct") |> 
    select(term_id, gene_symbol = text_label)

  list(
    go_cc = go_cc,
    go_bp = go_bp,
    go_mf = go_mf,
    re = re,
    kg = kg,
    wi = wi
  )
}


fetch_biocarta <- function() {
  lnk <- "https://maayanlab.cloud/static/hdfs/harmonizome/data/biocarta/gene_attribute_edges.txt.gz"

  raw <- read_tsv(lnk, skip = 1, show_col_types = FALSE) |> 
    select(gene_symbol = GeneSym, term_id = `Pathway ID`, term_name = Pathway) |> 
    mutate(term_id = as.character(term_id))

  terms <- raw |> 
    select(term_id, term_name) |> 
    distinct()
  mapping <- raw |> 
    select(gene_symbol, term_id) |> 
    distinct()

  list(
    terms = terms,
    mapping = mapping
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




# Written by ChatGPT
fetch_uniprot_fasta <- function(
  genes,
  organism_id = 9606,
  reviewed = TRUE,
  max_hits_per_gene = 5
) {
  genes <- unique(as.character(genes))
  base_url <- "https://rest.uniprot.org/uniprotkb/search"

  one_gene <- function(gene) {
    query <- sprintf(
      "gene_exact:%s AND organism_id:%s%s",
      gene,
      organism_id,
      if (reviewed) " AND reviewed:true" else ""
    )

    url <- paste0(
      base_url,
      "?query=", utils::URLencode(query, reserved = TRUE),
      "&format=tsv",
      "&fields=accession,id,gene_primary,protein_name,organism_name,length,sequence",
      "&size=", max_hits_per_gene
    )

    tab <- tryCatch(
      utils::read.delim(url, sep = "\t", quote = "", stringsAsFactors = FALSE),
      error = function(e) NULL
    )

    if (is.null(tab) || nrow(tab) == 0) {
      warning(sprintf("No UniProt match found for gene symbol: %s", gene))
      return(NULL)
    }

    if (nrow(tab) > 1) {
      warning(sprintf(
        "Multiple UniProt matches found for %s; using first accession: %s",
        gene, tab$Entry[1]
      ))
    }

    hit <- tab[1, ]

    header <- sprintf(
      ">%s|%s|%s OS=%s GN=%s PE=canonical",
      hit$Entry,
      hit$`Entry.Name`,
      hit$`Protein.names`,
      hit$Organism,
      hit$`Gene.Names..primary.`
    )

    seq_wrapped <- paste(
      substring(
        hit$Sequence,
        seq(1, nchar(hit$Sequence), by = 60),
        pmin(seq(60, nchar(hit$Sequence) + 59, by = 60), nchar(hit$Sequence))
      ),
      collapse = "\n"
    )

    paste(header, seq_wrapped, sep = "\n")
  }

  records <- Filter(Negate(is.null), lapply(genes, one_gene))
  paste(records, collapse = "\n")
}