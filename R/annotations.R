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
#' Fetch UniProt protein sequences for SignalP
#'
#' Fetch canonical/displayed UniProt sequences in batches, verify that every
#' requested gene maps to exactly one entry, and create SignalP-safe FASTA
#' records. Ambiguous, missing, or invalid records are never selected
#' arbitrarily.
#'
#' @param genes Gene symbols to retrieve.
#' @param organism_id NCBI taxonomy identifier. Defaults to human (9606).
#' @param reviewed Restrict results to reviewed UniProtKB/Swiss-Prot entries.
#' @param chunk_size Maximum number of gene symbols in each UniProt request.
#' @param strict If `TRUE`, stop when any gene is missing, ambiguous, or has a
#'   sequence that SignalP cannot accept. If `FALSE`, warn and omit such genes.
#'
#' @return A list containing `fasta`, a single SignalP-ready FASTA string;
#'   `mapping`, one row per selected gene; `missing`, a character vector;
#'   `ambiguous`, all matches for genes with multiple entries; and `invalid`,
#'   entries whose sequences do not satisfy SignalP input requirements.
fetch_uniprot_fasta <- function(
  genes,
  organism_id = 9606L,
  reviewed = TRUE,
  chunk_size = 100L,
  strict = TRUE
) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Package `httr2` is required to query UniProt.", call. = FALSE)
  }

  genes <- unique(trimws(as.character(genes)))
  if (!length(genes) || anyNA(genes) || any(genes == "")) {
    stop("`genes` must contain non-empty, non-NA gene symbols.", call. = FALSE)
  }
  if (
    length(organism_id) != 1L ||
      is.na(organism_id) ||
      !grepl("^[0-9]+$", as.character(organism_id))
  ) {
    stop("`organism_id` must be a single NCBI taxonomy identifier.", call. = FALSE)
  }
  if (!is.logical(reviewed) || length(reviewed) != 1L || is.na(reviewed)) {
    stop("`reviewed` must be TRUE or FALSE.", call. = FALSE)
  }
  if (
    !is.numeric(chunk_size) ||
      length(chunk_size) != 1L ||
      is.na(chunk_size) ||
      !is.finite(chunk_size) ||
      chunk_size < 1L ||
      chunk_size != floor(chunk_size) ||
      chunk_size > .Machine$integer.max
  ) {
    stop("`chunk_size` must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be TRUE or FALSE.", call. = FALSE)
  }

  chunk_size <- as.integer(chunk_size)
  chunks <- split(genes, ceiling(seq_along(genes) / chunk_size))
  fields <- paste(
    c("accession", "id", "gene_primary", "gene_names", "protein_name", "length", "sequence"),
    collapse = ","
  )

  quote_query_value <- function(x) {
    x <- gsub("\\\\", "\\\\\\\\", x, fixed = TRUE)
    x <- gsub("\"", "\\\\\"", x, fixed = TRUE)
    paste0("\"", x, "\"")
  }

  fetch_chunk <- function(chunk) {
    gene_query <- paste0(
      "gene_exact:",
      vapply(chunk, quote_query_value, character(1)),
      collapse = " OR "
    )
    query <- sprintf(
      "(%s) AND organism_id:%s%s",
      gene_query,
      organism_id,
      if (reviewed) " AND reviewed:true" else ""
    )

    response <- httr2::request("https://rest.uniprot.org/uniprotkb/stream") |>
      httr2::req_url_query(query = query, format = "tsv", fields = fields) |>
      httr2::req_retry(max_tries = 3) |>
      httr2::req_timeout(60) |>
      httr2::req_perform()

    body <- httr2::resp_body_string(response)
    tab <- utils::read.delim(
      text = body,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    required_columns <- c(
      "Entry", "Entry Name", "Gene Names", "Protein names", "Length", "Sequence"
    )
    absent_columns <- setdiff(required_columns, names(tab))
    if (length(absent_columns)) {
      stop(
        "Unexpected UniProt response; missing column(s): ",
        paste(absent_columns, collapse = ", "),
        call. = FALSE
      )
    }
    if (!nrow(tab)) {
      return(NULL)
    }

    gene_names <- strsplit(tab[["Gene Names"]], "[[:space:]]+")
    mapped <- lapply(seq_len(nrow(tab)), function(i) {
      requested <- chunk[
        toupper(chunk) %in% toupper(gene_names[[i]])
      ]
      if (!length(requested)) {
        return(NULL)
      }
      data.frame(
        gene_symbol = requested,
        accession = tab[["Entry"]][i],
        entry_name = tab[["Entry Name"]][i],
        protein_name = tab[["Protein names"]][i],
        length = tab[["Length"]][i],
        sequence = tab[["Sequence"]][i],
        stringsAsFactors = FALSE
      )
    })
    mapped <- Filter(Negate(is.null), mapped)
    if (!length(mapped)) NULL else do.call(rbind, mapped)
  }

  matches <- Filter(Negate(is.null), lapply(chunks, fetch_chunk))
  if (length(matches)) {
    matches <- unique(do.call(rbind, matches))
  } else {
    matches <- data.frame(
      gene_symbol = character(),
      accession = character(),
      entry_name = character(),
      protein_name = character(),
      length = integer(),
      sequence = character(),
      stringsAsFactors = FALSE
    )
  }

  hit_count <- table(factor(matches$gene_symbol, levels = genes))
  missing <- genes[hit_count == 0L]
  ambiguous_genes <- genes[hit_count > 1L]
  ambiguous <- matches[matches$gene_symbol %in% ambiguous_genes, , drop = FALSE]
  selected_genes <- genes[hit_count == 1L]
  mapping <- matches[
    match(selected_genes, matches$gene_symbol),
    ,
    drop = FALSE
  ]

  allowed_amino_acids <- "^[ACDEFGHIKLMNPQRSTVWYXUBZO]+$"
  invalid_reason <- ifelse(
    mapping$length < 10L,
    "shorter than 10 amino acids",
    ifelse(
      mapping$length > 10000L,
      "longer than 10,000 amino acids",
      ifelse(
        !grepl(allowed_amino_acids, mapping$sequence),
        "contains characters not accepted by SignalP",
        NA_character_
      )
    )
  )
  invalid <- mapping[!is.na(invalid_reason), , drop = FALSE]
  invalid$reason <- invalid_reason[!is.na(invalid_reason)]
  mapping <- mapping[is.na(invalid_reason), , drop = FALSE]

  summarize_genes <- function(x, n = 10L) {
    shown <- paste(utils::head(x, n), collapse = ", ")
    if (length(x) > n) paste0(shown, sprintf(" (+%d more)", length(x) - n)) else shown
  }

  problems <- character()
  if (length(missing)) {
    problems <- c(problems, paste0("missing: ", summarize_genes(missing)))
  }
  if (length(ambiguous_genes)) {
    problems <- c(problems, paste0("ambiguous: ", summarize_genes(ambiguous_genes)))
  }
  if (nrow(invalid)) {
    problems <- c(problems, paste0("invalid for SignalP: ", summarize_genes(invalid$gene_symbol)))
  }
  if (length(problems)) {
    message <- paste(problems, collapse = "; ")
    if (strict) {
      stop(
        "UniProt retrieval did not produce one valid sequence per gene (",
        message,
        "). Re-run with `strict = FALSE` to return valid records and inspect the report.",
        call. = FALSE
      )
    }
    warning("Omitting UniProt records (", message, ").", call. = FALSE)
  }

  wrap_sequence <- function(sequence, width = 60L) {
    starts <- seq.int(1L, nchar(sequence), by = width)
    paste0(substring(sequence, starts, starts + width - 1L), collapse = "\n")
  }
  safe_header_component <- function(x) {
    gsub("[^A-Za-z0-9_.:-]", "_", x)
  }

  if (nrow(mapping)) {
    headers <- paste(
      safe_header_component(mapping$accession),
      safe_header_component(mapping$entry_name),
      safe_header_component(mapping$gene_symbol),
      sep = "|"
    )
    records <- paste0(
      ">", headers, "\n",
      vapply(mapping$sequence, wrap_sequence, character(1))
    )
    fasta <- paste(records, collapse = "\n")
  } else {
    fasta <- ""
  }

  list(
    fasta = fasta,
    mapping = mapping,
    missing = missing,
    ambiguous = ambiguous,
    invalid = invalid
  )
}
