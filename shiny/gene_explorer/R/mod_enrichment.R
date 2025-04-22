# Module FUNCTIONAL ENRICHMENT
#
# Creates a table with functional enrichment for a selection of features from the
# correlation table.
#
# Input:
#    state$sel_enrichment - selection of feature IDs to use in functional enrichment.
#
# Uses:
#    data_set$fterms - a list of fterm objects prepared by fterm
#

require(shiny)
require(bslib)
require(bsicons)

require(SeuratObject)
require(DT)
require(readr)
require(dplyr)
require(tibble)
require(fenr)


# ----- UI definitions -----

mod_enrichment_ui <- function(id) {
  ns <- shiny::NS(id)

  ontology <- shiny::selectInput(
    inputId = ns("ontology"),
    label = "Ontology",
    choices = ""
  )

  fdr_limit <- shiny::numericInput(
    inputId = ns("fdr_limit"),
    label = "FDR limit",
    value = 0.05,
    min = 0,
    max = 1,
    step = 0.01
  )

  gear <- gear_icon(
    ontology,
    fdr_limit
  )

  info <- info_icon("enrichment")

  download <- shiny::uiOutput(ns("download_table"))

  bslib::card(
    bslib::card_header(
      "Functional enrichment",
      shiny::span(info, download, gear),
      class = "d-flex justify-content-between"
    ),
    bslib::card_body(
      DT::dataTableOutput(
        outputId = ns("enrichment")
      )
    )
  )

}

# ----- Server logic -----

mod_enrichment_server <- function(id, state, data) {

  server <- function(input, output, session) {

    # Update dummy ontology selections from data
    ontologies <- names(data$fterms)
    shiny::observe({
      shiny::updateSelectInput(
        session = session,
        inputId = "ontology",
        choices = ontologies
      )
    })

    # Create and return a functional enrichment table using fenr.
    make_table <- shiny::reactive({
      sel_ids <- state$sel_enrichment
      sample <- state$sample
      shiny::req(!is.null(sel_ids) & length(sel_ids) > 1, sample)
      all_ids <- SeuratObject::Features(data$seu[[sample]])
      make_functional_enrichment(sel_ids, all_ids, data$fterms[[input$ontology]], input$fdr_limit)
    })

    # Output DT table
    output$enrichment <- DT::renderDataTable({
      fe <- make_table()
      if(!is.null(fe) & !("Error" %in% colnames(fe))) {
        fe <- fe |>
          dplyr::select(TermID, Name, n = n_with_sel, OR, ids)
        DT::datatable(
          fe,
          options = list(paging = FALSE, dom = "t"),
          style = "bootstrap",
          selection = "single",
          rownames = FALSE
        ) |>
          DT::formatStyle(columns = colnames(fe), fontSize = "75%")
      }
    })

    # Display download icon only when there are data to show
    output$download_table <- shiny::renderUI({
      shiny::req(state$sel_enrichment)
      download_link(id)
    })

    # Download handler
    output$handle_download <- shiny::downloadHandler(
      filename = "enrichment.csv",
      content = function(file) {
        make_table() |>
          readr::write_csv(file)
      }
    )
  }

  shiny::moduleServer(id, server)
}


# ----- Functions used in the module -----


#' Make functional enrichment table
#'
#' @param sel_ids Selected feature IDs
#' @param all_ids All feature IDs (background)
#' @param trms An fterm object with functional term data
#' @param fdr_limit FDR limit for selection of terms
#' @param max_points Maximum number of points selectable
#'
#' @return A tibble with significant functional terms
#' @noRd
make_functional_enrichment <- function(sel_ids, all_ids, trms, fdr_limit, max_points = 3000) {
  n <- length(sel_ids)

  fe <- NULL
  if(n > 1 && n <= max_points) {
    fe <- fenr::functional_enrichment(all_ids, sel_ids, trms)
    if(!is.null(fe)) {
      fe <- fe |>
        dplyr::filter(p_adjust < fdr_limit) |>
        dplyr::select(TermID = term_id, Name = term_name, N_with, n_with_sel, OR = odds_ratio, ids, p_value, FDR = p_adjust)
    }
  } else if (n > max_points) {
    fe <- tibble::tibble(Error = stringr::str_glue("Only {max_points} points can be selected."))
  }
  fe
}
