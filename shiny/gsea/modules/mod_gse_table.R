# GSE table
# Creates a table with GSE results

# ----- UI definitions -----

mod_gse_table_ui <- function(id) {
  ns <- NS(id)

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
    max = 1
  )

  gear <- gear_icon(
    ontology,
    fdr_limit
  )

  info <- info_icon("gse_table")

  download <- shiny::uiOutput(ns("download_table"))

  bslib::card(
    bslib::card_header(
      "GSE table",
      shiny::span(info, download, gear),
      class = "d-flex justify-content-between"
    ),
    DT::dataTableOutput(
      outputId = ns("gse_table")
    )
  )
  
}

# ----- Server logic -----

mod_gse_table_server <- function(id, data_set, state) {
  
  server <- function(input, output, session) {
    

    # Update dummy ontology selections from data
    ontologies <- data_set$gse$ontology |> unique()
    shiny::observe({
      shiny::updateSelectInput(
        session = session,
        inputId = "ontology",
        choices = ontologies
      )
    })

    # Prepare updated feature info table
    make_table <- reactive({
      ctr <- state$contrast
      ont <- input$ontology
      set_name <- state$set_name
      req(ctr, ont)
      data_set$gse |>
        dplyr::filter(ontology == ont & experiment == set_name & contrast == ctr & fdr <= input$fdr_limit) |> 
        dplyr::select(TermId = term_id, Name = term_name, nes) |> 
        dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~signif(.x, 2))) |> 
        dplyr::mutate(TermId = as.character(TermId))
    })
    
    # Wait for row selection in the GSE, pass it to app state
    observeEvent(input$gse_table_rows_selected, ignoreNULL = FALSE, {
      rows_sel <- input$gse_table_rows_selected
      set_name <- state$set_name
      ctr <- state$contrast
      ont <- input$ontology
      req(ctr, ont)
      ids <- NULL
      if(!is.null(rows_sel)) {
        gt <- make_table()
        termid <- gt[rows_sel, ]$TermId
        r <- data_set$gse |> 
          dplyr::filter(ontology == ont & experiment == set_name &  contrast == ctr & term_id == termid)
        if(nrow(r) > 0) {
          ids <- r$leading_edge[[1]]
        }
      }
      state$sel_term <- ids
    })
    
    # feature_info table
    output$gse_table <- DT::renderDataTable({
      gt <- make_table()
      DT::datatable(
        gt,
        options = list(paging = FALSE),
        style = "bootstrap",
        selection = "single",
        rownames = FALSE
      ) |> 
        DT::formatStyle(columns = colnames(gt), fontSize = '80%')
    })

    # Display download icon only when there are data to show
    output$download_table <- shiny::renderUI({
      shiny::req(state$contrast, state$set_name)
      download_link(id)
    })

    # Download handler
    output$handle_download <- shiny::downloadHandler(
      filename = "gse_table.csv",
      content = function(file) {
        make_table() |>
          readr::write_csv(file)
      }
    )


  }


  
  moduleServer(id, server)
}
