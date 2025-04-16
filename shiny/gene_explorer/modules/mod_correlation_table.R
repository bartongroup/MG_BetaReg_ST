# Correlation table
# Creates a table with genes correlated with the selected gene

require(shiny)
require(bslib)
require(bsicons)

require(DT)


# ----- UI definitions -----

mod_correlation_table_ui <- function(id) {
  ns <- NS(id)

  fdr_limit <- shiny::numericInput(
    inputId = ns("fdr_limit"),
    label = "FDR limit",
    value = 0.05,
    min = 0,
    max = 1
  )

  correlation_limit <- shiny::numericInput(
    inputId = ns("correlation_limit"),
    label = "Correlation limit",
    value = 0.5,
    min = -1,
    max = 1
  )

  gear <- gear_icon(
    fdr_limit,
    correlation_limit
  )

  info <- info_icon("correlation_table")

  download <- shiny::uiOutput(ns("download_table"))

  bslib::card(
    bslib::card_header(
      "Correlated genes",
      shiny::span(info, download, gear),
      class = "d-flex justify-content-between"
    ),
    DT::dataTableOutput(
      outputId = ns("correlation_table")
    )
  )
  
}

# ----- Server logic -----

mod_correlation_table_server <- function(id, state) {
  
  server <- function(input, output, session) {
    
    # Prepare updated feature info table
    make_table <- reactive({
      sample <- state$sample
      gene <- state$gene_input
      req(sample, gene)
      cr <- COR[[sample]] |> dplyr::filter(p_adj < input$fdr_limit & correlation > input$correlation_limit)

      # For space saving the correlation matrix is triangular, so the gene of interest could be
      # in either column
      cr1 <- cr |>
        dplyr::filter(gene_1 == gene) |>
        dplyr::select(gene_symbol = gene_2, correlation) |>
        dplyr::collect()
      cr2 <- cr |>
        dplyr::filter(gene_2 == gene) |>
        dplyr::select(gene_symbol = gene_1, correlation) |>
        dplyr::collect()
      prc <- percentage_expressed(DATA[[sample]])

      dplyr::bind_rows(cr1, cr2) |> 
        dplyr::filter(gene_symbol != gene) |> 
        dplyr::arrange(desc(correlation)) |> 
        dplyr::left_join(prc, by = "gene_symbol")
    })
    
    # Wait for row selection in the table, pass it to app state
    observeEvent(input$correlation_table_rows_selected, ignoreNULL = FALSE, {
      rows_sel <- input$correlation_table_rows_selected
      if(!is.null(rows_sel)) {
        ct <- make_table()
        state$gene_cor <- ct[rows_sel, ]$gene_symbol  
      } else {
        state$gene_cor <- NULL
      }
    })
    
    # feature_info table
    output$correlation_table <- DT::renderDataTable({
      ct <- make_table() |> 
        dplyr::mutate(
          correlation = signif(correlation, digits = 2),
          perc = signif(perc, digits = 2)
        ) |> 
        dplyr::left_join(GENE_INFO, by = "gene_symbol") |> 
        dplyr::relocate(correlation, .after = description) |> 
        dplyr::relocate(perc, .after = correlation) |> 
        dplyr::rename(Gene = gene_symbol, Description = description, r = correlation, `%` = perc)
      DT::datatable(
        ct,
        options = list(paging = FALSE),
        style = "bootstrap",
        selection = "single",
        rownames = FALSE
      ) |> 
        DT::formatStyle(columns = colnames(ct), fontSize = '80%')
    })

    # Display download icon only when there are data to show
    output$download_table <- shiny::renderUI({
      shiny::req(state$sample, state$gene)
      download_link(id)
    })

    # Download handler
    output$handle_download <- shiny::downloadHandler(
      filename = "correlation_table.csv",
      content = function(file) {
        make_table() |>
          readr::write_csv(file)
      }
    )

  }
  moduleServer(id, server)
}
