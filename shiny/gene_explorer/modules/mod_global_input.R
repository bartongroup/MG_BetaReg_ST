# Module for global input for all other modules

require(shiny)

# ----- UI definitions -----

mod_global_input_ui <- function(id) {
  ns <- shiny::NS(id)

  sample_info <- info_icon("sample")
  gene_set_info <- info_icon("gene_set")
  gene_info <- info_icon("gene")

  shiny::tagList(
    shiny::selectInput(
      inputId = ns("sample"),
      label = shiny::span("Sample", sample_info),
      choices = NULL
    ),
    shiny::selectInput(
      inputId = ns("gene_set"),
      label = shiny::span("Gene set", gene_set_info),
      choices = NULL
    ),
    shiny::selectInput(
      inputId = ns("gene"),
      label = shiny::span("Gene", gene_info),
      choices = NULL
    ),
    shiny::selectInput(
      inputId = ns("x_var"),
      label = "X-axis variable",
      choices = c("Cluster" = "cluster", "Cell type" = "cell_type")
    )
  )
}

# ----- Server logic -----

mod_global_input_server <- function(id, state) {
  contrast <- NULL

  get_samples <- function() {
    names(DATA)
  }

  get_gene_sets <- function() {
    names(GENES)
  }

  get_genes <- function(gene_set) {
    GENES[[gene_set]]
  }

  server <- function(input, output, session) {

    # Update dummy sample selection
    shiny::observe({
      shiny::updateSelectInput(
        session = session,
        inputId = "sample",
        choices = get_samples()
      )
    })

    # Update dummy gene set selection
    shiny::observe({
      shiny::updateSelectInput(
        session = session,
        inputId = "gene_set",
        choices = get_gene_sets()
      )
    })

    # update dummy gene selection
    shiny::observe({
      shiny::updateSelectInput(
        session = session,
        inputId = "gene",
        choices = get_genes(names(GENES)[1])
      )
    })

    # Observe sample, copy to state
    shiny::observeEvent(input$sample, {
      sample <- input$sample
      shiny::req(sample)
      state$sample <- sample
    })

    # Observe gene set, copy to state, update genes
    shiny::observeEvent(input$gene_set, {
      gene_set <- input$gene_set
      shiny::req(gene_set)
      state$gene_set <- gene_set
      genes <- get_genes(gene_set)
      shiny::updateSelectInput(
        session = session,
        inputId = "gene",
        choices = genes
      )
    })

    # Observe gene, selection, copy to state
    shiny::observeEvent(input$gene, state$gene_input <- input$gene)

    # Observe x-variable, copy to state
    shiny::observeEvent(input$x_var, state$x_var <- input$x_var)
  }

  shiny::moduleServer(id, server)
}
