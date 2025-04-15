# Module DOT PLOT
#
# Creates dot plot for a selection of genes
#
# Input:
#    state$sample
#    state$gene_set
#
# Uses:
#    data_set - user data
#

require(shiny)
require(bslib)
require(bsicons)
require(htmltools)

require(Seurat)

require(ggplot2)
require(dplyr)
require(tibble)
require(forcats)

# ----- UI definitions -----

mod_dot_plot_ui <- function(id) {
  ns <- shiny::NS(id)

  x_variable <- shiny::selectInput(
    inputId = ns("x_var"),
    label = "X-axis variable",
    choices = c("cluster", "cell_type")
  )

  info <- info_icon("dot_plot")

  download <- shiny::uiOutput(ns("download_plot"))

  gear <- gear_icon(
    x_variable
  )

  bslib::card(
    bslib::card_header(
      "Dot plot",
      shiny::span(info, download, gear),
      class = "d-flex justify-content-between"
    ),
    bslib::card_body(
      min_height = 250,
      shiny::plotOutput(
        outputId = ns("dot_plot"),
        width = "100%",
        height = "300px"
      )
    )
  )

}

# ----- Server logic -----

mod_dot_plot_server <- function(id, data_set, state) {

  server <- function(input, output, session) {

    # Function to make plot
    make_plot <- shiny::reactive({
      gene_set <- state$gene_set
      sample <- state$sample
      meta_col <- input$x_var
      shiny::req(gene_set, sample, meta_col)
      genes <- GENES[[gene_set]]
      DATA[[sample]] |> 
        sh_plot_gene_dots(genes, meta_column = meta_col)
    })

    # Output plot
    output$dot_plot <- shiny::renderPlot({
      make_plot()
    })

    output$download_plot <- shiny::renderUI({
      download_link(id)
    })

    # Download handler
    output$handle_download <- shiny::downloadHandler(
      filename = "dot_plot.pdf",
      content = function(file) {
        ggplot2::ggsave(file, plot = make_plot(), device = "pdf", width = 6, height = 6)
      }
    )
  }

  shiny::moduleServer(id, server)
}


# ----- Functions used in the module -----

th <- ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank()
  )


sh_plot_gene_dots <- function(seu, gns, meta_column = "cluster", label_size = 10) {
  SeuratObject::Idents(seu) <- meta_column
  
  dp <- Seurat::DotPlot(seu, features = gns)
  d <- dp$data |> 
    tibble::as_tibble(rownames = "gene_symbol")

  if(meta_column == "cluster") {
    d <- d |> 
      mutate(
        i = id |> as.character() |> as.numeric(),
        id = forcats::fct_reorder(id, i)
      )
  }
   
  ggplot(d, aes(x = id, y = features.plot, colour = avg.exp.scaled, size = pct.exp)) +
    th +
    theme(
      axis.text.y = element_text(size = label_size),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    geom_point() +
    scale_colour_gradient(low = "grey90", high = "blue") +
    scale_radius(range = c(1, 6)) +
    labs(x = meta_column, y = "Gene", size = "Percent expressed", colour = "Average expresion")
}

