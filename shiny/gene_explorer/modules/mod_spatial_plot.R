# Module SPATIAL PLOT
#
# Creates spatial plot for a gene
#
# Input:
#    state$sample
#    state$gene
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

# ----- UI definitions -----

mod_spatial_plot_ui <- function(id) {
  ns <- shiny::NS(id)

  info <- info_icon("spatial_plot")

  download <- shiny::uiOutput(ns("download_plot"))

  bslib::card(
    bslib::card_header(
      "Spatial plot",
      shiny::span(info, download),
      class = "d-flex justify-content-between"
    ),
    bslib::card_body(
      min_height = 250,
      max_height = 500,
      shiny::plotOutput(
        outputId = ns("spatial_plot"),
        width = "100%",
        height = "300px"
      )
    )
  )

}

# ----- Server logic -----

mod_spatial_plot_server <- function(id, data_set, state) {

  server <- function(input, output, session) {

    # Function to make plot
    make_plot <- shiny::reactive({
      gene <- state$gene
      sample <- state$sample
      shiny::req(gene, sample)
      DATA[[sample]] |> 
        sh_plot_spatial(gene)
    })

    # Output plot
    output$spatial_plot <- shiny::renderPlot({
      make_plot()
    })

    output$download_plot <- shiny::renderUI({
      download_link(id)
    })

    # Download handler
    output$handle_download <- shiny::downloadHandler(
      filename = "spatial_plot.pdf",
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

thm <- th +
  ggplot2::theme(
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  )




sh_plot_spatial <- function(seu, genes, point_size = 0.7) {

  # Create x-y coordinates for spatial
  xy <- SeuratObject::GetTissueCoordinates(seu) |> 
    tibble::as_tibble(rownames = "cell_id") |>
    dplyr::select(-cell) |> 
    # Rearrange cooridnates to look like Newcastle report
    dplyr::rename(
      t = x,
      x = y
    ) |> 
    dplyr::rename(y = t) |> 
    dplyr::mutate(y = -y)

  # Get gene expression per cell
  ex <-  SeuratObject::GetAssayData(seu, layer = "data")[genes, , drop = FALSE] |> 
    t() |> 
    tibble::as_tibble(rownames = "cell_id") |> 
    tidyr::pivot_longer(-cell_id, names_to = "gene_symbol", values_to = "expression")
  
  d <- inner_join(xy, ex, by = "cell_id") 

  pl <- d |> 
    ggplot(aes(x = x, y = y, colour = expression)) +
      thm +
      theme(legend.position = "bottom") +
      geom_point(size = point_size, shape = 15, stroke = 0.1) +
      scale_colour_distiller(palette = "Purples", direction = 1) +
      labs(x = NULL, y = NULL)
    
  
  return(pl)
}

