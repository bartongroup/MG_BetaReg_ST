# Module SPATIAL PLOT
#
# Creates spatial plot for a gene
#
# Input:
#    state$sample
#    state$gene_plot
#
# Uses:
#    data$seu - user data
#

require(shiny)
require(bslib)
require(bsicons)
require(htmltools)

require(Seurat)

require(ggplot2)
require(ggbeeswarm)
require(dplyr)
require(tibble)
require(forcats)

require(cowplot)

# ----- UI definitions -----

mod_expression_plot_ui <- function(id) {
  ns <- shiny::NS(id)

  info <- info_icon("expression_plot")

  layer <- shiny::radioButtons(
    inputId = ns("layer"),
    label = "Data layer",
    choices = c("Normalised" = "data", "Counts" = "counts"),
    selected = "data"
  )
  
  spatial_point_size <- shiny::sliderInput(
    inputId = ns("spatial_point_size"),
    label = "Point size",
    min = 0.1,
    max = 3,
    value = 0.7,
    step = 0.02
  )

  palette <- shiny::selectInput(
    inputId = ns("palette"),
    label = "Colour palette",
    choices = c(
      "Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", 
      "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", 
      "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"
    ),
    selected = "Purples"
  )

  gear <- gear_icon(
    layer,
    spatial_point_size,
    palette
  )

  download <- shiny::uiOutput(ns("download_plot"))

  bslib::card(
    bslib::card_header(
      "Expression plots",
      shiny::span(info, download, gear),
      class = "d-flex justify-content-between"
    ),
    bslib::card_body(
      min_height = 600,
      max_height = 1000,
      shiny::textOutput(
        outputId = ns("gene_info")
      ),
      shiny::plotOutput(
        outputId = ns("expression_plot"),
        width = "100%",
        height = "800px"
      )
    )
  )

}

# ----- Server logic -----

mod_expression_plot_server <- function(id, state, data) {

  server <- function(input, output, session) {

    # Gene information
    output$gene_info <- shiny::renderText({
      gene <- state$gene_plot
      shiny::req(gene)
      tb <- data$gene_info |>
        dplyr::filter(gene_symbol == gene)
      if(nrow(tb) == 1) {
        txt <- paste(tb$gene_symbol, tb$description)
      } else {
        txt <- gene
      }
      txt
    })

    # Function to make plot
    make_plot <- shiny::reactive({
      gene <- state$gene_plot
      sample <- state$sample
      group_var <- state$x_var
      shiny::req(gene, sample, group_var)
      pl1 <- data$seu[[sample]] |> 
        sh_plot_spatial(gene, point_size = input$spatial_point_size, dist_palette = input$palette)
      pl2 <- data$seu[[sample]] |> 
        sh_plot_expression(gene, group_var, layer = input$layer)
      cowplot::plot_grid(pl1, pl2, ncol = 1, rel_heights = c(1, 0.65))
    })

    # Output plot
    output$expression_plot <- shiny::renderPlot({
      make_plot()
    })

    output$download_plot <- shiny::renderUI({
      download_link(id)
    })

    # Download handler
    output$handle_download <- shiny::downloadHandler(
      filename = "expression_plot.pdf",
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

get_expression <- function(seu, gene, layer = "data") {
  SeuratObject::GetAssayData(seu, layer = layer)[gene, , drop = FALSE] |> 
    t() |> 
    tibble::as_tibble(rownames = "cell_id") |> 
    tidyr::pivot_longer(-cell_id, names_to = "gene_symbol", values_to = "expression")
}


sh_plot_spatial <- function(seu, gene, point_size = 0.7, dist_palette = "Purples") {

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
  ex <-  get_expression(seu, gene)
  
  d <- dplyr::inner_join(xy, ex, by = "cell_id") 

  pl <- d |> 
    ggplot(aes(x = x, y = y, colour = expression)) +
    thm +
    theme(legend.position = "bottom") +
    geom_point(size = point_size, shape = 15, stroke = 0.1) +
    scale_colour_distiller(palette = dist_palette, direction = 1) +
    labs(x = NULL, y = NULL, title = NULL)
    
  return(pl)
}



sh_plot_expression <- function(seu, gene, group_var, layer) {
  
  mt <- seu@meta.data |> 
    tibble::as_tibble(rownames = "cell_id") |> 
    dplyr::mutate(grp = get(group_var)) |> 
    dplyr::select(cell_id, grp)

  if(group_var == "cluster") {
    mt <- mt |> 
      dplyr::mutate(
        i = grp |> as.character() |> as.numeric(),
        grp = forcats::fct_reorder(grp, i)
      )
  } else {
    mt <- mt |> 
      dplyr::mutate(
        grp = grp |> forcats::as_factor() |> forcats::fct_relevel(c("MPC", "EP", "Beta")) |> forcats::fct_relevel("Unknown", after = Inf),
        i = grp |> as.integer()
      )
  }
  
  ex <- get_expression(seu, gene, layer) |> 
    dplyr::left_join(mt, by = "cell_id")

  groups <- levels(ex$grp)
  m <- ex |> 
    dplyr::group_by(i, grp) |> 
    dplyr::summarise(
      mn = mean(expression[expression > 0], na.rm = TRUE),
      md = median(expression[expression > 0], na.rm = TRUE),
      p = length(which(expression > 0)) / n()
    )
  
  mrg <- 0.2
  max_y <- max(ex$expression)
  min_y <- -mrg * max_y

  m <- m |>
    dplyr::mutate(y = (1 - p) * min_y) |> 
    tidyr::drop_na()
  
  ex <- ex |> 
    dplyr::filter(expression > 0)

  ylab <- ifelse(layer == "data", "Normalised data", "Counts")
  
  q_wid <- 0.3
  m_wid <- 0.3
  ex |> 
    ggplot(aes(x = i, y = expression)) +
    th +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    geom_quasirandom(size = 0.5, colour = "grey60", width = q_wid) +
    geom_segment(data = m, aes(x = i - m_wid, xend = i + m_wid, y = md, yend = md), linewidth = 1) +
    geom_col(data = m, aes(y = min_y), fill = "brown", width = m_wid) +
    geom_col(data = m, aes(y = y), fill = "grey80", width = m_wid) +
    geom_hline(yintercept = 0, colour = "grey70") +
    scale_x_continuous(breaks = seq_along(groups), labels = groups) +
    scale_y_continuous(limits = c(min_y, max_y), expand = expansion(mult = c(0, 0.03))) +
    labs(x = group_var, y = ylab)
}
