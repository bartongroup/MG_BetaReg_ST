# GENE EXPLORER
#
# Main UI
#

require(bslib)
require(shiny)

version <- "1.0"

bslib::page_sidebar(
  theme = bslib::bs_theme(bootswatch = "spacelab", spacer = "0.8rem") |>
    bslib::bs_add_rules(".info-pop { max-width: 500px; background-color: #f0f9e8; }"),
  title = "Beta cell regeneration project",

  sidebar = bslib::sidebar(
    title = "Gene explorer",
    mod_global_input_ui("global_input")
    #shiny::tags$hr(),
    #shiny::tags$span(style = "font-size: 0.7em; color: #191970", version)
  ),

  bslib::layout_columns(
    col_widths = c(3, 5, 4),
    mod_dot_plot_ui("dot_plot"),
    mod_expression_plot_ui("expression_plot"),
    bslib::layout_column_wrap(
      width = 1,
      mod_correlation_table_ui("correlation_table"),
      mod_enrichment_ui("enrichment")
    )
  )
)
