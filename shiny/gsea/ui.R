# GENE SET ENRICHMENT EXPLORER
#
# Main UI
#

version <- "1.0"

bslib::page_sidebar(
  theme = bslib::bs_theme(bootswatch = "spacelab", spacer = "0.8rem") |>
    bslib::bs_add_rules(".info-pop { max-width: 500px; background-color: #f0f9e8; }"),
  title = "DEXDASH",

  sidebar = bslib::sidebar(
    title = "GSE explorer",
    mod_global_input_ui("global_input"),
    shiny::tags$hr(),
    shiny::tags$span(style = "font-size: 0.7em; color: #191970", version)
  ),

  bslib::layout_column_wrap(
    width = 1/3,
    bslib::layout_column_wrap(
      width = 1,
      mod_volma_plot_ui("volma_plot"),
      mod_feature_plot_ui("feature_plot")
    ),
    mod_feature_info_ui("feature_info"),
    mod_gse_table_ui("gse_table"),
  )
)
