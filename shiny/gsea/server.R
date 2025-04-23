# GENE SET ENRICHMENT EXPLORER
#
# Main server
#

server <- function(input, output, session) {
  # Prevents RStudio from crashing when Shiny window closed manually
  session$onSessionEnded(function() {
    stopApp()
  })
  
  # Initialise app state, reactive object for communication between modules
  # Note: an undefined reactive value defaults to NULL, so we can start with empty list
  app_state <- new_app_state(x_variable = "cell_type", colour_variable = "group")
  
  # server logic: modules
  mod_global_input_server("global_input", DATA, app_state)
  mod_gse_table_server("gse_table", DATA, app_state)
  mod_volma_plot_server("volma_plot", DATA, app_state)
  mod_feature_info_server("feature_info", DATA, app_state)
  mod_feature_plot_server("feature_plot", DATA, app_state)
  mod_communication_server("communication", DATA, app_state)
}
