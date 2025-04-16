# GENE EXPLORER
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
  app_state <- new_app_state()
  
  # server logic: modules
  mod_global_input_server("global_input", app_state)
  mod_communication_server("communication", app_state)
  mod_dot_plot_server("dot_plot", app_state)
  mod_expression_plot_server("expression_plot", app_state)
  mod_correlation_table_server("correlation_table", app_state)
}
